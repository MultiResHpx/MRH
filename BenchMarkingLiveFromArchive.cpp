#include <time.h>
#include <ctime>
#include <vector>
#include <algorithm>
#include <math.h>
#include <limits>
#include <cmath>
#include <cstdlib>
#include "lsconstants.h"
#include <stdlib.h>
#include <string>
#include "MultiResHpx_Map.h"

#include <Windows.h>

#include "geom_utils.h"
#include "healpix_map.h"

using namespace::std;

#define TA1_DEBUG true
#define NEIGHBORS_DEBUG true
#define STRIP_DEBUG true
#define POLY_DEBUG true
#define POLY_ACC_DEBUG true
#define DISC_ACC_DEBUG true
#define RANDGEN1 false
#define EMPTY -99999
#define SANITYCHECK
#define ATTEMPTS 5
#define TWOPI 6.2831853071796
#define MAXHPXORDER 12

#define DO_DISC
#define DO_POLY
#define DO_STRIP
#define DO_NEIGHBORS


int HPX_CRIT_COUNT = 0;
int MRH_CRIT_COUNT = 0;
int MRH_COVER_MAP_SIZE = 0;
int HPX_COVER_MAP_SIZE = 0;
int MRH_UPSEARCH_COUNT = 0;
int MRH_COVER_MAP_CELL_RES = 0;

#ifdef _WIN32
#include <Windows.h>
#else
#include <sys/time.h>
#include <ctime>
#endif

/* Remove if already defined */
typedef long long int64; typedef unsigned long long uint64;

/* Returns the amount of milliseconds elapsed since the UNIX epoch. Works on both
 * windows and linux. */

uint64 GetTimeMs64()
{
#ifdef _WIN32
 /* Windows */
 FILETIME ft;
 LARGE_INTEGER li;

 /* Get the amount of 100 nano seconds intervals elapsed since January 1, 1601 (UTC) and copy it
  * to a LARGE_INTEGER structure. */
 GetSystemTimeAsFileTime(&ft);
 li.LowPart = ft.dwLowDateTime;
 li.HighPart = ft.dwHighDateTime;

 uint64 ret = li.QuadPart;
 ret -= 116444736000000000LL; /* Convert from file time to UNIX epoch time. */
 ret /= 10000; /* From 100 nano seconds (10^-7) to 1 millisecond (10^-3) intervals */

 return ret;
#else
 /* Linux */
 struct timeval tv;

 gettimeofday(&tv, NULL);

 uint64 ret = tv.tv_usec;
 /* Convert from micro seconds (10^-6) to milliseconds (10^-3) */
 ret /= 1000;

 /* Adds the seconds (10^0) after converting them to milliseconds (10^-3) */
 ret += (tv.tv_sec * 1000);

 return ret;
#endif
}

bool AlmostEqual(double a,double b) {
	if( fabs(a-b) < std::numeric_limits<double>::epsilon() ) { 
		return true;
	}
	return false;
	//return std::fabs(a-b) < std::numeric_limits<double>::epsilon()*std::abs(a+b)*ulp  || std::abs(a-b) < std::numeric_limits<double>::min();
}




class DiscType
{
public:
	DiscType(){pt.phi=0.0;pt.theta=0.0;radius=0.0;};
	DiscType(pointing p,double r){pt=p;radius=r;};
	~DiscType(){};
	pointing pt;
	double radius;
};

class PolyType
{
public:
	PolyType(){pts.clear();};
	PolyType(std::vector<pointing> p){pts=p;};
	~PolyType(){};
	std::vector<pointing> pts;
};

class StripType
{
public:
	StripType(){theta1 = 0.0; theta2 = 0.0;};
	StripType(double t1, double t2){theta1=t1;theta2=t2;};
	~StripType(){};
	double theta1;
	double theta2;
};

class NeighborType
{
public:
	NeighborType(){pt.phi = 0.0; pt.theta = 0.0;};
	NeighborType(pointing _pt){pt = _pt;};
	~NeighborType(){};
	pointing pt;
};

class Measurement 
{
public:
	Measurement(){rec = EMPTY;};

	~Measurement(){};

	bool equals(const Measurement &other);
	
	// Data
	int64 rec; //unique record number of measurement
	pointing pt; //spatial location of measurement
	
	int	data1; 
	int data2;
	int data3;
	int data4;

	double val1;
	double val2;
	double val3;
	double val4;
};

inline bool Measurement::equals(const Measurement &other){
	if ( this->rec == other.rec ) { return true; }
	return false;
}


std::vector<Measurement> ReadRandomMeasurements(ifstream& fp,int64& hpxOrder)
{
	std::vector<Measurement> measurements;
	Measurement m;
	int64 numRecs;

	measurements.clear();
             
	// First line is number of records and HPX Order of dataset (maximum resolution or minimum
	// spatial distance between data points.
	fp >> numRecs >> hpxOrder;  
	skip_line(fp);

	// Skip header line
	skip_line(fp);

	// Next lines are Measurement records
	for(int i = 0; i < numRecs; i++)
	{
		fp >> m.rec >> m.pt.phi >> m.pt.theta 
			>> m.data1 >> m.data2 >> m.data3 >> m.data4
			>> m.val1 >> m.val2 >> m.val3 >> m.val4;
	 
		measurements.push_back(m);

		// Skip the rest of the line
		skip_line(fp);

		if(i%1000 == 0) {
			cout << "Read " << i << " records...\n";
		}

	}
	return( measurements );
}


std::vector<DiscType> ReadRandomDiscs(ifstream& fp)
{
	std::vector<DiscType> discs;
	DiscType d;
	int numRecs;

	discs.clear();
             
	// First line is number of records
	fp >> numRecs;  skip_line(fp);

	// Skip header line
	skip_line(fp);

	// Next lines are Disc records
	for(int i = 0; i < numRecs; i++)
	{
		fp >> d.pt.phi >> d.pt.theta >> d.radius;
	 
		discs.push_back(d);

		// Skip the rest of the line
		skip_line(fp);

	}
	return( discs );
}

std::vector<PolyType> ReadRandomPolys(ifstream& fp)
{
	std::vector<PolyType> polys;
	PolyType p;
	pointing pt;
	int numRecs,numPoints;
	double nP,nT;

	polys.clear();
             
	// First line is number of records
	fp >> numRecs;  skip_line(fp);

	// Skip header line
	skip_line(fp);

	// Next lines are Poly records
	for(int i = 0; i < numRecs; i++)
	{
		fp >> numPoints;
		p.pts.clear();
		for(int j = 0; j < numPoints; j++ )
		{
			fp >> nP >> nT;
			pt.phi = nP; pt.theta = nT;
			p.pts.push_back(pt);
		}
		polys.push_back(p);

		// Skip the rest of the line
		skip_line(fp);

	}
	return( polys );
}

std::vector<StripType> ReadRandomStrips(ifstream& fp)
{
	std::vector<StripType> strips;
	StripType s;
	int numRecs;

	strips.clear();
             
	// First line is number of records
	fp >> numRecs;  skip_line(fp);

	// Skip header line
	skip_line(fp);

	// Next lines are Strip records
	for(int i = 0; i < numRecs; i++)
	{
		fp >> s.theta1 >> s.theta2;
	 
		strips.push_back(s);

		// Skip the rest of the line
		skip_line(fp);

	}
	return( strips );
}


std::vector<NeighborType> ReadRandomNeighbors(ifstream& fp)
{
	std::vector<NeighborType> neighbors;
	NeighborType n;
	int numRecs;

	neighbors.clear();
             
	// First line is number of records
	fp >> numRecs;  skip_line(fp);

	// Skip header line
	skip_line(fp);

	// Next lines are Neighbor Index records
	for(int i = 0; i < numRecs; i++)
	{
		fp >> n.pt.phi >> n.pt.theta;
	 
		neighbors.push_back(n);

		// Skip the rest of the line
		skip_line(fp);

	}
	return( neighbors );
}


//######### COMPARE HPX & MRH QUERY RESULTS #########

void AnalyzeResults
(
 std::vector<Measurement> foundHPX,
 std::vector<Measurement> foundMRH,
 int& Matches, 
 int& HPXUnique,
 int& MRHUnique,
 bool verbose
 )
{
	bool matchedYN = true;
	pointing pt;
	float p,z;
	int j,k;
	Matches = 0;
	HPXUnique = 0;
	MRHUnique = 0;

	if(verbose){
		cout << "\nHPX Found Measurement Indices:\n";
		cout << "#,REC,PHI,THETA,P,Z\n";
	}

	// Compare each HPX reported data index to list of
	// those found by MRH. Will get number of HPX-MRH matches
	// as well as unique HPX found.
	for( j = 0; j < foundHPX.size(); j++ ) {
		matchedYN = false;
		for( k = 0; k < foundMRH.size(); k++ ) {
			// Check for match
			if( foundMRH[k].equals(foundHPX[j]) ) {
				matchedYN = true;
			}
		}
		// Measurement index found in HPX but NOT MRH!
		if( matchedYN == false ) {
			HPXUnique++;
		}
		if(verbose){
			// Output the HPX query found data points
			pt = foundHPX[j].pt;
			p = pt.phi/pi;
			z = cos(pt.theta);
			cout << j+1 << "," 
				 << foundHPX[j].rec << "," 
				 << foundHPX[j].pt.phi << ","
				 << foundHPX[j].pt.theta << ","
				 << p << "," 
				 << z << endl;
		}
	}
	
	// Now go the other way and compare the MRH found list to HPX found
	// list to discover those data indices ONLY found in MRH.
	if(verbose){
		cout << "\nMRH Found Measurement Indices:\n";
		cout << "#,REC,PHI,THETA,P,Z\n";
	}

	for( j = 0; j < foundMRH.size(); j++ ) {
		matchedYN = false;
		for( k = 0; k < foundHPX.size(); k++ ) {
			// Check for match
			if( foundMRH[j].equals(foundHPX[k]) ) {
				Matches++;
				matchedYN = true;
			}
		}
		// Measurement index found in MRH but NOT HPX!
		if( matchedYN == false ) {
			MRHUnique++;
		}
		if(verbose){
			pt = foundMRH[j].pt;
			p = pt.phi/pi;
			z = cos(pt.theta);
			cout << j+1 << "," 
				 << foundMRH[j].rec << ","
				 << foundMRH[j].pt.phi << ","
				 << foundMRH[j].pt.theta << ","
				 << p << "," 
				 << z << endl;
		}
	}
	if(verbose){
		cout << "\nMatches: " << Matches
			<< " MRHOnly: " << MRHUnique
			<< " HPXOnly: " << HPXUnique << endl;
	}
}


void MRHBenchMarkFINAL
(
 std::string archiveFile,
 std::string mapArchiveFile, 
 std::vector<std::string> queryFilesList, 
 std::string outFile,
 std::string outFile2,
 std::string outPath,
 int numTrials,
 int verbose,
 int maxDepthMRH
)
{
	ofstream fp;
	ifstream ip;
	int i,j,k;
	int64 order;
	uint64 start_s,stop_s;
	uint64 start_test,stop_test;
	std::vector<Measurement> ms;
	std::vector<Measurement> found;
	std::vector<DiscType> discs;
	std::vector<PolyType> polys;
	std::vector<StripType> strips;
	std::vector<NeighborType> neighbors;
	pointing pt,pt0,p0proj;
	std::vector<pointing> p;
	int64 hpxid;
	fix_arr<int64,8> pixsetN;
	rangeset<int64> pixset;
	std::vector<int64> foundV;
	int	totalMrhFound,totalHpxFound,totalMatches,totalMrhUnique,totalHpxUnique;
	int nMRHnodes,nHPXnodes;
	int Matches = 0;
	int MRHUnique = 0;
	int HPXUnique = 0;
	std::vector<Measurement> foundMRH,foundHPX;
	double x,y,azi,rk;
	double avgQueryTime = 0.0;
	double avgPtInsertTime = 0.0;
	double avgQueryTrials = 0.0;
	double TotalTimeCurQuery = 0.0;
	double TotalTimeAllQuery = 0.0;

	GeographicLib::Gnomonic g(GeographicLib::Geodesic(1.0,0.0));


// BEGIN: RECORD RUNTIME FOR FULL BENCHMARK CODE
	start_test = GetTimeMs64();

// Open Results Output File
	cout << "Open Results Output File: " << outFile.c_str() << " ...\n";
	fp.open(outFile.c_str());
	if(fp.fail()){
	  cout << "Unable to open " << outFile.c_str() << "!" << endl;
	  exit(1);     
	}

// Load Disc Query File
	cout << "Begin Loading Disc Queries from: " << queryFilesList[0].c_str() << " ...\n";
	ip.open(queryFilesList[0].c_str());
	if(ip.fail()){
	  cout << "Unable to open " << queryFilesList[0].c_str() << "!" << endl;
	  exit(1);     
	}
	discs = ReadRandomDiscs(ip);
	ip.close();
	cout << "	Completed loading Disc Queries File!\n\n";

// Load Poly Query File
	cout << "Begin Loading Poly Queries from: " << queryFilesList[1].c_str() << " ...\n";
	ip.open(queryFilesList[1].c_str());
	if(ip.fail()){
	  cout << "Unable to open " << queryFilesList[1].c_str() << "!" << endl;
	  exit(1);     
	}
	polys = ReadRandomPolys(ip);
	ip.close();
	cout << "	Completed loading Poly Queries File!\n\n";

// Load Strip Query File
	cout << "Begin Loading Strip Queries from: " << queryFilesList[2].c_str() << " ...\n";
	ip.open(queryFilesList[2].c_str());
	if(ip.fail()){
	  cout << "Unable to open " << queryFilesList[2].c_str() << "!" << endl;
	  exit(1);     
	}
	strips = ReadRandomStrips(ip);
	ip.close();
	cout << "	Completed loading Strip Queries File!\n\n";

// Load Neighbor Query File
	cout << "Begin Loading Neighbors Queries from: " << queryFilesList[3].c_str() << " ...\n";
	ip.open(queryFilesList[3].c_str());
	if(ip.fail()){
	  cout << "Unable to open " << queryFilesList[3].c_str() << "!" << endl;
	  exit(1);     
	}
	neighbors = ReadRandomNeighbors(ip);
	ip.close();
	cout << "	Completed loading Strip Queries File!\n\n";

// Construct MRH from Archive Files
	cout << "##########################################" << endl;
	cout << "#### MRH BUILD FROM ARCHIVE BENCHMARK ####" << endl;
	cout << "##########################################" << endl << endl;
	
	fp << "##########################################" << "\n";
	fp << "#### MRH BUILD FROM ARCHIVE BENCHMARK ####" << "\n";
	fp << "##########################################" << "\n" << "\n";

	cout << "Begin Build MRH From Archive...\n";
	fp << "Begin Build MRH From Archive...\n";
	MultiResHpx_Map<Measurement> mMRH(maxDepthMRH,NEST);

	order = maxDepthMRH;
	
	// START Benchmark
	start_s=GetTimeMs64();
	mMRH.LoadFromFile(archiveFile);
	mMRH.LoadMapFromArchive(mapArchiveFile);
	// STOP Benchmark
	stop_s=GetTimeMs64();

	// Report Trial Results
	cout << std::setprecision(20) << "	Completed MRH Build From Archive in " << double(stop_s-start_s)  << " ms!\n";

	// Write out Trial Result to output file
	fp << std::setprecision(20) << "	Completed MRH Build From Archive in " << double(stop_s-start_s)  << " ms!\n";

#ifdef DO_DISC
//##### BEGIN DISC QUERY BENCHMARK #####
	cout << endl << endl;
	cout << "##################################" << endl;
	cout << "#### MRH DISC QUERY BENCHMARK ####" << endl;
	cout << "##################################" << endl << endl;

	fp << "\n\n";
	fp << "##################################" << "\n";
	fp << "#### MRH DISC QUERY BENCHMARK ####" << "\n";
	fp << "##################################" << "\n" << "\n";

	// START Benchmark
	cout << "Begin " << discs.size() << " Disc Queries..." << endl;
	TotalTimeAllQuery = 0.0;
	fp << "Query #, Num Trials, TotalTime, MeanTime\n";
	for( i = 0; i < discs.size(); i++ )
	{
		TotalTimeCurQuery = 0.0;
		start_s=GetTimeMs64();
		for( j = 0; j < numTrials; j++)
		{          
			found = mMRH.QueryDisc(discs[i].pt,discs[i].radius);
		}
		// STOP Benchmark
		stop_s=GetTimeMs64();
		avgQueryTime = double(stop_s-start_s)/double(numTrials);
		TotalTimeCurQuery = (double(stop_s-start_s));
		TotalTimeAllQuery += TotalTimeCurQuery;

		// Report Query Results, Table format
		fp << i+1 << "," << numTrials << "," << TotalTimeCurQuery << "," << avgQueryTime << "\n";

		// Report Query Results
		cout << "	Completed " << numTrials << " trials of Disc Query " << i+1 << " in " << TotalTimeCurQuery 
			<< " ms! Average Query: " << avgQueryTime << " ms\n\n";
	}
	cout << "End " << discs.size() << " Disc Queries!" << endl;
	cout << endl << endl;

	// Report Results of All Trials
	cout << std::setprecision(20) << "\nCompleted ALL Trials of " << numTrials*discs.size() << " Disc Queries in " << TotalTimeAllQuery 
		<< " ms! Average Query: " << TotalTimeAllQuery/double(numTrials*discs.size()) << " ms\n";
	fp << std::setprecision(20) << "ALL," << numTrials*discs.size() << "," << TotalTimeAllQuery << "," << TotalTimeAllQuery/double(numTrials*discs.size()) << "\n";
//##### END DISC QUERY BENCHMARK #####
#endif

#ifdef DO_POLY
//##### BEGIN POLY QUERY BENCHMARK #####
	cout << endl << endl;
	cout << "##################################" << endl;
	cout << "#### MRH POLY QUERY BENCHMARK ####" << endl;
	cout << "##################################" << endl << endl;
	
	fp << "\n\n";
	fp << "##################################" << "\n";
	fp << "#### MRH POLY QUERY BENCHMARK ####" << "\n";
	fp << "##################################" << "\n" << "\n";

	// START Benchmark
	cout << "Begin " << polys.size() << " Poly Queries..." << endl;
	fp << "Query #, Num Trials, TotalTime, MeanTime\n";
	TotalTimeAllQuery = 0.0;
	for( i = 0; i < polys.size(); i++ )
	{
		TotalTimeCurQuery = 0.0;
		start_s=GetTimeMs64();
		for( j = 0; j < numTrials; j++)
		{          
			found = mMRH.QueryPolygon(polys[i].pts);
		}
		// STOP Benchmark
		stop_s=GetTimeMs64();
		avgQueryTime = double(stop_s-start_s)/double(numTrials);
		TotalTimeCurQuery = (double(stop_s-start_s));
		TotalTimeAllQuery += TotalTimeCurQuery;

		// Report Query Results, Table format
		fp << i+1 << "," << numTrials << "," << TotalTimeCurQuery << "," << avgQueryTime << "\n";

		// Report Trial Results
		cout << "	Completed " << numTrials << " trials of Poly Query " << i+1 << " in " << TotalTimeCurQuery 
			<< " ms! Average Query: " << avgQueryTime << " ms\n\n";
	}
	cout << "End " << polys.size() << " Poly Queries!" << endl;

	// Report Results of All Trials
	cout << std::setprecision(20) << "\nCompleted ALL Trials of " << numTrials*polys.size() << " Poly Queries in " << TotalTimeAllQuery 
		<< " ms! Average Query: " <<TotalTimeAllQuery/double(numTrials*polys.size()) << " ms\n";
	fp << std::setprecision(20) << "ALL," << numTrials*polys.size() << "," << TotalTimeAllQuery << "," << TotalTimeAllQuery/double(numTrials*polys.size()) << "\n";

//##### END POLY QUERY BENCHMARK #####
#endif


#ifdef DO_STRIP
//##### BEGIN STRIP QUERY BENCHMARK #####
	cout << endl << endl;
	cout << "##################################" << endl;
	cout << "#### MRH STRIP QUERY BENCHMARK ####" << endl;
	cout << "##################################" << endl << endl;

	fp << "\n\n";
	fp << "##################################" << "\n";
	fp << "#### MRH STRIP QUERY BENCHMARK ####" << "\n";
	fp << "##################################" << "\n" << "\n";

	// START Benchmark
	cout << "Begin " << strips.size() << " Strip Queries..." << endl;
	fp << "Query #, Num Trials, TotalTime, MeanTime\n";
	TotalTimeAllQuery = 0.0;
	for( i = 0; i < strips.size(); i++ )
	{
		TotalTimeCurQuery = 0.0;
		start_s=GetTimeMs64();
		for( j = 0; j < numTrials; j++)
		{          
			found = mMRH.QueryStrip(strips[i].theta1,strips[i].theta2);
		}
		// STOP Benchmark
		stop_s=GetTimeMs64();
		avgQueryTime = double(stop_s-start_s)/double(numTrials);
		TotalTimeCurQuery = (double(stop_s-start_s));
		TotalTimeAllQuery += TotalTimeCurQuery;

		// Report Query Results, Table format
		fp << i+1 << "," << numTrials << "," << TotalTimeCurQuery << "," << avgQueryTime << "\n";

		// Report Trial Results
		cout << "	Completed " << numTrials << " trials of Strip Query " << i+1 << " in " << TotalTimeCurQuery 
			<< " ms! Average Query: " << avgQueryTime << " ms\n\n";
	}
	cout << "End " << strips.size() << " Strip Queries!" << endl;

	// Report Results of All Trials
	cout << std::setprecision(20) << "\nCompleted ALL Trials of " << numTrials*strips.size() << " Strip Queries in " << TotalTimeAllQuery 
		<< " ms! Average Query: " << TotalTimeAllQuery/double(numTrials*strips.size()) << " ms\n";
	fp << std::setprecision(20) << "ALL," << numTrials*strips.size() << "," << TotalTimeAllQuery << "," << TotalTimeAllQuery/double(numTrials*strips.size()) << "\n";
	
//##### BEGIN STRIP QUERY BENCHMARK #####
#endif


#ifdef DO_NEIGHBORS
//##### BEGIN NEIGHBOR QUERY BENCHMARK #####
	cout << endl << endl;
	cout << "#######################################" << endl;
	cout << "#### MRH NEIGHBORS QUERY BENCHMARK ####" << endl;
	cout << "#######################################" << endl << endl;

	fp << "\n\n";
	fp << "#######################################" << "\n";
	fp << "#### MRH NEIGHBORS QUERY BENCHMARK ####" << "\n";
	fp << "#######################################" << "\n" << "\n";

	// START Benchmark
	cout << "Begin " << neighbors.size() << " Neighbor Queries..." << endl;
	fp << "Query #, Num Trials, TotalTime, MeanTime\n";
	TotalTimeAllQuery = 0.0;
	for( i = 0; i < neighbors.size(); i++ )
	{
		TotalTimeCurQuery = 0.0;
		start_s=GetTimeMs64();
		for( j = 0; j < numTrials; j++)
		{          
			found = mMRH.Neighbors(neighbors[i].pt,order);
		}
		// STOP Benchmark
		stop_s=GetTimeMs64();
		avgQueryTime = double(stop_s-start_s)/double(neighbors.size());
		TotalTimeCurQuery = (double(stop_s-start_s));
		TotalTimeAllQuery += TotalTimeCurQuery;

		// Report Query Results, Table format
		fp << i+1 << "," << numTrials << "," << TotalTimeCurQuery << "," << avgQueryTime << "\n";

		// Report Trial Results
		cout << "	Completed " << numTrials << " trials of Neighbor Query " << i+1 << " in " << TotalTimeCurQuery
			<< " ms! Average Query: " << avgQueryTime << " ms\n\n";
	}
	cout << "End " << neighbors.size() << " Neighbor Queries!" << endl;

	// Report Results of All Trials
	cout << std::setprecision(20) << "\nCompleted ALL Trials of " << numTrials*neighbors.size() << " Neighbor Queries in " << TotalTimeAllQuery 
		<< " ms! Average Query: " << TotalTimeAllQuery/double(numTrials*neighbors.size()) << " ms\n";
	fp << std::setprecision(20) << "ALL," << numTrials*neighbors.size() << "," << TotalTimeAllQuery << "," << TotalTimeAllQuery/double(numTrials*neighbors.size()) << "\n";
//##### END NEIGHBOR QUERY BENCHMARK #####
#endif

// END: RECORD RUNTIME FOR FULL BENCHMARK CODE
	stop_test = GetTimeMs64();

	cout << std::setprecision(20) << "\n####### TOTAL RUNTIME: " << double(stop_test-start_test) << " ms #######\n";
	fp << std::setprecision(20) << "\n####### TOTAL RUNTIME: " << double(stop_test-start_test) << " ms #######\n";

	fp.close();

}


void HPXBenchMarkFINAL(std::string dataInputFile, std::vector<std::string> queryFilesList, std::string outFile, int numTrials)
{
	ofstream fp;
	ifstream ip;
	int i,j,k,l;
	bool over_horizon = false;
	int64 order,hpxid;
	fix_arr<int64,8> pixsetN;
	rangeset<int64> pixset;
	std::vector<int64> foundV;

	uint64 start_s,stop_s;
	uint64 start_test,stop_test;
	std::vector<Measurement> ms;
	std::vector<Measurement> found;
	std::vector<DiscType> discs;
	std::vector<PolyType> polys;
	std::vector<StripType> strips;
	std::vector<NeighborType> neighbors;
	pointing pt,pt0,p0proj;
	std::vector<pointing> p;
	double avgQueryTime = 0.0;
	double avgPtInsertTime = 0.0;
	double avgQueryTrials = 0.0;
	double TotalTimeCurQuery = 0.0;
	double TotalTimeAllQuery = 0.0;
	double x,y,azi,rk;

	GeographicLib::Gnomonic g(GeographicLib::Geodesic(1.0,0.0));


// BEGIN: RECORD RUNTIME FOR FULL BENCHMARK CODE
	start_test = GetTimeMs64();


// Load Measurements File
// Open and Parse Data Input File
	cout << "Begin Loading Data Measurements from: " << dataInputFile.c_str() << " ...\n";
	ip.open(dataInputFile.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << dataInputFile.c_str() << "!" << endl;
	  exit(1);     
	}
	ms = ReadRandomMeasurements(ip,order);
	ip.close();
	cout << "	Completed loading Data Measurements File!\n\n";

// Open Results Output File
	cout << "Open Results Output File: " << outFile.c_str() << " ...\\nn";
	fp.open(outFile.c_str());
	if(fp.fail()){
	  cout << "Unable to open " << outFile.c_str() << "!" << endl;
	  exit(1);     
	}

// Load Disc Query File
	cout << "Begin Loading Disc Queries from: " << queryFilesList[0].c_str() << " ...\\nn";
	ip.open(queryFilesList[0].c_str());
	if(ip.fail()){
	  cout << "Unable to open " << queryFilesList[0].c_str() << "!" << endl;
	  exit(1);     
	}
	discs = ReadRandomDiscs(ip);
	ip.close();
	cout << "	Completed loading Disc Queries File!\n\n";

// Load Poly Query File
	cout << "Begin Loading Poly Queries from: " << queryFilesList[1].c_str() << " ...\n";
	ip.open(queryFilesList[1].c_str());
	if(ip.fail()){
	  cout << "Unable to open " << queryFilesList[1].c_str() << "!" << endl;
	  exit(1);     
	}
	polys = ReadRandomPolys(ip);
	ip.close();
	cout << "	Completed loading Poly Queries File!\n\n";

// Load Neighbor Query File
	cout << "Begin Loading Neighbors Queries from: " << queryFilesList[2].c_str() << " ...\n";
	ip.open(queryFilesList[2].c_str());
	if(ip.fail()){
	  cout << "Unable to open " << queryFilesList[2].c_str() << "!" << endl;
	  exit(1);     
	}
	neighbors = ReadRandomNeighbors(ip);
	ip.close();
	cout << "	Completed loading Strip Queries File!\n\n";

// POINT INSERTION BENCHMARK
// Build HPX Data Structure
	cout << "#######################################" << endl;
	cout << "#### HPX POINT INSERTION BENCHMARK ####" << endl;
	cout << "#######################################" << endl << endl;
	
	fp << "#######################################" << "\n";
	fp << "#### HPX POINT INSERTION BENCHMARK ####" << "\n";
	fp << "#######################################" << "\n" << "\n";

	cout << "Begin " << ms.size() << " MRH Point Insertions...\n";
	fp << "Begin " << ms.size() << " MRH Point Insertions...\n";
	Healpix_Map<Measurement> mHPX(order,NEST);
	// START Benchmark
	start_s=GetTimeMs64();
	for( i = 0; i < ms.size(); i++)
	{
		// Compute HEALPix index of next GIS longitude,latitude,data index tuple
		hpxid = mHPX.ang2pix(ms[i].pt);
		// Insert data index associated with
		mHPX[hpxid] = ms[i]; 
	}
	// STOP Benchmark
	stop_s=GetTimeMs64();
	avgPtInsertTime = double(stop_s-start_s)/double(ms.size());

	// Report Trial Results
	cout << std::setprecision(20) << "	Completed " << ms.size() << " Point Insertions in " << double(stop_s-start_s) 
		<< " ms! Average HPX Point Insertion: " << avgPtInsertTime << " ms\n";

	// Write out Trial Result to output file
	fp << std::setprecision(20) << "	Completed " << ms.size() << " Point Insertions in " << double(stop_s-start_s) 
		<< " ms! Average HPX Point Insertion: " << avgPtInsertTime << " ms\n";

#ifdef DO_DISC
//##### BEGIN DISC QUERY BENCHMARK #####
	cout << endl << endl;
	cout << "##################################" << endl;
	cout << "#### HPX DISC QUERY BENCHMARK ####" << endl;
	cout << "##################################" << endl << endl;

	fp << "\n\n";
	fp << "##################################" << "\n";
	fp << "#### HPX DISC QUERY BENCHMARK ####" << "\n";
	fp << "##################################" << "\n" << "\n";

	// START Benchmark
	cout << "Begin " << discs.size() << " Disc Queries..." << endl;
	//fp << "Begin " << discs.size() << " Disc Queries..." << endl;
	fp << "Query #, Num Trials, TotalTime, MeanTime\n";
	TotalTimeAllQuery = 0.0;
	for( i = 0; i < discs.size(); i++ )
	{
		TotalTimeCurQuery = 0.0;
		start_s=GetTimeMs64();
		for( j = 0; j < numTrials; j++)
		{          
			mHPX.query_disc_inclusive(discs[i].pt,discs[i].radius,pixset);
			
			// Convert cell rangeset to list of individual
			// cells.
			pixset.toVector(foundV);

			// Spin through HPX's found indices keep only the
			// non empty indices that pass the point-in-disc test
			for( k=0; k < foundV.size(); k++ ) {
				if( mHPX[foundV[k]].rec != EMPTY ) {
					if( IsPointInDisc(discs[i].pt,discs[i].radius,mHPX[foundV[k]].pt) ){
						found.push_back(mHPX[foundV[k]]);
					}
				} 
			}
		}
		// STOP Benchmark
		stop_s=GetTimeMs64();
		avgQueryTime = double(stop_s-start_s)/double(numTrials);
		TotalTimeCurQuery = (double(stop_s-start_s));
		TotalTimeAllQuery += TotalTimeCurQuery;

		// Report Query Results, Table format
		fp << i+1 << "," << numTrials << "," << TotalTimeCurQuery << "," << avgQueryTime << "\n";

		// Report Query Results
		cout << "	Completed " << numTrials << " trials of Disc Query " << i+1 << " in " << TotalTimeCurQuery 
			<< " ms! Average Query: " << avgQueryTime << " ms\n\n";
	}
	cout << "End " << discs.size() << " Disc Queries!" << endl;
	cout << endl << endl;

	// Report Results of All Trials
	cout << std::setprecision(20) << "\nCompleted ALL Trials of " << numTrials*discs.size() << " Disc Queries in " << TotalTimeAllQuery 
		<< " ms! Average Query: " << TotalTimeAllQuery/double(numTrials*discs.size()) << " ms\n";
	fp << std::setprecision(20) << "ALL," << numTrials*discs.size() << "," << TotalTimeAllQuery << "," << TotalTimeAllQuery/double(numTrials*discs.size()) << "\n";
//##### END DISC QUERY BENCHMARK #####
#endif

#ifdef DO_POLY
//##### BEGIN POLY QUERY BENCHMARK #####
	cout << endl << endl;
	cout << "##################################" << endl;
	cout << "#### HPX POLY QUERY BENCHMARK ####" << endl;
	cout << "##################################" << endl << endl;
	
	fp << "\n\n";
	fp << "##################################" << "\n";
	fp << "#### HPX POLY QUERY BENCHMARK ####" << "\n";
	fp << "##################################" << "\n" << "\n";

	// START Benchmark
	cout << "Begin " << polys.size() << " Poly Queries..." << endl;
	fp << "Query #, Num Trials, TotalTime, MeanTime\n";
	TotalTimeAllQuery = 0.0;
	for( i = 0; i < polys.size(); i++ )
	{
		TotalTimeCurQuery = 0.0;
		start_s=GetTimeMs64();
		for( j = 0; j < numTrials; j++)
		{          
			mHPX.query_polygon_inclusive(polys[i].pts,pixset,1);
			
			// Convert cell rangeset to list of individual
			// cells.
			pixset.toVector(foundV);

			// Spin through HPX's found indices keep only the
			// non empty indices that pass the point-in-disc test
			for( k=0; k < foundV.size(); k++ ) {
				if( mHPX[foundV[k]].rec != EMPTY ) {
					// Convert point to GIS degrees longitude,latitude
					pt0 = HPXtoGIS(mHPX[foundV[k]].pt);
					pt0.phi *= rad2degr;
					pt0.theta *= rad2degr;
					// Do Gnomonic Projection of point onto 2D flat surface
					g.Forward(pt0.theta,pt0.phi,pt0.theta,pt0.phi,x,y,azi,rk);
					p0proj.theta = x;
					p0proj.phi = y;

					// Use GeographicLib to do ellipsoidal gnomonic projection of query
					// polygon so can then use standard point-in-polygon test.
					p.clear();
					over_horizon = false;
					for( l = 0; l < polys[i].pts.size(); l++ ) {
						pt = HPXtoGIS(polys[i].pts[l]);
						pt.phi *= rad2degr;
						pt.theta *= rad2degr;
						g.Forward(pt0.theta,pt0.phi,pt.theta,pt.phi,x,y,azi,rk);
						p.push_back(pointing(y,x));
						// Flag over the horizon point
						if(rk < 0 ) {
							over_horizon = true;
							l = polys[i].pts.size(); // exit the loop
						}

					}

					// If found point is flagged as over the horizon then reject, else
					// apply point-in-polygon test to filter out data points outside
					// polygon query
					if( over_horizon == false ) {
						if( IsPointInPoly2(p0proj,p) ){
							found.push_back(mHPX[foundV[k]]);
						}
					}
				} 
			}
		}
		// STOP Benchmark
		stop_s=GetTimeMs64();
		avgQueryTime = double(stop_s-start_s)/double(numTrials);
		TotalTimeCurQuery = (double(stop_s-start_s));
		TotalTimeAllQuery += TotalTimeCurQuery;

		// Report Query Results, Table format
		fp << i+1 << "," << numTrials << "," << TotalTimeCurQuery << "," << avgQueryTime << "\n";

		// Report Trial Results
		cout << "	Completed " << numTrials << " trials of Poly Query " << i+1 << " in " << TotalTimeCurQuery 
			<< " ms! Average Query: " << avgQueryTime << " ms\n\n";
	}
	cout << "End " << polys.size() << " Poly Queries!" << endl;

	// Report Results of All Trials
	cout << std::setprecision(20) << "\nCompleted ALL Trials of " << numTrials*polys.size() << " Poly Queries in " << TotalTimeAllQuery 
		<< " ms! Average Query: " <<TotalTimeAllQuery/double(numTrials*polys.size()) << " ms\n";
	fp << std::setprecision(20) << "ALL," << numTrials*polys.size() << "," << TotalTimeAllQuery << "," << TotalTimeAllQuery/double(numTrials*polys.size()) << "\n";
//##### END POLY QUERY BENCHMARK #####
#endif

#ifdef DO_NEIGHBORS
//##### BEGIN NEIGHBOR QUERY BENCHMARK #####
	cout << endl << endl;
	cout << "#######################################" << endl;
	cout << "#### HPX NEIGHBORS QUERY BENCHMARK ####" << endl;
	cout << "#######################################" << endl << endl;

	fp << "\n\n";
	fp << "#######################################" << "\n";
	fp << "#### HPX NEIGHBORS QUERY BENCHMARK ####" << "\n";
	fp << "#######################################" << "\n" << "\n";

	// START Benchmark
	cout << "Begin " << neighbors.size() << " Neighbor Queries..." << endl;
	fp << "Query #, Num Trials, TotalTime, MeanTime\n";
	TotalTimeAllQuery = 0.0;
	for( i = 0; i < neighbors.size(); i++ )
	{
		TotalTimeCurQuery = 0.0;
		start_s=GetTimeMs64();
		for( j = 0; j < numTrials; j++)
		{          
			hpxid = mHPX.ang2pix(neighbors[i].pt);
			mHPX.neighbors(hpxid,pixsetN);
			
			// Spin through HPX's found indices keep only the
			// non empty indices that pass the point-in-disc test
			for( k=0; k < 8; k++ ) {
				if( mHPX[pixsetN[k]].rec != EMPTY ) {
					found.push_back(mHPX[pixsetN[k]]);
				} 
			}
		}
		// STOP Benchmark
		stop_s=GetTimeMs64();
		avgQueryTime = double(stop_s-start_s)/double(neighbors.size());
		TotalTimeCurQuery = (double(stop_s-start_s));
		TotalTimeAllQuery += TotalTimeCurQuery;

		// Report Query Results, Table format
		fp << i+1 << "," << numTrials << "," << TotalTimeCurQuery << "," << avgQueryTime << "\n";

		// Report Trial Results
		cout << "	Completed " << numTrials << " trials of Neighbor Query " << i+1 << " in " << TotalTimeCurQuery
			<< " ms! Average Query: " << avgQueryTime << " ms\n\n";
	}
	cout << "End " << neighbors.size() << " Neighbor Queries!" << endl;

	// Report Results of All Trials
	cout << std::setprecision(20) << "\nCompleted ALL Trials of " << numTrials*neighbors.size() << " Neighbor Queries in " << TotalTimeAllQuery 
		<< " ms! Average Query: " << TotalTimeAllQuery/double(numTrials*neighbors.size()) << " ms\n";
	fp << std::setprecision(20) << "ALL," << numTrials*neighbors.size() << "," << TotalTimeAllQuery << "," << TotalTimeAllQuery/double(numTrials*neighbors.size()) << "\n";
//##### END NEIGHBOR QUERY BENCHMARK #####
#endif

// END: RECORD RUNTIME FOR FULL BENCHMARK CODE
	stop_test = GetTimeMs64();

	cout << std::setprecision(20) << "\n####### TOTAL RUNTIME: " << double(stop_test-start_test) << " ms #######\n";
	fp << std::setprecision(20) << "\n####### TOTAL RUNTIME: " << double(stop_test-start_test) << " ms #######\n";

	fp.close();
}





void HPXStripBenchMarkFINAL(std::string dataInputFile, std::string stripQueryFile, std::string outFile, int numTrials)
{
	ofstream fp;
	ifstream ip;
	int i,j,k;
	int64 order,hpxid;
	rangeset<int64> pixset;
	std::vector<int64> foundV;

	uint64 start_s,stop_s;
	uint64 start_test,stop_test;
	std::vector<Measurement> ms;
	std::vector<Measurement> found;
	std::vector<StripType> strips;
	pointing pt;
	double avgQueryTime = 0.0;
	double avgPtInsertTime = 0.0;
	double avgQueryTrials = 0.0;
	double TotalTimeCurQuery = 0.0;
	double TotalTimeAllQuery = 0.0;

// BEGIN: RECORD RUNTIME FOR FULL BENCHMARK CODE
	start_test = GetTimeMs64();

// Load Measurements File
// Open and Parse Data Input File
	cout << "Begin Loading Data Measurements from: " << dataInputFile.c_str() << " ...\n";
	ip.open(dataInputFile.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << dataInputFile.c_str() << "!" << endl;
	  exit(1);     
	}
	ms = ReadRandomMeasurements(ip,order);
	ip.close();
	cout << "	Completed loading Data Measurements File!\n\n";

// Open Results Output File
	cout << "Open Results Output File: " << outFile.c_str() << " ...\n";
	fp.open(outFile.c_str());
	if(fp.fail()){
	  cout << "Unable to open " << outFile.c_str() << "!" << endl;
	  exit(1);     
	}

	cout << "Create new HPX Data Structure that's indexed using RING scheme...\n";
	fp << "Create new HPX Data Structure that's indexed using RING scheme...\n";

	cout << "Begin adding " << ms.size() << " records to the HPX Data Structure...\n";
	fp << "Begin adding " << ms.size() << " records to the HPX Data Structure...\n";
	Healpix_Map<Measurement> mHPX2(order,RING);
	for( i = 0; i < ms.size(); i++)
	{
		if( i%1000 == 0 ) {
			cout << "	Inserting Point " << i+1 << " of " << ms.size() << endl;
			fp << "	Inserting Point " << i+1 << " of " << ms.size() << "\n";
		}
		// Compute HEALPix index of next GIS longitude,latitude,data index tuple
		hpxid = mHPX2.ang2pix(ms[i].pt);

		// Insert data index associated with
		mHPX2[hpxid] = ms[i]; 
	}
	cout << "Completed adding " << ms.size() << " records to the HPX Data Structure!\n";
	fp << "Completed adding " << ms.size() << " records to the HPX Data Structure!\n";
	
	// Load Strip Query File
	cout << "Begin Loading Strip Queries from: " << stripQueryFile.c_str() << " ...\n";
	ip.open(stripQueryFile.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << stripQueryFile.c_str() << "!" << endl;
	  exit(1);     
	}
	strips = ReadRandomStrips(ip);
	ip.close();
	cout << "	Completed loading Strip Queries File!\n\n";


#ifdef DO_STRIP
//##### BEGIN STRIP QUERY BENCHMARK #####
	cout << endl << endl;
	cout << "##################################" << endl;
	cout << "#### HPX STRIP QUERY BENCHMARK ####" << endl;
	cout << "##################################" << endl << endl;

	fp << "\n\n";
	fp << "##################################" << "\n";
	fp << "#### HPX STRIP QUERY BENCHMARK ####" << "\n";
	fp << "##################################" << "\n" << "\n";

	// START Benchmark
	cout << "Begin " << strips.size() << " Strip Queries..." << endl;
	fp << "Query #, Num Trials, TotalTime, MeanTime\n";
	TotalTimeAllQuery = 0.0;
	for( i = 0; i < strips.size(); i++ )
	{
		TotalTimeCurQuery = 0.0;
		start_s=GetTimeMs64();
		for( j = 0; j < numTrials; j++)
		{          
			mHPX2.query_strip(strips[i].theta1,strips[i].theta2,true,pixset);
			
			// Convert cell rangeset to list of individual
			// cells.
			pixset.toVector(foundV);

			// Spin through HPX's found indices keep only the
			// non empty indices that pass the point-in-disc test
			for( k=0; k < foundV.size(); k++ ) {
				if( mHPX2[foundV[k]].rec != EMPTY ) {
					if( IsPointInStrip(strips[i].theta1,strips[i].theta2,mHPX2[foundV[k]].pt) ){
						found.push_back(mHPX2[foundV[k]]);
					}
				} 
			}
		}
		// STOP Benchmark
		stop_s=GetTimeMs64();
		avgQueryTime = double(stop_s-start_s)/double(numTrials);
		TotalTimeCurQuery = (double(stop_s-start_s));
		TotalTimeAllQuery += TotalTimeCurQuery;

		// Report Query Results, Table format
		fp << i+1 << "," << numTrials << "," << TotalTimeCurQuery << "," << avgQueryTime << "\n";

		// Report Trial Results
		cout << "	Completed " << numTrials << " trials of Strip Query " << i+1 << " in " << TotalTimeCurQuery 
			<< " ms! Average Query: " << avgQueryTime << " ms\n\n";
	}
	cout << "End " << strips.size() << " Strip Queries!" << endl;

	// Report Results of All Trials
	cout << std::setprecision(20) << "\nCompleted ALL Trials of " << numTrials*strips.size() << " Strip Queries in " << TotalTimeAllQuery 
		<< " ms! Average Query: " << TotalTimeAllQuery/double(numTrials*strips.size()) << " ms\n";
	fp << std::setprecision(20) << "ALL," << numTrials*strips.size() << "," << TotalTimeAllQuery << "," << TotalTimeAllQuery/double(numTrials*strips.size()) << "\n";
//##### BEGIN STRIP QUERY BENCHMARK #####
#endif

// END: RECORD RUNTIME FOR FULL BENCHMARK CODE
	stop_test = GetTimeMs64();

	cout << std::setprecision(20) << "\n####### TOTAL RUNTIME: " << double(stop_test-start_test) << " ms #######\n";
	fp << std::setprecision(20) << "\n####### TOTAL RUNTIME: " << double(stop_test-start_test) << " ms #######\n";

	fp.close();


}

//
// Just compare the MRH vs HPX Query Outputs!
//
void MRHvsHPXComparisionFINAL
(
 std::string archiveFile, 
 std::string mapArchiveFile,
 std::vector<std::string> queryFilesList, 
 std::string outFile2,
 std::string outPath,
 int verbose,
 int maxDepthMRH
)
{
	ofstream fp;
	ifstream ip;
	int i,j,k;
	int64 order;
	uint64 start_s,stop_s;
	uint64 start_test,stop_test;
	std::vector<Measurement> ms;
	std::vector<Measurement> found;
	std::vector<DiscType> discs;
	std::vector<PolyType> polys;
	std::vector<StripType> strips;
	std::vector<NeighborType> neighbors;
	pointing pt,pt0,p0proj;
	std::vector<pointing> p;
	int64 hpxid;
	fix_arr<int64,8> pixsetN;
	rangeset<int64> pixset;
	std::vector<int64> foundV;
	int	totalMrhFound,totalHpxFound,totalMatches,totalMrhUnique,totalHpxUnique;
	int64 totalMrhCritCount,totalHpxCritCount,totalMrhUpsearchCount;
	int Matches = 0;
	int MRHUnique = 0;
	int HPXUnique = 0;
	int nMRHnodes = 0;
	int nHPXnodes = 0;
	bool over_horizon = false;
	std::vector<Measurement> foundMRH,foundHPX;
	double x,y,azi,rk;
	double avgQueryTime = 0.0;
	double avgPtInsertTime = 0.0;
	double avgQueryTrials = 0.0;
	double TotalTimeCurQuery = 0.0;
	double TotalTimeAllQuery = 0.0;

	GeographicLib::Gnomonic g(GeographicLib::Geodesic(1.0,0.0));

// Load Measurements File
	cout << "Begin Loading Data Measurements from: " << queryFilesList[0].c_str() << " ...\n";
	ip.open(queryFilesList[0].c_str());
	if(ip.fail()){
	  cout << "Unable to open " << queryFilesList[0].c_str() << "!" << endl;
	  exit(1);     
	}
	ms = ReadRandomMeasurements(ip,order);
	ip.close();
	cout << "	Completed loading Data Measurements File!\n\n";


// Load Disc Query File
	cout << "Begin Loading Disc Queries from:\n " << queryFilesList[1].c_str() << " ...\n";
	ip.open(queryFilesList[1].c_str());
	if(ip.fail()){
	  cout << "Unable to open " << queryFilesList[1].c_str() << "!" << endl;
	  exit(1);     
	}
	discs = ReadRandomDiscs(ip);
	ip.close();
	cout << "	Completed loading Disc Queries File!\n\n";

// Load Poly Query File
	cout << "Begin Loading Poly Queries from:\n " << queryFilesList[2].c_str() << " ...\n";
	ip.open(queryFilesList[2].c_str());
	if(ip.fail()){
	  cout << "Unable to open " << queryFilesList[2].c_str() << "!" << endl;
	  exit(1);     
	}
	polys = ReadRandomPolys(ip);
	ip.close();
	cout << "	Completed loading Poly Queries File!\n\n";

// Load Strip Query File
	cout << "Begin Loading Strip Queries from:\n " << queryFilesList[3].c_str() << " ...\n";
	ip.open(queryFilesList[3].c_str());
	if(ip.fail()){
	  cout << "Unable to open " << queryFilesList[3].c_str() << "!" << endl;
	  exit(1);     
	}
	strips = ReadRandomStrips(ip);
	ip.close();
	cout << "	Completed loading Strip Queries File!\n\n";

// Load Neighbor Query File
	cout << "Begin Loading Neighbors Queries from:\n " << queryFilesList[4].c_str() << " ...\n";
	ip.open(queryFilesList[4].c_str());
	if(ip.fail()){
	  cout << "Unable to open " << queryFilesList[4].c_str() << "!" << endl;
	  exit(1);     
	}
	neighbors = ReadRandomNeighbors(ip);
	ip.close();
	cout << "	Completed loading Neighbor Queries File!\n\n";

	// For outputing query results to individual files if
	// verbose is set.
	ofstream qp;



// Open Results Output File
	cout << "Open Results Output File: " << outFile2.c_str() << " ...\n\n";
	fp.open(outFile2.c_str());
	if(fp.fail()){
	  cout << "Unable to open " << outFile2.c_str() << "!" << endl;
	  exit(1);     
	}

// Build MultiResHpx Data Structure from Archive
	MultiResHpx_Map<Measurement> mMRH(maxDepthMRH,NEST);
	cout << "Building MRH Data Structure From Archive Files...\n";
	
	cout << "LOAD MLQ FOREST ARCHIVE: \n\t" << archiveFile.c_str() << endl;
	mMRH.LoadFromFile(archiveFile);
	cout << "\tMLQ FOREST BUILD COMPLETE!\n";
	
	cout << "LOAD MRH MAP ARCHIVE: \n\t" << mapArchiveFile.c_str() << endl;
	mMRH.LoadMapFromArchive(mapArchiveFile);
	cout << "\tMRH MAP BUILD COMPLETE!\n";



//##########################################################################################################
//####### MRH VERSUS HPX QUERY OUTPUT COMPARISON  ##########################################################
//##########################################################################################################

	cout << endl;
	cout << "########################################################################\n";
	cout << "##### Now perform MRH versus HPX query output comparison analysis! #####\n";
	cout << "########################################################################\n\n";

// BEGIN: RECORD RUNTIME FOR FULL BENCHMARK CODE
	start_test = GetTimeMs64();

// Build HPX Data Structure
	Healpix_Map<Measurement> mHPX(order,NEST);
	cout << "Begin " << ms.size() << " HPX Point Insertions...\n";
	fp << "Begin " << ms.size() << " HPX Point Insertions...\n";
	for( i = 0; i < ms.size(); i++)
	{
		if( i%1000 == 0 ) {
			cout << "	Inserting Point " << i << " of " << ms.size() << endl;
		}
		// Compute HEALPix index of next GIS longitude,latitude,data index tuple
		hpxid = mHPX.ang2pix(ms[i].pt);
		// Insert data index associated with
		mHPX[hpxid] = ms[i]; 
	}
	cout << "Completed " << ms.size() << " HPX Point Insertions!\n";
	fp << "Completed " << ms.size() << " HPX Point Insertions!\n";

#ifdef DO_DISC
//###### BEGIN COMPARE DISC QUERY RESULTS ######
	cout << "\n\n##### BEGIN DISC QUERY RESULTS COMPARISON #####\n";
	fp << "\n\n##### BEGIN DISC QUERY RESULTS COMPARISON #####\n";
	cout << "Begin " << discs.size() << " Disc Queries..." << endl;
#ifdef CRITCOUNT
	fp << "Query#, Matches, MRH Only, HPX Only, MRH Crit, HPX Crit, MRH Cvr Map, HPX Cvr Map, MRH Upsearch, MRH Cvr Map Cell Res" << "\n";
#else
	fp << "Query#, Matches, MRH Only, HPX Only," << "\n";
#endif
	totalMrhFound = 0;totalHpxFound=0;totalMatches=0;totalMrhUnique=0;totalHpxUnique=0;totalMrhCritCount=0;totalHpxCritCount=0;totalMrhUpsearchCount=0;
	for( i = 0; i < discs.size(); i++)
	{          

		// Open new output file containing data points, query definition
		// and found query points.
		if(verbose) {
			std::ostringstream stream;
			stream << outPath << "/discQuery" << i+1 << ".csv";
			std::string queryOut = stream.str();
			qp.open(queryOut.c_str());
			if(qp.fail()){
				cout << "Unable to open " << queryOut.c_str() << "!" << endl;
				exit(1);     
			}
			// First output the data points
			qp << "Phi,Theta,P,Z,Rec\n";
			for(j=0; j < ms.size(); j++) {
			   // Write out data location to query comparision file
			   qp << ms[j].pt.phi << "," << ms[j].pt.theta << ","
				  << ms[j].pt.phi/pi << "," << cos(ms[j].pt.theta) << ","
				  << ms[j].rec << "\n";
			}
			qp << "\n\n";

			// Next output the query definition
			qp << "Query #" << i+1 << "\n";
			qp << "Phi,Theta,P,Z,Radius\n";
			qp << discs[i].pt.phi << "," << discs[i].pt.theta << ","
			   << discs[i].pt.phi/pi << "," << cos(discs[i].pt.theta) << ","
			   << discs[i].radius << "\n\n";
		}

		Matches = 0; MRHUnique = 0; HPXUnique = 0;
		foundMRH.clear(); foundHPX.clear(); foundV.clear();
#ifdef CRITCOUNT
        HPX_CRIT_COUNT = 0;
		MRH_CRIT_COUNT = 0;
		mHPX.ResetCritCount();
        MRH_UPSEARCH_COUNT = 0;
		mMRH.ResetUpSearchCount();
#endif

		// Do MRH Disc Query
		foundMRH = mMRH.QueryDisc(discs[i].pt,discs[i].radius);
#ifdef CRITCOUNT
		MRH_CRIT_COUNT = mMRH.GetCritCount();
		MRH_COVER_MAP_SIZE = mMRH.GetCoverMapSize();
		MRH_UPSEARCH_COUNT = mMRH.GetUpSearchCount();
		MRH_COVER_MAP_CELL_RES = mMRH.GetCoverMapCellRes();
#endif

		// Do HPX Disc Query
		mHPX.query_disc_inclusive(discs[i].pt,discs[i].radius,pixset);
		// Convert cell rangeset to list of individual
		// cells.
		pixset.toVector(foundV);

#ifdef CRITCOUNT
		HPX_CRIT_COUNT = mHPX.GetCritCount();
		HPX_COVER_MAP_SIZE = foundV.size();
#endif		

		// Spin through HPX's found indices keep only the
		// non empty indices that pass the point-in-disc test
		for( j=0; j < foundV.size(); j++ ) {
			if( mHPX[foundV[j]].rec != EMPTY ) {
				if( IsPointInDisc(discs[i].pt,discs[i].radius,mHPX[foundV[j]].pt) ){
					foundHPX.push_back(mHPX[foundV[j]]);
				}
			}
#ifdef CRITCOUNT
            HPX_CRIT_COUNT += 1;
#endif
		}
		AnalyzeResults(foundHPX,foundMRH,Matches,HPXUnique,MRHUnique,false);
		totalMatches += Matches;
		totalMrhUnique += MRHUnique;
		totalHpxUnique += HPXUnique;
		totalMrhFound += foundMRH.size();
		totalHpxFound += foundHPX.size();
#ifdef CRITCOUNT
		totalMrhCritCount += MRH_CRIT_COUNT;
		totalHpxCritCount += HPX_CRIT_COUNT;
		totalMrhUpsearchCount += MRH_UPSEARCH_COUNT;
#endif

#ifdef CRITCOUNT
		cout << "Query " << i+1 << " Matches: " << Matches << " MRH Only: " << MRHUnique << " HPX Only: " << HPXUnique 
			<< " MRH Crit Count: " << MRH_CRIT_COUNT << " HPX Crit Count: " << HPX_CRIT_COUNT << " MRH Cvr Map: " << MRH_COVER_MAP_SIZE 
			<< " HPX Cvr Map: " << HPX_COVER_MAP_SIZE << " MRH Upsearch Count: " << MRH_UPSEARCH_COUNT << " MRH Cvr Map Cell Res: " << MRH_COVER_MAP_CELL_RES << endl;
		fp << i+1 << "," << Matches << "," << MRHUnique << "," << HPXUnique 
			<< "," << MRH_CRIT_COUNT << "," << HPX_CRIT_COUNT << "," << MRH_COVER_MAP_SIZE << "," << HPX_COVER_MAP_SIZE 
			<< "," << MRH_UPSEARCH_COUNT << "," << MRH_COVER_MAP_CELL_RES << "\n";
#else
		cout << "Query " << i+1 << " Matches: " << Matches << " MRH Only: " << MRHUnique << " HPX Only: " << HPXUnique << endl;
		
		// Output Trial Results in Table format
		fp << i+1 << "," << Matches << "," << MRHUnique << "," << HPXUnique << "\n";
#endif
		// Output Query and Query results to output file
		// Write Out Query # and the Query in phi,theta,p,z
		if(verbose) {
			// Now write out the MRH found results
			qp << "Query # " << i+1 << " MRH Found Results\n";
			qp << "Phi,Theta,P,Z,Rec\n";
			for( j=0; j < foundMRH.size(); j++) {
				qp << foundMRH[j].pt.phi << "," << foundMRH[j].pt.theta << ","
				   << foundMRH[j].pt.phi/pi << "," << cos(foundMRH[j].pt.theta) << ","
				   << foundMRH[j].rec << "\n";
			}
			qp << "\n";
			// Now write out the HPX found results
			qp << "Query # " << i+1 << " HPX Found Results\n";
			qp << "Phi,Theta,P,Z,Rec\n";
			for( j=0; j < foundHPX.size(); j++) {
				qp << foundHPX[j].pt.phi << "," << foundHPX[j].pt.theta << ","
				   << foundHPX[j].pt.phi/pi << "," << cos(foundHPX[j].pt.theta) << ","
				   << foundHPX[j].rec << "\n";	
			}
			qp << "\n\n";
			qp.close();
		}
	}
	cout << "	Finished Computing MRH vs HPX Disc Query Accuracy Test!\n";

	// Compute/get statistics of the final data structures 
	nMRHnodes = mMRH.NumNodes();
	nHPXnodes = mHPX.Npix();
    // Get Min,Max,Average Depth of QuadTrees of each data stucture
	//int maxDepthMRH = mMRH.MaxDepth();
	int minDepthMRH = mMRH.MinDepth();
	int avgDepthMRH = mMRH.AvgDepth();

	cout << "Number Points: " << ms.size() << endl;
	cout << "Number Disc Queries: " << discs.size() << endl;
	cout << "Total MRH Nodes: " << nMRHnodes << " Total HPX Nodes: " << nHPXnodes << endl;
	cout << "Memory Size of MRH: " << mMRH.SizeInMemKb() << " Kb" << " Size of HPX: " << (nHPXnodes*sizeof(Measurement))/1024 << " Kb" << endl;
	cout << "MRH Min Depth: " << minDepthMRH << " Max Depth: " << maxDepthMRH << " Avg Depth: " << avgDepthMRH << endl;
	cout << "MLQ\t#Nodes\tDepth\tAvgDepth\n";
	for( i = 0; i < 12; i++ ) {
		cout << i << "\t" 
			 << mMRH.NumNodesAtMLQ(i) << "\t"
			 << mMRH.DepthAtMLQ(i) << "\t"
			 << mMRH.AvgDepthAtMLQ(i) << "\n";
	}
	cout << "HPX Depth: " << order << endl;
	cout << "Total Matches: " << totalMatches << endl;
	cout << "Total MRH Found: " << totalMrhFound << " Total MRH Unique: " << totalMrhUnique << endl;
	cout << "Total HPX Found: " << totalHpxFound << " Total HPX Unique: " << totalHpxUnique << endl;
#ifdef CRITCOUNT
	cout << "Total MRH Crit Count: " << totalMrhCritCount << endl;
	cout << "Total HPX Crit Count: " << totalHpxCritCount << endl;
	cout << "Total Upsearch Count: " << totalMrhUpsearchCount << endl;
#endif

	fp << "\nSummary:\n";
	fp << "Number Points: " << ms.size() << "\n";
	fp << "Number Disc Queries: " << discs.size() << "\n";
	fp << "Total MRH Nodes: " << nMRHnodes << " Total HPX Nodes: " << nHPXnodes << "\n";
	fp << "Memory Size of MRH: " << mMRH.SizeInMemKb() << " Kb" << " Size of HPX: " << (nHPXnodes*sizeof(Measurement))/1024 << " Kb\n";
	fp << "MRH Min Depth: " << minDepthMRH << " Max Depth: " << maxDepthMRH << " Avg Depth: " << avgDepthMRH << "\n";
	fp << "MLQ\t#Nodes\tDepth\tAvgDepth\n";
	for( i = 0; i < 12; i++ ) {
		fp << i << "\t" 
			 << mMRH.NumNodesAtMLQ(i) << "\t"
			 << mMRH.DepthAtMLQ(i) << "\t"
			 << mMRH.AvgDepthAtMLQ(i) << "\n";
	}
	fp << "HPX Depth: " << order << "\n";
	fp << "Total Matches: " << totalMatches << endl;
	fp << "Total MRH Found: " << totalMrhFound << " Total MRH Unique: " << totalMrhUnique << "\n";
	fp << "Total HPX Found: " << totalHpxFound << " Total HPX Unique: " << totalHpxUnique << "\n";
#ifdef CRITCOUNT
	fp << "Total MRH Crit Count: " << totalMrhCritCount << "\n";
	fp << "Total HPX Crit Count: " << totalHpxCritCount << "\n";
	fp << "Total Upsearch Count: " << totalMrhUpsearchCount << "\n";
#endif
	cout << "\n##### END DISC QUERY RESULTS COMPARISON #####\n";
	fp << "\n##### END DISC QUERY RESULTS COMPARISON #####\n";
//###### END COMPARE DISC QUERY RESULTS ######
#endif



#ifdef DO_POLY
//###### BEGIN COMPARE POLY QUERY RESULTS ######
    cout << "\n\n##### BEGIN POLY QUERY RESULTS COMPARISON #####\n";
	fp << "\n\n##### BEGIN POLY QUERY RESULTS COMPARISON #####\n";
	cout << "Begin " << polys.size() << " Poly Queries..." << endl;
#ifdef CRITCOUNT
	fp << "Query#, Matches, MRH Only, HPX Only, MRH Crit, HPX Crit, MRH Cvr Map, HPX Cvr Map, MRH Upsearch, MRH Cvr Map Cell Res" << "\n";
#else
	fp << "Query#, Matches, MRH Only, HPX Only," << "\n";
#endif
	totalMrhFound = 0;totalHpxFound=0;totalMatches=0;totalMrhUnique=0;totalHpxUnique=0;totalMrhCritCount=0;totalHpxCritCount=0,totalMrhUpsearchCount=0;
	for( i = 0; i < polys.size(); i++)
	{          
		// Open new output file containing data points, query definition
		// and found query points.
		if(verbose) {
			std::ostringstream stream;
			stream << outPath << "/polyQuery" << i+1 << ".csv";
			std::string queryOut = stream.str();
			qp.open(queryOut.c_str());
			if(qp.fail()){
				cout << "Unable to open " << queryOut.c_str() << "!" << endl;
				exit(1);     
			}
			// First output the data points
			qp << "Phi,Theta,P,Z,Rec\n";
			for(j=0; j < ms.size(); j++) {
			   // Write out data location to query comparision file
			   qp << ms[j].pt.phi << "," << ms[j].pt.theta << ","
				  << ms[j].pt.phi/pi << "," << cos(ms[j].pt.theta) << ","
				  << ms[j].rec << "\n";
			}
			qp << "\n\n";

			// Next output the query definition
			qp << "Query #" << i+1 << "\n";
			qp << "Phi,Theta,P,Z,...\n";
			for( j=0; j < polys[i].pts.size(); j++) {
			   qp << polys[i].pts[j].phi << "," << polys[i].pts[j].theta << ","
				  << polys[i].pts[j].phi/pi << "," << cos(polys[i].pts[j].theta) << "\n";
			}
			qp << polys[i].pts[0].phi << "," << polys[i].pts[0].theta << ","
				  << polys[i].pts[0].phi/pi << "," << cos(polys[i].pts[0].theta) << "\n";
			qp << "\n\n";
		}


		Matches = 0; MRHUnique = 0; HPXUnique = 0;
		foundMRH.clear(); foundHPX.clear(); foundV.clear();
#ifdef CRITCOUNT
        HPX_CRIT_COUNT = 0;
		MRH_CRIT_COUNT = 0;
		MRH_UPSEARCH_COUNT = 0;
		MRH_COVER_MAP_CELL_RES = 0;
		mHPX.ResetCritCount();
		mMRH.ResetUpSearchCount(); 
#endif
		// Do MRH Poly Query
		foundMRH = mMRH.QueryPolygon(polys[i].pts);
#ifdef CRITCOUNT
		MRH_CRIT_COUNT = mMRH.GetCritCount();
		MRH_COVER_MAP_SIZE = mMRH.GetCoverMapSize();
		MRH_UPSEARCH_COUNT = mMRH.GetUpSearchCount();
		MRH_COVER_MAP_CELL_RES = mMRH.GetCoverMapCellRes();
#endif

		// Do HPX Poly Query
		mHPX.query_polygon_inclusive(polys[i].pts,pixset,1);
		// Convert cell rangeset to list of individual
		// cells.
		pixset.toVector(foundV);

#ifdef CRITCOUNT
		HPX_CRIT_COUNT = mHPX.GetCritCount();
		HPX_COVER_MAP_SIZE = foundV.size();
#endif		

		// Spin through HPX's found indices keep only the
		// non empty indices that pass the point-in-disc test
		//for( j=0; j < foundV.size(); j++ ) {
		//	if( mHPX[foundV[j]].rec != EMPTY ) {
		//		if( IsPointInPoly2(mHPX[foundV[j]].pt,polys[i].pts) ){
		//			foundHPX.push_back(mHPX[foundV[j]]);
		//		}
		//	} 
		//}
		// Spin through HPX's found indices keep only the
		// non empty indices that pass the point-in-poly test
		for( j=0; j < foundV.size(); j++ ) {
			if( mHPX[foundV[j]].rec != EMPTY ) {
				// Convert point to GIS degrees longitude,latitude
				pt0 = HPXtoGIS(mHPX[foundV[j]].pt);
				pt0.phi *= rad2degr;
				pt0.theta *= rad2degr;
				// Do Gnomonic Projection of point onto 2D flat surface
				g.Forward(pt0.theta,pt0.phi,pt0.theta,pt0.phi,x,y,azi,rk);
				p0proj.theta = x;
				p0proj.phi = y;

				// Use GeographicLib to do ellipsoidal gnomonic projection of query
				// polygon so can then use standard point-in-polygon test.
				p.clear();
				over_horizon = false;
				for( k = 0; k < polys[i].pts.size(); k++ ) {
					pt = HPXtoGIS(polys[i].pts[k]);
					pt.phi *= rad2degr;
					pt.theta *= rad2degr;
					g.Forward(pt0.theta,pt0.phi,pt.theta,pt.phi,x,y,azi,rk);
					p.push_back(pointing(y,x));

					// Flag over the horizon point
					if(rk < 0 ) {
						over_horizon = true;
						k = polys[i].pts.size(); // exit the loop
					}
#ifdef CRITCOUNT
                    HPX_CRIT_COUNT += 1;
#endif
				}

				// If found point is flagged as over the horizon then reject, else
				// apply point-in-polygon test to filter out data points outside
				// polygon query
				if( over_horizon == false ) {
					if( IsPointInPoly2(p0proj,p) ){
						foundHPX.push_back(mHPX[foundV[j]]);
					}
				}
			} 
#ifdef CRITCOUNT
            HPX_CRIT_COUNT += 1;
#endif
		}
		AnalyzeResults(foundHPX,foundMRH,Matches,HPXUnique,MRHUnique,false);
		totalMatches += Matches;
		totalMrhUnique += MRHUnique;
		totalHpxUnique += HPXUnique;
		totalMrhFound += foundMRH.size();
		totalHpxFound += foundHPX.size();
#ifdef CRITCOUNT
		totalMrhCritCount += MRH_CRIT_COUNT;
		totalHpxCritCount += HPX_CRIT_COUNT;
		totalMrhUpsearchCount += MRH_UPSEARCH_COUNT;
#endif


#ifdef CRITCOUNT
		cout << "Query " << i+1 << " Matches: " << Matches << " MRH Only: " << MRHUnique << " HPX Only: " << HPXUnique 
			<< " MRH Crit Count: " << MRH_CRIT_COUNT << " HPX Crit Count: " << HPX_CRIT_COUNT << " MRH Cvr Map: " << MRH_COVER_MAP_SIZE 
			<< " HPX Cvr Map: " << HPX_COVER_MAP_SIZE << " MRH Upsearch Count: " << MRH_UPSEARCH_COUNT << " MRH Cvr Map Cell Res: " << MRH_COVER_MAP_CELL_RES << endl;
		fp << i+1 << "," << Matches << "," << MRHUnique << "," << HPXUnique 
			<< "," << MRH_CRIT_COUNT << "," << HPX_CRIT_COUNT << "," << MRH_COVER_MAP_SIZE << "," << HPX_COVER_MAP_SIZE 
			<< "," << MRH_UPSEARCH_COUNT << "," << MRH_COVER_MAP_CELL_RES << "\n";
#else
		cout << "Query " << i+1 << " Matches: " << Matches << " MRH Only: " << MRHUnique << " HPX Only: " << HPXUnique << endl;
		
		// Output Trial Results in Table format
		fp << i+1 << "," << Matches << "," << MRHUnique << "," << HPXUnique << "\n";
#endif

		// Output Query and Query results to output file
		// Write Out Query # and the Query in phi,theta,p,z
		if(verbose) {
			// Now write out the MRH found results
			qp << "Query # " << i+1 << " MRH Found Results\n";
			qp << "Phi,Theta,P,Z,Rec\n";
			for( j=0; j < foundMRH.size(); j++) {
				qp << foundMRH[j].pt.phi << "," << foundMRH[j].pt.theta << ","
				   << foundMRH[j].pt.phi/pi << "," << cos(foundMRH[j].pt.theta) << ","
				   << foundMRH[j].rec << "\n";
			}
			qp << "\n";
			// Now write out the HPX found results
			qp << "Query # " << i+1 << " HPX Found Results\n";
			qp << "Phi,Theta,P,Z,Rec\n";
			for( j=0; j < foundHPX.size(); j++) {
				qp << foundHPX[j].pt.phi << "," << foundHPX[j].pt.theta << ","
				   << foundHPX[j].pt.phi/pi << "," << cos(foundHPX[j].pt.theta) << ","
				   << foundHPX[j].rec << "\n";	
			}
			qp << "\n\n";
			qp.close();
		}

	}
	cout << "	Finished Computing MRH vs HPX Poly Query Accuracy Test!\n";

	// Compute/get statistics of the final data structures 
	nMRHnodes = mMRH.NumNodes();
	nHPXnodes = mHPX.Npix();
    // Get Min,Max,Average Depth of QuadTrees of each data stucture
	maxDepthMRH = mMRH.MaxDepth();
	minDepthMRH = mMRH.MinDepth();
	avgDepthMRH = mMRH.AvgDepth();

	cout << "Number Points: " << ms.size() << endl;
	cout << "Number Poly Queries: " << polys.size() << endl;
	cout << "Total MRH Nodes: " << nMRHnodes << " Total HPX Nodes: " << nHPXnodes << endl;
	cout << "Memory Size of MRH: " << mMRH.SizeInMemKb() << " Kb" << " Size of HPX: " << (nHPXnodes*sizeof(Measurement))/1024 << " Kb" << endl;
	cout << "MRH Min Depth: " << minDepthMRH << " Max Depth: " << maxDepthMRH << " Avg Depth: " << avgDepthMRH << endl;
	cout << "MLQ\t#Nodes\tDepth\tAvgDepth\n";
	for( i = 0; i < 12; i++ ) {
		cout << i << "\t" 
			 << mMRH.NumNodesAtMLQ(i) << "\t"
			 << mMRH.DepthAtMLQ(i) << "\t"
			 << mMRH.AvgDepthAtMLQ(i) << "\n";
	}
	cout << "HPX Depth: " << order << endl;
	cout << "Total Matches: " << totalMatches << endl;
	cout << "Total MRH Found: " << totalMrhFound << " Total MRH Unique: " << totalMrhUnique << endl;
	cout << "Total HPX Found: " << totalHpxFound << " Total HPX Unique: " << totalHpxUnique << endl;
#ifdef CRITCOUNT
	cout << "Total MRH Crit Count: " << totalMrhCritCount << endl;
	cout << "Total HPX Crit Count: " << totalHpxCritCount << endl;
	cout << "Total MRH Upsearch Count: " << totalMrhUpsearchCount << endl;
#endif

	fp << "\nSummary:\n";
	fp << "Number Points: " << ms.size() << "\n";
	fp << "Number Poly Queries: " << polys.size() << "\n";
	fp << "Total MRH Nodes: " << nMRHnodes << " Total HPX Nodes: " << nHPXnodes << "\n";
	fp << "Memory Size of MRH: " << mMRH.SizeInMemKb() << " Kb" << " Size of HPX: " << (nHPXnodes*sizeof(Measurement))/1024 << " Kb\n";
	fp << "MRH Min Depth: " << minDepthMRH << " Max Depth: " << maxDepthMRH << " Avg Depth: " << avgDepthMRH << "\n";
	fp << "MLQ\t#Nodes\tDepth\tAvgDepth\n";
	for( i = 0; i < 12; i++ ) {
		fp << i << "\t" 
			 << mMRH.NumNodesAtMLQ(i) << "\t"
			 << mMRH.DepthAtMLQ(i) << "\t"
			 << mMRH.AvgDepthAtMLQ(i) << "\n";
	}
	fp << "HPX Depth: " << order << "\n";
	fp << "Total Matches: " << totalMatches << "\n";
	fp << "Total MRH Found: " << totalMrhFound << " Total MRH Unique: " << totalMrhUnique << "\n";
	fp << "Total HPX Found: " << totalHpxFound << " Total HPX Unique: " << totalHpxUnique << "\n";	
#ifdef CRITCOUNT
	fp << "Total MRH Crit Count: " << totalMrhCritCount << "\n";
	fp << "Total HPX Crit Count: " << totalHpxCritCount << "\n";
	fp << "Total MRH Upsearch Count: " << totalMrhUpsearchCount << "\n";
#endif

	cout << "\n##### END POLY QUERY RESULTS COMPARISON #####\n";
	fp << "\n##### END POLY QUERY RESULTS COMPARISON #####\n";

//###### END COMPARE POLY QUERY RESULTS ######
#endif




#ifdef DO_NEIGHBORS
//###### BEGIN COMPARE NEIGHBOR QUERY RESULTS ######

	cout << "\n\n##### BEGIN NEIGHBORS QUERY RESULTS COMPARISON #####\n";
	fp << "\n\n##### BEGIN NEIGHBORS QUERY RESULTS COMPARISON #####\n";

	cout << "Begin " << neighbors.size() << " Neighbor Queries..." << endl;
#ifdef CRITCOUNT
	fp << "Query#, Matches, MRH Only, HPX Only, MRH Crit, HPX Crit, MRH Upsearch" << "\n";
#else
	fp << "Query#, Matches, MRH Only, HPX Only," << "\n";
#endif
	totalMrhFound = 0;totalHpxFound=0;totalMatches=0;totalMrhUnique=0;totalHpxUnique=0;totalMrhCritCount=0;totalHpxCritCount=0;totalMrhUpsearchCount=0;
	for( i = 0; i < neighbors.size(); i++)
	{          


		// Open new output file containing data points, query definition
		// and found query points.
		if(verbose) {
			std::ostringstream stream;
			stream << outPath << "/neighborQuery" << i+1 << ".csv";
			std::string queryOut = stream.str();
			qp.open(queryOut.c_str());
			if(qp.fail()){
				cout << "Unable to open " << queryOut.c_str() << "!" << endl;
				exit(1);     
			}
			// First output the data points
			qp << "Phi,Theta,P,Z,Rec\n";
			for(j=0; j < ms.size(); j++) {
			   // Write out data location to query comparision file
			   qp << ms[j].pt.phi << "," << ms[j].pt.theta << ","
				  << ms[j].pt.phi/pi << "," << cos(ms[j].pt.theta) << ","
				  << ms[j].rec << "\n";
			}
			qp << "\n\n";

			// Write out the query definition
			qp << "Query #" << i+1 << "\n";
			qp << "Phi,Theta,P,Z,...\n";
			   qp << neighbors[i].pt.phi << "," << neighbors[i].pt.theta << ","
				  << neighbors[i].pt.phi/pi << "," << cos(neighbors[i].pt.theta) << "\n";
			qp << "\n";


		}

		Matches = 0; MRHUnique = 0; HPXUnique = 0;
		foundMRH.clear(); foundHPX.clear(); foundV.clear();
#ifdef CRITCOUNT
        HPX_CRIT_COUNT = 0;
		MRH_CRIT_COUNT = 0;
		MRH_UPSEARCH_COUNT = 0;
		MRH_COVER_MAP_CELL_RES = 0;
		mHPX.ResetCritCount();
		mMRH.ResetUpSearchCount();
#endif
		// Do MRH Neighbor Query
		foundMRH = mMRH.Neighbors(neighbors[i].pt,order);
#ifdef CRITCOUNT
		MRH_CRIT_COUNT = mMRH.GetCritCount();
		MRH_UPSEARCH_COUNT = mMRH.GetUpSearchCount();
		MRH_COVER_MAP_CELL_RES = mMRH.GetCoverMapCellRes();
#endif

		// Do HPX Neighbor Query
		hpxid = mHPX.ang2pix(neighbors[i].pt);
		mHPX.neighbors(hpxid,pixsetN);
#ifdef CRITCOUNT
		HPX_CRIT_COUNT = mHPX.GetCritCount();
#endif	

		// Spin through HPX's found indices keep only the
		// non empty indices that pass the point-in-disc test
		for( j=0; j < 8; j++ ) {
			if( mHPX[pixsetN[j]].rec != EMPTY ) {
				foundHPX.push_back(mHPX[pixsetN[j]]);
			} 
#ifdef CRITCOUNT
            HPX_CRIT_COUNT += 1;
#endif
		}

		AnalyzeResults(foundHPX,foundMRH,Matches,HPXUnique,MRHUnique,false);
		totalMatches += Matches;
		totalMrhUnique += MRHUnique;
		totalHpxUnique += HPXUnique;
		totalMrhFound += foundMRH.size();
		totalHpxFound += foundHPX.size();
#ifdef CRITCOUNT
		totalMrhCritCount += MRH_CRIT_COUNT;
		totalHpxCritCount += HPX_CRIT_COUNT;
		totalMrhUpsearchCount += MRH_UPSEARCH_COUNT;
#endif

#ifdef CRITCOUNT
		cout << "Query " << i+1 << " Matches: " << Matches << " MRH Only: " << MRHUnique << " HPX Only: " << HPXUnique 
			<< " MRH Crit Count: " << MRH_CRIT_COUNT << " HPX Crit Count: " << HPX_CRIT_COUNT << " MRH Upsearch Count: " << MRH_UPSEARCH_COUNT << endl;
		fp << i+1 << "," << Matches << "," << MRHUnique << "," << HPXUnique 
			<< "," << MRH_CRIT_COUNT << "," << HPX_CRIT_COUNT << "," << MRH_UPSEARCH_COUNT << "\n";
#else
		cout << "Query " << i+1 << " Matches: " << Matches << " MRH Only: " << MRHUnique << " HPX Only: " << HPXUnique << endl;
		
		// Output Trial Results in Table format
		fp << i+1 << "," << Matches << "," << MRHUnique << "," << HPXUnique << "\n";
#endif

		// Output Query and Query results to output file
		// Write Out Query # and the Query in phi,theta,p,z
		if(verbose) {
			// Now write out the MRH found results
			qp << "Query # " << i+1 << " MRH Found Results\n";
			qp << "Phi,Theta,P,Z,Rec\n";
			for( j=0; j < foundMRH.size(); j++) {
				qp << foundMRH[j].pt.phi << "," << foundMRH[j].pt.theta << ","
				   << foundMRH[j].pt.phi/pi << "," << cos(foundMRH[j].pt.theta) << ","
				   << foundMRH[j].rec << "\n";
			}
			qp << "\n";
			// Now write out the HPX found results
			qp << "Query # " << i+1 << " HPX Found Results\n";
			qp << "Phi,Theta,P,Z,Rec\n";
			for( j=0; j < foundHPX.size(); j++) {
				qp << foundHPX[j].pt.phi << "," << foundHPX[j].pt.theta << ","
				   << foundHPX[j].pt.phi/pi << "," << cos(foundHPX[j].pt.theta) << ","
				   << foundHPX[j].rec << "\n";	
			}
			qp << "\n\n";	
			qp.close();
		}
	}
	cout << "	Finished Computing MRH vs HPX Neighbor Query Accuracy Test!\n";

	// Compute/get statistics of the final data structures 
	nMRHnodes = mMRH.NumNodes();
	nHPXnodes = mHPX.Npix();
    // Get Min,Max,Average Depth of QuadTrees of each data stucture
	maxDepthMRH = mMRH.MaxDepth();
	minDepthMRH = mMRH.MinDepth();
	avgDepthMRH = mMRH.AvgDepth();

	cout << "Number Points: " << ms.size() << endl;
	cout << "Number Neighbor Queries: " << neighbors.size() << endl;
	cout << "Total MRH Nodes: " << nMRHnodes << " Total HPX Nodes: " << nHPXnodes << endl;
	cout << "Memory Size of MRH: " << mMRH.SizeInMemKb() << " Kb" << " Size of HPX: " << (nHPXnodes*sizeof(Measurement))/1024 << " Kb" << endl;
	cout << "MRH Min Depth: " << minDepthMRH << " Max Depth: " << maxDepthMRH << " Avg Depth: " << avgDepthMRH << endl;
	cout << "MLQ\t#Nodes\tDepth\tAvgDepth\n";
	for( i = 0; i < 12; i++ ) {
		cout << i << "\t" 
			 << mMRH.NumNodesAtMLQ(i) << "\t"
			 << mMRH.DepthAtMLQ(i) << "\t"
			 << mMRH.AvgDepthAtMLQ(i) << "\n";
	}
	cout << "HPX Depth: " << order << endl;
	cout << "Total Matches: " << totalMatches << endl;
	cout << "Total MRH Found: " << totalMrhFound << " Total MRH Unique: " << totalMrhUnique << endl;
	cout << "Total HPX Found: " << totalHpxFound << " Total HPX Unique: " << totalHpxUnique << endl;
#ifdef CRITCOUNT
	cout << "Total MRH Crit Count: " << totalMrhCritCount << endl;
	cout << "Total HPX Crit Count: " << totalHpxCritCount << endl;
	cout << "Total MRH Upsearch Count: " << totalMrhUpsearchCount << endl;
#endif

	fp << "\nSummary:\n";
	fp << "Number Points: " << ms.size() << "\n";
	fp << "Number Neighbor Queries: " << neighbors.size() << "\n";
	fp << "Total MRH Nodes: " << nMRHnodes << " Total HPX Nodes: " << nHPXnodes << "\n";
	fp << "Memory Size of MRH: " << mMRH.SizeInMemKb() << " Kb" << " Size of HPX: " << (nHPXnodes*sizeof(Measurement))/1024 << " Kb\n";
	fp << "MRH Min Depth: " << minDepthMRH << " Max Depth: " << maxDepthMRH << " Avg Depth: " << avgDepthMRH << "\n";
	fp << "MLQ\t#Nodes\tDepth\tAvgDepth\n";
	for( i = 0; i < 12; i++ ) {
		fp << i << "\t" 
			 << mMRH.NumNodesAtMLQ(i) << "\t"
			 << mMRH.DepthAtMLQ(i) << "\t"
			 << mMRH.AvgDepthAtMLQ(i) << "\n";
	}
	fp << "HPX Depth: " << order << "\n";
	fp << "Total Matches: " << totalMatches << "\n";
	fp << "Total MRH Found: " << totalMrhFound << " Total MRH Unique: " << totalMrhUnique << "\n";
	fp << "Total HPX Found: " << totalHpxFound << " Total HPX Unique: " << totalHpxUnique << "\n";
#ifdef CRITCOUNT
	fp << "Total MRH Crit Count: " << totalMrhCritCount << "\n";
	fp << "Total HPX Crit Count: " << totalHpxCritCount << "\n";
	fp << "Total MRH Upsearch Count: " << totalMrhUpsearchCount << "\n";
#endif

	cout << "\n##### END NEIGHBORS QUERY RESULTS COMPARISON #####\n";
	fp << "\n##### END NEIGHBORS QUERY RESULTS COMPARISON #####\n";

//###### END COMPARE NEIGHBOR QUERY RESULTS ######
#endif




#ifdef DO_STRIP
//###### BEGIN COMPARE STRIP QUERY RESULTS ######

	// Build HPX Data Structure
	cout << "\n\n##### BEGIN STRIP QUERY RESULTS COMPARISON #####\n";
	fp << "\n\n##### BEGIN STRIP QUERY RESULTS COMPARISON #####\n";

	cout << "Create new HPX Data Structure that's indexed using RING scheme...\n";
	fp << "Create new HPX Data Structure that's indexed using RING scheme...\n";

	cout << "Begin adding " << ms.size() << " records to the HPX Data Structure...\n";
	fp << "Begin adding " << ms.size() << " records to the HPX Data Structure...\n";
	Healpix_Map<Measurement> mHPX2(order,RING);
	for( i = 0; i < ms.size(); i++)
	{
		if( i%1000 == 0 ) {
			cout << "	Inserting Point " << i+1 << " of " << ms.size() << endl;
			fp << "	Inserting Point " << i+1 << " of " << ms.size() << "\n";
		}
		// Compute HEALPix index of next GIS longitude,latitude,data index tuple
		hpxid = mHPX2.ang2pix(ms[i].pt);

		// Insert data index associated with
		mHPX2[hpxid] = ms[i]; 
	}
	cout << "Completed adding " << ms.size() << " records to the HPX Data Structure!\n";
	fp << "Completed adding " << ms.size() << " records to the HPX Data Structure!\n";

	cout << "Begin " << strips.size() << " Strip Queries..." << endl;
#ifdef CRITCOUNT
	fp << "Query#, Matches, MRH Only, HPX Only, MRH Crit, HPX Crit, MRH Upsearch" << "\n";
#else
	fp << "Query#, Matches, MRH Only, HPX Only," << "\n";
#endif
	totalMrhFound = 0;totalHpxFound=0;totalMatches=0;totalMrhUnique=0;totalHpxUnique=0;totalMrhCritCount=0;totalHpxCritCount=0;totalMrhUpsearchCount=0;
	for( i = 0; i < strips.size(); i++)
	{          

		std::ostringstream stream;
		stream << outPath << "/stripQuery" << i+1 << ".csv";
		std::string queryOut = stream.str();
		qp.open(queryOut.c_str());
		if(qp.fail()){
			cout << "Unable to open " << queryOut.c_str() << "!" << endl;
			exit(1);     
		}
		// First output the data points
		qp << "Phi,Theta,P,Z,Rec\n";
		for(j=0; j < ms.size(); j++) {
		   // Write out data location to query comparision file
		   qp << ms[j].pt.phi << "," << ms[j].pt.theta << ","
			  << ms[j].pt.phi/pi << "," << cos(ms[j].pt.theta) << ","
			  << ms[j].rec << "\n";
		}
		qp << "\n\n";

		// Output Query and Query results to output file
		// Write Out Query # and the Query in phi,theta,p,z
		qp << "Query #" << i+1 << "\n";
		qp << "Theta1,Theta2,Z1,Z2\n";
		qp << strips[i].theta1 << "," << strips[i].theta2 << ","
		   << cos(strips[i].theta1) << "," << cos(strips[i].theta2) << "\n";
		qp << "\n\n";

		Matches = 0; MRHUnique = 0; HPXUnique = 0;
		foundMRH.clear(); foundHPX.clear(); foundV.clear();
#ifdef CRITCOUNT
        HPX_CRIT_COUNT = 0;
		MRH_CRIT_COUNT = 0;
		MRH_UPSEARCH_COUNT = 0;
		mHPX.ResetCritCount();
		mMRH.ResetUpSearchCount();
#endif

		// Do MRH Strip Query
		foundMRH = mMRH.QueryStrip(strips[i].theta1,strips[i].theta2);
#ifdef CRITCOUNT
		MRH_CRIT_COUNT = mMRH.GetCritCount();
		MRH_UPSEARCH_COUNT = mMRH.GetUpSearchCount();
#endif

		// Do HPX Strip Query
		mHPX2.query_strip(strips[i].theta1,strips[i].theta2,true,pixset);
#ifdef CRITCOUNT
		HPX_CRIT_COUNT = mHPX.GetCritCount();
#endif

		// Convert cell rangeset to list of individual
		// cells.
		pixset.toVector(foundV);

		// Spin through HPX's found indices keep only the
		// non empty indices that pass the point-in-disc test
		for( j=0; j < foundV.size(); j++ ) {
			if( mHPX2[foundV[j]].rec != EMPTY ) {
				if( IsPointInStrip(strips[i].theta1,strips[i].theta2,mHPX2[foundV[j]].pt) ){
					foundHPX.push_back(mHPX2[foundV[j]]);
				}
			} 
#ifdef CRITCOUNT
            HPX_CRIT_COUNT += 1;
#endif
		}
		AnalyzeResults(foundHPX,foundMRH,Matches,HPXUnique,MRHUnique,false);
		totalMatches += Matches;
		totalMrhUnique += MRHUnique;
		totalHpxUnique += HPXUnique;
		totalMrhFound += foundMRH.size();
		totalHpxFound += foundHPX.size();
#ifdef CRITCOUNT
		totalMrhCritCount += MRH_CRIT_COUNT;
		totalHpxCritCount += HPX_CRIT_COUNT;
		totalMrhUpsearchCount += MRH_UPSEARCH_COUNT;
#endif

#ifdef CRITCOUNT
		cout << "Query " << i+1 << " Matches: " << Matches << " MRH Only: " << MRHUnique << " HPX Only: " << HPXUnique 
			<< " MRH Crit Count: " << MRH_CRIT_COUNT << " HPX Crit Count: " << HPX_CRIT_COUNT << " MRH Upsearch Count: " << MRH_UPSEARCH_COUNT << endl;
		fp << i+1 << "," << Matches << "," << MRHUnique << "," << HPXUnique 
			<< "," << MRH_CRIT_COUNT << "," << HPX_CRIT_COUNT << "," << MRH_UPSEARCH_COUNT << "\n";
#else
		cout << "Query " << i+1 << " Matches: " << Matches << " MRH Only: " << MRHUnique << " HPX Only: " << HPXUnique << endl;
		
		// Output Trial Results in Table format
		fp << i+1 << "," << Matches << "," << MRHUnique << "," << HPXUnique << "\n";
#endif

		// Output Query and Query results to output file
		// Write Out Query # and the Query in phi,theta,p,z
		if(verbose) {
			// Now write out the MRH found results
			qp << "Query # " << i+1 << " MRH Found Results\n";
			qp << "Phi,Theta,P,Z,Rec\n";
			for( j=0; j < foundMRH.size(); j++) {
				qp << foundMRH[j].pt.phi << "," << foundMRH[j].pt.theta << ","
				   << foundMRH[j].pt.phi/pi << "," << cos(foundMRH[j].pt.theta) << ","
				   << foundMRH[j].rec << "\n";
			}
			qp << "\n\n";
			// Now write out the HPX found results
			qp << "Query # " << i+1 << " HPX Found Results\n";
			qp << "Phi,Theta,P,Z,Rec\n";
			for( j=0; j < foundHPX.size(); j++) {
				qp << foundHPX[j].pt.phi << "," << foundHPX[j].pt.theta << ","
				   << foundHPX[j].pt.phi/pi << "," << cos(foundHPX[j].pt.theta) << ","
				   << foundHPX[j].rec << "\n";	
			}
			qp << "\n\n";	
			qp.close();
		}
	}
	cout << "	Finished Computing MRH vs HPX Strip Query Accuracy Test!\n";

	// Compute/get statistics of the final data structures 
	nMRHnodes = mMRH.NumNodes();
	nHPXnodes = mHPX2.Npix();
    // Get Min,Max,Average Depth of QuadTrees of each data stucture
	maxDepthMRH = mMRH.MaxDepth();
	minDepthMRH = mMRH.MinDepth();
	avgDepthMRH = mMRH.AvgDepth();

	cout << "Number Points: " << ms.size() << endl;
	cout << "Number Strip Queries: " << strips.size() << endl;
	cout << "Total MRH Nodes: " << nMRHnodes << " Total HPX Nodes: " << nHPXnodes << endl;
	cout << "Memory Size of MRH: " << mMRH.SizeInMemKb() << " Kb" << " Size of HPX: " << (nHPXnodes*sizeof(Measurement))/1024 << " Kb" << endl;
	cout << "MRH Min Depth: " << minDepthMRH << " Max Depth: " << maxDepthMRH << " Avg Depth: " << avgDepthMRH << endl;
	cout << "MLQ\t#Nodes\tDepth\tAvgDepth\n";
	for( i = 0; i < 12; i++ ) {
		cout << i << "\t" 
			 << mMRH.NumNodesAtMLQ(i) << "\t"
			 << mMRH.DepthAtMLQ(i) << "\t"
			 << mMRH.AvgDepthAtMLQ(i) << "\n";
	}
	cout << "HPX Depth: " << order << endl;
	cout << "Total Matches: " << totalMatches << endl;
	cout << "Total MRH Found: " << totalMrhFound << " Total MRH Unique: " << totalMrhUnique << endl;
	cout << "Total HPX Found: " << totalHpxFound << " Total HPX Unique: " << totalHpxUnique << endl;
#ifdef CRITCOUNT
	cout << "Total MRH Crit Count: " << totalMrhCritCount << endl;
	cout << "Total HPX Crit Count: " << totalHpxCritCount << endl;
	cout << "Total MRH Upsearch Count: " << totalMrhUpsearchCount << endl;
#endif

	fp << "\nSummary:\n";
	fp << "Number Points: " << ms.size() << "\n";
	fp << "Number Strip Queries: " << strips.size() << "\n";
	fp << "Total MRH Nodes: " << nMRHnodes << " Total HPX Nodes: " << nHPXnodes << "\n";
	fp << "Memory Size of MRH: " << mMRH.SizeInMemKb() << " Kb" << " Size of HPX: " << (nHPXnodes*sizeof(Measurement))/1024 << " Kb\n";
	fp << "MRH Min Depth: " << minDepthMRH << " Max Depth: " << maxDepthMRH << " Avg Depth: " << avgDepthMRH << "\n";
	fp << "MLQ\t#Nodes\tDepth\tAvgDepth\n";
	for( i = 0; i < 12; i++ ) {
		fp << i << "\t" 
			 << mMRH.NumNodesAtMLQ(i) << "\t"
			 << mMRH.DepthAtMLQ(i) << "\t"
			 << mMRH.AvgDepthAtMLQ(i) << "\n";
	}
	fp << "HPX Depth: " << order << "\n";
	fp << "Total Matches: " << totalMatches << "\n";
	fp << "Total MRH Found: " << totalMrhFound << " Total MRH Unique: " << totalMrhUnique << "\n";
	fp << "Total HPX Found: " << totalHpxFound << " Total HPX Unique: " << totalHpxUnique << "\n";
#ifdef CRITCOUNT
	fp << "Total MRH Crit Count: " << totalMrhCritCount << "\n";
	fp << "Total HPX Crit Count: " << totalHpxCritCount << "\n";
	fp << "Total MRH Upsearch Count: " << totalMrhUpsearchCount << "\n";
#endif

	cout << "\n##### END STRIP QUERY RESULTS COMPARISON #####\n";
	fp << "\n##### END STRIP QUERY RESULTS COMPARISON #####\n";

//###### END COMPARE STRIP QUERY RESULTS ######
#endif

// END: RECORD RUNTIME FOR FULL BENCHMARK CODE
	stop_test = GetTimeMs64();

	cout << std::setprecision(20) << "\n####### TOTAL RUNTIME: " << double(stop_test-start_test) << " ms #######\n";
	fp << std::setprecision(20) << "\n####### TOTAL RUNTIME: " << double(stop_test-start_test) << " ms #######\n";


	fp.close();	

}


void BuildMRHFromArchive
(
 std::string archiveFile,
 std::string mapArchiveFile,
 std::string confirmFile,
 bool verbose,
 int maxDepthMRH
)

{
	int i = 0;
	// Create an MRH data structure constructed from
	// archive file!
	MultiResHpx_Map<Measurement> mMRH(maxDepthMRH,NEST);

	cout << "Now create MRH data structure from archive file!\n";

	cout << "Building MRH Quadtree Forest from archive...\n";
	mMRH.LoadFromFile(archiveFile);
	cout << "\tCOMPLETE!\n";

	// Now read in archived MRH Map array 
	cout << "Building MRH Map from archive...\n";
	mMRH.LoadMapFromArchive(mapArchiveFile);
	cout << "\tCOMPLETE!\n";

	// Output contents of MRH Map
	ofstream cf;
	cf.open(confirmFile.c_str());


	// Test integrity of MLQs
	////Print all trees
	for( i = 0; i < 12; i++ ) 
	{
		cout << "\n\n######################\n";
		cout << "#### BASE CELL " << i << " ####\n";
		cout << "######################" << "\n";
		mMRH.PrintTreeAtIndex(i); 
	}

	// Test integrity of MRH Map
	cout << "Size of recovered MRH Map: " << mMRH.NumRec() << endl;

	cf << "##### CONTENTS OF RECOVERED MRH MAP #####\n";
	for( i = 0; i < mMRH.NumRec(); i++)
	{
		cf << mMRH[i].data1 << "\t"
			<< mMRH[i].data2 << "\t"
			<< mMRH[i].data3 << "\t"
			<< mMRH[i].data4 << "\t"
			<< mMRH[i].pt.phi << "\t"
			<< mMRH[i].pt.theta << "\t"
			<< mMRH[i].val1 << "\t"
			<< mMRH[i].val2 << "\t"
			<< mMRH[i].val3 << "\t"
			<< mMRH[i].val4 << "\n";
	}
	cf.close();
}


void MRHBuildAndArchive
(
 std::string dataFile, 
 std::string archiveFile,
 int maxDepthMRH
)
{
	int i;
	int64 order;
	ifstream ip;
	std::vector<Measurement> ms;
	// Open and Parse Data File
	ip.open(dataFile.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << dataFile.c_str() << "!" << endl;
	  exit(1);     
	}
	ms = ReadRandomMeasurements(ip,order);
	ip.close();
	MultiResHpx_Map<Measurement> mMRH(maxDepthMRH,NEST);
	cout << "First load data points into MRH data structure...\n";
	for( i = 0; i < ms.size(); i++)
	{
		if( i%1000 == 0 ) {
			cout << "	Inserting Point " << i << " of " << ms.size() << endl;
		}
		mMRH.AddRecord(ms[i],ms[i].pt);
	}
	cout << endl << endl;

	// Archive the MRH Map 
	cout << "Archive the MRH Map...\n\n";
	std::string archiveMapFile;
	archiveMapFile = archiveFile + ".map";
	mMRH.SaveMapToArchive(archiveMapFile);
	cout << "MRH Map saved to archive file: " << archiveMapFile.c_str() << endl;

	// Archive the MRH MLQ forest 
	cout << "Archive the MRH MLQ forest...\n\n";
	mMRH.SaveToFile(archiveFile);
	cout << "MRH MLQ forest saved to archive file: " << archiveFile.c_str() << endl;

}




int main(int64 argc, char* argv[])

{
 	srand(time(NULL));

	int testNum, numTrials, verbose01,maxDepthMRH;
	bool verbose;
	std::string outFile,outFile2,outPath,dataFile,archiveFile,mapArchiveFile,queryFile;
	std::string confirmFile;
	std::vector<std::string> queryFilesList;

	testNum = atoi(argv[1]);

	cout << "Running Test Number: " << testNum << endl << endl;

	// Benchmark Test: HPX Insertion,Disc,Poly,Neighbors Query Benchmarks
	if( testNum == 1 ) {
		dataFile = argv[2];
		queryFilesList.push_back(argv[3]); // Disc Query File
		queryFilesList.push_back(argv[4]); // Poly Query File
		queryFilesList.push_back(argv[5]); // Neighbors Query File
		outFile = argv[6];
		numTrials = atoi(argv[7]);
		HPXBenchMarkFINAL(dataFile,queryFilesList,outFile,numTrials);
	}

	//Benchmark Test: HPX Strip Query Benchmark (because must be built with RING scheme).
	if( testNum == 2 ) {
		dataFile = argv[2];
		queryFile = argv[3];
		outFile = argv[4];
		numTrials = atoi(argv[5]);
		HPXStripBenchMarkFINAL(dataFile,queryFile,outFile,numTrials);
	}

	// Benchmark Test: MRH Insertion,Disc,Poly,Strip,Neighbors Query Benchmarks
	// AND
	// Benchmark Test: Compare Results of Disc,Poly,Strip and Neighbors Queries
	if( testNum == 3 ) {
		verbose = false;
		archiveFile = argv[2];  //MRH Archive File
		mapArchiveFile = argv[3];  //MRH Map Archive File
		queryFilesList.push_back(argv[4]); // Disc Query File
		queryFilesList.push_back(argv[5]); // Poly Query File
		queryFilesList.push_back(argv[6]); // Strip Query File
		queryFilesList.push_back(argv[7]); // Neighbors Query File
		outFile = argv[8];
		outPath = argv[9];
		numTrials = atoi(argv[10]);
		verbose01 = atoi(argv[11]);
		if(verbose01 == 1) verbose = true;
		maxDepthMRH = atoi(argv[12]);
		MRHBenchMarkFINAL(archiveFile,mapArchiveFile,queryFilesList,outFile,outFile2,outPath,numTrials,verbose,maxDepthMRH);
	}

	// Benchmark Test: Compare Results of Disc,Poly,Strip and Neighbors Queries
	if( testNum == 4 ) {
		verbose = false;
		archiveFile = argv[2];  //MRH Archive File
		mapArchiveFile = argv[3];  //MRH Map Archive File
		queryFilesList.push_back(argv[4]); // Measurement File
		queryFilesList.push_back(argv[5]); // Disc Query File
		queryFilesList.push_back(argv[6]); // Poly Query File
		queryFilesList.push_back(argv[7]); // Strip Query File
		queryFilesList.push_back(argv[8]); // Neighbors Query File
		outFile = argv[9];
		outPath = argv[10];
		verbose01 = atoi(argv[11]);
		maxDepthMRH = atoi(argv[12]);
		if(verbose01 == 1) verbose = true;
		MRHvsHPXComparisionFINAL(archiveFile,mapArchiveFile,queryFilesList,outFile,outPath,verbose,maxDepthMRH);
	}

	// Simply build MRH data structure from archive and report structure
	if( testNum == 0 ) {
		verbose = false;
		archiveFile = argv[2];  //MRH Archive File
		mapArchiveFile = argv[3];  //MRH Map Archive File
		confirmFile = argv[4];
		maxDepthMRH = atoi(argv[5]);
		verbose = atoi(argv[6]);
		if(verbose01 == 1) verbose = true;
		BuildMRHFromArchive(archiveFile,mapArchiveFile,confirmFile,verbose,maxDepthMRH);

	}

	// Build MRH Data Structure Then Archive it!
	if( testNum == -1 ) {
		dataFile = argv[2];
		outFile = argv[3];
		maxDepthMRH = atoi(argv[4]);
		MRHBuildAndArchive(dataFile,outFile,maxDepthMRH);
	}
}

