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



void buildRotArbAxis
( 
 double ang, 
 double ax, 
 double ay, 
 double az, 
 double (&mat)[3][3]
)
{
   // Construct Rotation Matrix
   double c = cos(D2R*ang);
   double s = sin(D2R*ang);
	   mat[0][0] = c + (1-c)*ax*ax;
	   mat[0][1] = (1-c)*ax*ay - s*az;
	   mat[0][2] = (1-c)*ax*az + s*ay;
	   mat[1][0] = (1-c)*ax*ay +s*az;
	   mat[1][1] = c + (1-c)*ay*ay;
	   mat[1][2] = (1-c)*ay*az - s*ax;
	   mat[2][0] = (1-c)*ax*az - s*ay;
	   mat[2][1] = (1-c)*ay*az + s*ax;
	   mat[2][2] = c + (1-c)*az*az;
}

void matrixMultiply
(
 double (mat)[3][3], 
 double x1, 
 double y1, 
 double z1, 
 double& x2, 
 double& y2,
 double& z2
)
{
   x2 = mat[0][0]*x1 + mat[1][0]*y1 + mat[2][0]*z1;
   y2 = mat[0][1]*x1 + mat[1][1]*y1 + mat[2][1]*z1;
   z2 = mat[0][2]*x1 + mat[1][2]*y1 + mat[2][2]*z1;
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

//######### PROCESS NOAA WEATHER STATION LOCATION FILE ########
void ProcessNOAA(std::string inputFile,std::string outputFile)
{
	// GIS Longitude (Phi) Range: [-180.0,+180.0] degrees, relative to Prime Meridian.
	// GIS Latitude (Theta) Range: [-90.0,+90.0] degrees, relative to Equator.
	// HPX Longitude (Phi) Range: [0.0,360.0] degrees, relative to Prime Meridian and progressing East.
	// HPX Colatitude (Theta) Range: [0.0,180.0] degrees, relative to North Pole and progressing South to
	// South Pole.
	int MIN_INT = -999; int MAX_INT = 999;
	double MIN_DOUBLE = -9.0; double MAX_DOUBLE = 9.0;
	int COOP,LATD,LATM,LATS,LONGD,LONGM,LONGS,ELEV,numRecs,order;
	int i,j,k,data1,data2,data3,data4;
	double val1,val2,val3,val4;
	int duplicates = 0;
	double nPhi,nTheta,minDist,curDist,pixres,temp_pixres;
	bool orderSet = false;
	bool duplicateYN;
	pointing ptA,ptB;
	std::vector<double> PHI,THETA;
	std::vector<int> COOPS;
	std::vector<int> LATDS,LATMS,LATSS,LONGDS,LONGMS,LONGSS,ELEVS;
	pointing pt1,pt2;
	ofstream ofp;
	ifstream ifp;

	// HPX Resolution Table
	std::vector<double> HpxResTable(30);
	HpxResTable[0] = 1.023326707946480; //order = 0
	HpxResTable[1] = 0.511663353973244; //order = 1
	HpxResTable[2] = 0.255831676986622; //order = 2
	HpxResTable[3] = 0.127915838493311; //order = 3
	HpxResTable[4] = 0.063957919246656; //order = 4
	HpxResTable[5] = 0.031978959623328; //order = 5
	HpxResTable[6] = 0.015989479811664; //order = 6
	HpxResTable[7] = 0.007994739905832; //order = 7
	HpxResTable[8] = 0.003997369952916; //order = 8
	HpxResTable[9] = 0.001998684976458; //order = 9
	HpxResTable[10] = 0.000999342488229; //order = 10
	HpxResTable[11] = 0.000499671244114; //order = 11
	HpxResTable[12] = 0.000249835622057; //order = 12
	HpxResTable[13] = 0.000124917811029; //order = 13
	HpxResTable[14] = 0.000062458905514; //order = 14
	HpxResTable[15] = 0.000031229452757; //order = 15
	HpxResTable[16] = 0.000015614726379; //order = 16
	HpxResTable[17] = 0.000007807363189; //order = 17
	HpxResTable[18] = 0.000003903681595; //order = 18
	HpxResTable[19] = 0.000001951840797; //order = 19
	HpxResTable[20] = 0.000000975920399; //order = 20
	HpxResTable[21] = 0.000000487960199; //order = 21
	HpxResTable[22] = 0.000000243980100; //order = 22
	HpxResTable[23] = 0.000000121990050; //order = 23
	HpxResTable[24] = 0.000000060995025; //order = 24
	HpxResTable[25] = 0.000000030497512; //order = 25
	HpxResTable[26] = 0.000000015248756; //order = 26
	HpxResTable[27] = 0.000000007624378; //order = 27
	HpxResTable[28] = 0.000000003812189; //order = 28
	HpxResTable[29] = 0.000000001906095; //order = 29

	ifp.open(inputFile.c_str());
	if(ifp.fail()){
	  cout << "Unable to open " << inputFile.c_str() << "!" << endl;
	  exit(1);     
	}

	// First line is number of records
	ifp >> numRecs;  
	skip_line(ifp); //skip rest of line

	// Skip header line
	skip_line(ifp);

	// Next lines are Measurement records
	duplicates = 0;
	cout << "Parsing NOAA file, checking for duplicate records..." << endl;
	for( i = 0; i < numRecs; i++)
	{
		ifp >> COOP >> LATD >> LATM >> LATS >> LONGD >> LONGM >> LONGS >> ELEV;
	 
 	    // Skip the rest of the line
		skip_line(ifp);

		// Check for duplicate Longitude & Latitude components 
		duplicateYN = false;
		for( j = 0; j < LONGDS.size(); j++ )
		{
			// Check for duplicate record, remove record if found.
			if(	LONGD == LONGDS[j] && LONGM == LONGMS[j] && LONGS == LONGSS[j] &&
				LATD == LATDS[j] && LATM == LATMS[j] && LATS == LATSS[j] ) 
			{
				duplicateYN = true;
				duplicates++;
			} 
		}
		if( duplicateYN == false ) {
			COOPS.push_back(COOP);
			LATDS.push_back(LATD);
			LATMS.push_back(LATM);
			LATSS.push_back(LATS);
			LONGDS.push_back(LONGD);
			LONGMS.push_back(LONGM);
			LONGSS.push_back(LONGS);
			ELEVS.push_back(ELEV);
		}
	}

	cout << "Duplicate records discovered: " << duplicates << endl;

	//Update numRecs
	numRecs = LONGDS.size();

	// Now must "scrub" the Longitude & Latitude pairs to eliminate duplicate
	// position data that would effect the required HPX resolution computation.
	cout << "Now converting GIS Longitude/Latitude to HPX Phi/Theta..." << endl;
	for( i = 0; i < numRecs; i++ )
	{
		// Convert D,M,S to fractional Degrees
		nPhi = double(LONGDS[i]);
		nPhi += (double(LONGMS[i])/60.0);
		nPhi += (double(LONGSS[i])/360.0);
       
		nTheta = double(LATDS[i]);
		nTheta += (double(LATMS[i])/60.0);
		nTheta += (double(LATSS[i])/360.0);
       
		// Convert Degrees to Radians
		nPhi *= D2R;
		nTheta *= D2R;

		// Convert GIS to HPX coordinate system
		pt1.phi = nPhi;
		pt1.theta = nTheta;
		pt2 = GIStoHPX(pt1);
       
		PHI.push_back(pt2.phi);
		THETA.push_back(pt2.theta);
	}

	ifp.close();

	// Now compute the minimum point separation to get proper
	// HPX resolution to use (set in output file). 
	minDist = 99999999.0;
	cout << "Now computing minimum data point separation ..." << endl;
	for( i = 0; i < numRecs; i++)
	{
	  ptA.phi = PHI[i];
	  ptA.theta = THETA[i];

	  for( j = 0; j < numRecs; j++)
	  {
		   if( i != j )
		   {
		      ptB.phi = PHI[j];
			  ptB.theta = THETA[j];
			  curDist = RadialDist(ptA,ptB);
	   	      
			  if( !AlmostEqual(curDist,0.0) ) 
			  {
				  if(curDist < minDist)
				  {
      				   minDist = curDist;
  					   // Determine proper HPX Order to set based on minDist vs HPX resolution
  					   orderSet = false;
					   for( k = 0; k <= 29; k++ ) 
					   {
							temp_pixres = HpxResTable[k];
							if(temp_pixres < minDist && orderSet == false)
							{
								order = k;
								pixres = temp_pixres;
								orderSet = true;
							}
					   }
					   cout << std::setprecision(21) << "New minDist: " << minDist;
					   cout << " giving HPX order: " << order;
					   cout << std::setprecision(21) << " of resolution " << pixres << endl << endl;
					   if( order > MAXHPXORDER ) {
							order = MAXHPXORDER;
							cout << "Exceeded maximum supported resolution for platform! Setting to order " << MAXHPXORDER << endl << endl;
					   }
				  }
			  }
		   }
	  }
	}
    cout << "Final minDist: " << minDist 
	     << " Use HPX Order: " << order
		 << " with Resolution: " << pixres << endl;	


	// Open the output file
	ofp.open(outputFile.c_str());

	if(ofp.fail()){
	  cout << "Unable to open " << outputFile.c_str() << "!" << endl;
	  exit(1);     
	}       
       
	// Write out "Measurement" file header
	ofp << numRecs << "\t" << order << endl;
	ofp << "REC\tPHI\tTHETA\tD1\tD2\tD3\tD4\tV1\tV2\tV3\tV4\t\n";

	// Write out "Measurement" record to output file
	cout << "Lastly write out converted NOAA records in Measurement file format..." << endl;
	for( i = 0; i < numRecs; i++ ) {
		data1 = ELEVS[i];
		data2 = COOPS[i];
		data3 = MIN_INT + static_cast <int> (rand()) /( static_cast <double> (RAND_MAX/(MAX_INT-MIN_INT)));
		data4 = MIN_INT + static_cast <int> (rand()) /( static_cast <double> (RAND_MAX/(MAX_INT-MIN_INT)));
		val1 = MIN_DOUBLE + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(MAX_DOUBLE-MIN_DOUBLE)));
		val2 = MIN_DOUBLE + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(MAX_DOUBLE-MIN_DOUBLE)));
		val3 = MIN_DOUBLE + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(MAX_DOUBLE-MIN_DOUBLE)));
		val4 = MIN_DOUBLE + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(MAX_DOUBLE-MIN_DOUBLE)));

		ofp << i << "\t"
			<< std::setprecision(21) << PHI[i] << "\t" << THETA[i] << "\t"
			<< std::setprecision(1) << data1 << "\t" << data2 << "\t" << data3 << "\t" << data4 << "\t"
			<< std::setprecision(5) << val1 << "\t" << val2 << "\t" << val3 << "\t" << val4 << "\n";
	}
	ofp.close();
	cout << "NOAA file processing complete!" << endl;
}



//######### PROCESS SDSS III OBJECT QUERY FILE ########
void ProcessSDSS(std::string inputFile,std::string outputFile)
{
	// (RA) Right Ascension Range: [0.0,360.00] degrees, relative to Prime Meridian and progressing East.
	// (DEC) Declination: [-90.0,+90] degrees, relative to Celestial Equator.
	// HPX Longitude (Phi) Range: [0.0,360.0] degrees, relative to Prime Meridian and progressing East.
	// HPX Colatitude (Theta) Range: [0.0,180.0] degrees, relative to North Pole and progressing South to
	// South Pole.
	int _run,_rerun,_camcol,_field,_obj,_type,numRecs,order,count;
	int i,j,k;
	double _ra,_dec,_z;
	int64 _objid;
	int duplicates = 0;
	double nPhi,nTheta,minDist,curDist,pixres,temp_pixres;
	bool orderSet = false;
	bool duplicateYN;
	pointing ptA,ptB;
	std::vector<double> PHI,THETA;
	std::vector<int64> OBJID;
	std::vector<double> Z;
	pointing pt1,pt2;
	ofstream ofp;
	ifstream ifp;

	// HPX Resolution Table
	std::vector<double> HpxResTable(30);
	HpxResTable[0] = 1.023326707946480; //order = 0
	HpxResTable[1] = 0.511663353973244; //order = 1
	HpxResTable[2] = 0.255831676986622; //order = 2
	HpxResTable[3] = 0.127915838493311; //order = 3
	HpxResTable[4] = 0.063957919246656; //order = 4
	HpxResTable[5] = 0.031978959623328; //order = 5
	HpxResTable[6] = 0.015989479811664; //order = 6
	HpxResTable[7] = 0.007994739905832; //order = 7
	HpxResTable[8] = 0.003997369952916; //order = 8
	HpxResTable[9] = 0.001998684976458; //order = 9
	HpxResTable[10] = 0.000999342488229; //order = 10
	HpxResTable[11] = 0.000499671244114; //order = 11
	HpxResTable[12] = 0.000249835622057; //order = 12
	HpxResTable[13] = 0.000124917811029; //order = 13
	HpxResTable[14] = 0.000062458905514; //order = 14
	HpxResTable[15] = 0.000031229452757; //order = 15
	HpxResTable[16] = 0.000015614726379; //order = 16
	HpxResTable[17] = 0.000007807363189; //order = 17
	HpxResTable[18] = 0.000003903681595; //order = 18
	HpxResTable[19] = 0.000001951840797; //order = 19
	HpxResTable[20] = 0.000000975920399; //order = 20
	HpxResTable[21] = 0.000000487960199; //order = 21
	HpxResTable[22] = 0.000000243980100; //order = 22
	HpxResTable[23] = 0.000000121990050; //order = 23
	HpxResTable[24] = 0.000000060995025; //order = 24
	HpxResTable[25] = 0.000000030497512; //order = 25
	HpxResTable[26] = 0.000000015248756; //order = 26
	HpxResTable[27] = 0.000000007624378; //order = 27
	HpxResTable[28] = 0.000000003812189; //order = 28
	HpxResTable[29] = 0.000000001906095; //order = 29

	ifp.open(inputFile.c_str());
	if(ifp.fail()){
	  cout << "Unable to open " << inputFile.c_str() << "!" << endl;
	  exit(1);     
	}

	// First line is number of records
	ifp >> numRecs;  
	skip_line(ifp); //skip rest of line

	// Skip header line
	skip_line(ifp);

	// Next lines are Measurement records
	duplicates = 0;
	cout << "Parsing SDSSIII Query file..." << endl;
	for( i = 0; i < numRecs; i++)
	{
		if(i%1000 == 0) {
			cout << "Processing record " << i << " of " << numRecs << endl;
		}
		ifp >> _ra >> _dec >> _z >> _objid;
	 
 	    // Skip the rest of the line
		skip_line(ifp);

		// Convert RA/DEC to HPX PHI/THETA
		pt1.phi = _ra*D2R;
		pt1.theta = _dec*D2R;
		pt2 = RADECtoHPX(pt1);

		// Store record
		PHI.push_back(pt2.phi);
		THETA.push_back(pt2.theta);
		OBJID.push_back(_objid);
		Z.push_back(_z);
	}
	ifp.close();

	// Now compute the minimum point separation to get proper
	// HPX resolution to use (set in output file). 
	minDist = 99999999.0;
	cout << "Now computing minimum data point separation ..." << endl;
	order = 0;
	for( i = 0; i < numRecs; i++)
	{
	  ptA.phi = PHI[i];
	  ptA.theta = THETA[i];

	  for( j = 0; j < numRecs; j++)
	  {
		   if( i != j )
		   {

				// Compute the minimum point seperaton if haven't reached maximum resolution
				if( order < MAXHPXORDER ) {
					ptB.phi = PHI[j];
					ptB.theta = THETA[j];
					curDist = RadialDist(ptA,ptB);
					if(curDist < minDist) {
						minDist = curDist;
						cout << "New MinDist: " << minDist << " Cur Order: " << order << endl;
					}
					// Compute required HEALPix resolution
					for( k = 0; k <= 29; k++ ) 
					{
						if( HpxResTable[k] < minDist )
						{
							order = k;
							break;
						}
					}
				} else {
					order = MAXHPXORDER;
					j = numRecs;
					i = numRecs;
				}
		   }
	  }
	}

	// Open the output file
	ofp.open(outputFile.c_str());

	if(ofp.fail()){
	  cout << "Unable to open " << outputFile.c_str() << "!" << endl;
	  exit(1);     
	}       
       
	// Write out "Measurement" file header
	ofp << numRecs << "\t" << order << endl;
	ofp << "REC\tPHI\tTHETA\tD1\tD2\tD3\tD4\tV1\tV2\tV3\tV4\t\n";

	// Write out "Measurement" record to output file
	cout << "Lastly write out converted SDSSIII Query records in Measurement file format..." << endl;
	for( i = 0; i < numRecs; i++ ) {

		ofp << i << "\t"
			<< std::setprecision(21) << PHI[i] << "\t" << THETA[i] << "\t"
			<< std::setprecision(1) << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t"
			<< std::setprecision(10) << Z[i] << "\t" << 0.0 << "\t" << 0.0 << "\t" << 0.0 << "\n";
	}
	ofp.close();
	cout << "SDSSIII Query file processing complete!" << endl;
}


//######### PROCESS FRAG FLYOUT FILE ########
void ProcessFRAG(std::string inputFile,std::string outputFile)
{

	double _time,_x,_y,_z,_vx,_vy,_vz,_mass;
	int _fragnum,numRecs,order;
	int i,j,k;
	int64 hpxid;
	double nPhi,nTheta,minDist,curDist,pixres,temp_pixres;
	bool orderSet = false;
	bool duplicateYN;
	pointing ptA,ptB;
	std::vector<double> PHI,THETA;
	std::vector<double> TIME,VX,VY,VZ,MASS;
	vec3 _vec;
	pointing pt1,pt2;
	ofstream ofp;
	ifstream ifp;
	Healpix_Custom mHPX(29,NEST);

	// HPX Resolution Table
	std::vector<double> HpxResTable(30);
	HpxResTable[0] = 1.023326707946480; //order = 0
	HpxResTable[1] = 0.511663353973244; //order = 1
	HpxResTable[2] = 0.255831676986622; //order = 2
	HpxResTable[3] = 0.127915838493311; //order = 3
	HpxResTable[4] = 0.063957919246656; //order = 4
	HpxResTable[5] = 0.031978959623328; //order = 5
	HpxResTable[6] = 0.015989479811664; //order = 6
	HpxResTable[7] = 0.007994739905832; //order = 7
	HpxResTable[8] = 0.003997369952916; //order = 8
	HpxResTable[9] = 0.001998684976458; //order = 9
	HpxResTable[10] = 0.000999342488229; //order = 10
	HpxResTable[11] = 0.000499671244114; //order = 11
	HpxResTable[12] = 0.000249835622057; //order = 12
	HpxResTable[13] = 0.000124917811029; //order = 13
	HpxResTable[14] = 0.000062458905514; //order = 14
	HpxResTable[15] = 0.000031229452757; //order = 15
	HpxResTable[16] = 0.000015614726379; //order = 16
	HpxResTable[17] = 0.000007807363189; //order = 17
	HpxResTable[18] = 0.000003903681595; //order = 18
	HpxResTable[19] = 0.000001951840797; //order = 19
	HpxResTable[20] = 0.000000975920399; //order = 20
	HpxResTable[21] = 0.000000487960199; //order = 21
	HpxResTable[22] = 0.000000243980100; //order = 22
	HpxResTable[23] = 0.000000121990050; //order = 23
	HpxResTable[24] = 0.000000060995025; //order = 24
	HpxResTable[25] = 0.000000030497512; //order = 25
	HpxResTable[26] = 0.000000015248756; //order = 26
	HpxResTable[27] = 0.000000007624378; //order = 27
	HpxResTable[28] = 0.000000003812189; //order = 28
	HpxResTable[29] = 0.000000001906095; //order = 29

	ifp.open(inputFile.c_str());
	if(ifp.fail()){
	  cout << "Unable to open " << inputFile.c_str() << "!" << endl;
	  exit(1);     
	}

	// First line is number of records
	ifp >> numRecs;  
	skip_line(ifp); //skip rest of line

	// Skip header line
	skip_line(ifp);

	// Next lines are Measurement records
	cout << "Parsing Frag Flyout file, checking for duplicate records..." << endl;
	for( i = 0; i < numRecs; i++)
	{
		ifp >> _fragnum >> _time >> _x >> _y >> _z >> _vx >> _vy >> _vz >> _mass;
	 
 	    // Skip the rest of the line
		skip_line(ifp);

		// Convert velocity vector to HPX PHI/THETA
		_vec.x = _vx;
		_vec.y = _vy;
		_vec.z = _vz;
		hpxid = mHPX.vec2pix(_vec);
		pt1 = mHPX.pix2ang(hpxid);

		// Store record
		PHI.push_back(pt1.phi);
		THETA.push_back(pt1.theta);
		MASS.push_back(_mass);

	}
	ifp.close();

	// Now compute the minimum point separation to get proper
	// HPX resolution to use (set in output file). 
	minDist = 99999999.0;
	cout << "Now computing minimum data point separation ..." << endl;
	for( i = 0; i < numRecs; i++)
	{
	  ptA.phi = PHI[i];
	  ptA.theta = THETA[i];

	  for( j = 0; j < numRecs; j++)
	  {
		   if( i != j )
		   {
		      ptB.phi = PHI[j];
			  ptB.theta = THETA[j];
			  curDist = RadialDist(ptA,ptB);
	   	      
			  if( !AlmostEqual(curDist,0.0) ) 
			  {
				  if(curDist < minDist)
				  {
      				   minDist = curDist;
  					   // Determine proper HPX Order to set based on minDist vs HPX resolution
  					   orderSet = false;
					   for( k = 0; k <= 29; k++ ) 
					   {
							temp_pixres = HpxResTable[k];
							if(temp_pixres < minDist && orderSet == false)
							{
								order = k;
								pixres = temp_pixres;
								orderSet = true;
							}
					   }
					   cout << std::setprecision(21) << "New minDist: " << minDist;
					   cout << " giving HPX order: " << order;
					   cout << std::setprecision(21) << " of resolution " << pixres << endl << endl;
					   if( order > MAXHPXORDER ) {
							order = MAXHPXORDER;
							cout << "Exceeded maximum supported resolution for platform! Setting to order " << MAXHPXORDER << endl << endl;
					   }
				  }
			  }
		   }
	  }
	}
    cout << "Final minDist: " << minDist 
	     << " Use HPX Order: " << order
		 << " with Resolution: " << pixres << endl;	

	// Open the output file
	ofp.open(outputFile.c_str());

	if(ofp.fail()){
	  cout << "Unable to open " << outputFile.c_str() << "!" << endl;
	  exit(1);     
	}       
       
	// Write out "Measurement" file header
	ofp << numRecs << "\t" << order << endl;
	ofp << "REC\tPHI\tTHETA\tD1\tD2\tD3\tD4\tV1\tV2\tV3\tV4\t\n";

	// Write out "Measurement" record to output file
	cout << "Lastly write out converted Frag Flyout records in Measurement file format..." << endl;
	for( i = 0; i < numRecs; i++ ) {

		ofp << i << "\t"
			<< std::setprecision(21) << PHI[i] << "\t" << THETA[i] << "\t" 
			<< std::setprecision(1) << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t"
			<< std::setprecision(10) << MASS[i] << "\t" << 0.0 << "\t" << 0.0 << "\t" << 0.0 << "\n";
	}
	ofp.close();
	cout << "Frag Flyout file processing complete!" << endl;
}

//######### PROCESS MOON CRATER LOCATION FILE ########
void ProcessMOON(std::string inputFile,std::string outputFile)
{
	// LUNAR Longitude (Phi) Range: [-180.0,+180.0] degrees, relative to Prime Meridian.
	// LUNAR Latitude (Theta) Range: [-90.0,+90.0] degrees, relative to Equator.
	// HPX Longitude (Phi) Range: [0.0,360.0] degrees, relative to Prime Meridian and progressing East.
	// HPX Colatitude (Theta) Range: [0.0,180.0] degrees, relative to North Pole and progressing South to
	// South Pole.
	int MIN_INT = -999; int MAX_INT = 999;
	double MIN_DOUBLE = -9.0; double MAX_DOUBLE = 9.0;
	double LONG,LAT,R,D;
	int numRecs,order,DR,P;
	int i,j,k,data1,data2,data3,data4;
	double val1,val2,val3,val4;
	int duplicates = 0;
	double nPhi,nTheta,minDist,curDist,pixres,temp_pixres;
	bool orderSet = false;
	bool duplicateYN;
	pointing ptA,ptB;
	std::vector<double> PHI,THETA;
	std::vector<double> LONGS,LATS,RS,DS;
	std::vector<int> DRS,PS;
	pointing pt1,pt2;
	ofstream ofp;
	ifstream ifp;

	// HPX Resolution Table
	std::vector<double> HpxResTable(30);
	HpxResTable[0] = 1.023326707946480; //order = 0
	HpxResTable[1] = 0.511663353973244; //order = 1
	HpxResTable[2] = 0.255831676986622; //order = 2
	HpxResTable[3] = 0.127915838493311; //order = 3
	HpxResTable[4] = 0.063957919246656; //order = 4
	HpxResTable[5] = 0.031978959623328; //order = 5
	HpxResTable[6] = 0.015989479811664; //order = 6
	HpxResTable[7] = 0.007994739905832; //order = 7
	HpxResTable[8] = 0.003997369952916; //order = 8
	HpxResTable[9] = 0.001998684976458; //order = 9
	HpxResTable[10] = 0.000999342488229; //order = 10
	HpxResTable[11] = 0.000499671244114; //order = 11
	HpxResTable[12] = 0.000249835622057; //order = 12
	HpxResTable[13] = 0.000124917811029; //order = 13
	HpxResTable[14] = 0.000062458905514; //order = 14
	HpxResTable[15] = 0.000031229452757; //order = 15
	HpxResTable[16] = 0.000015614726379; //order = 16
	HpxResTable[17] = 0.000007807363189; //order = 17
	HpxResTable[18] = 0.000003903681595; //order = 18
	HpxResTable[19] = 0.000001951840797; //order = 19
	HpxResTable[20] = 0.000000975920399; //order = 20
	HpxResTable[21] = 0.000000487960199; //order = 21
	HpxResTable[22] = 0.000000243980100; //order = 22
	HpxResTable[23] = 0.000000121990050; //order = 23
	HpxResTable[24] = 0.000000060995025; //order = 24
	HpxResTable[25] = 0.000000030497512; //order = 25
	HpxResTable[26] = 0.000000015248756; //order = 26
	HpxResTable[27] = 0.000000007624378; //order = 27
	HpxResTable[28] = 0.000000003812189; //order = 28
	HpxResTable[29] = 0.000000001906095; //order = 29

	ifp.open(inputFile.c_str());
	if(ifp.fail()){
	  cout << "Unable to open " << inputFile.c_str() << "!" << endl;
	  exit(1);     
	}

	// First line is number of records
	ifp >> numRecs;  
	skip_line(ifp); //skip rest of line

	// Skip header line
	skip_line(ifp);

	// Next lines are Measurement records
	cout << "Parsing Moon Crater Positions file..." << endl;
	for( i = 0; i < numRecs; i++)
	{
		ifp >> LONG >> LAT >> R >> D >> DR >> P;
	 
 	    // Skip the rest of the line
		skip_line(ifp);

		LONGS.push_back(LONG);
		LATS.push_back(LAT);
		RS.push_back(R);
		DS.push_back(D);
		DRS.push_back(DR);
		PS.push_back(P);
	}


	cout << "Now converting GIS Longitude/Latitude to HPX Phi/Theta..." << endl;
	for( i = 0; i < numRecs; i++ )
	{      
		// Convert Degrees to Radians
		nPhi = D2R*LONGS[i];
		nTheta = D2R*LATS[i];

		// Convert GIS to HPX coordinate system
		pt1.phi = nPhi;
		pt1.theta = nTheta;
		pt2 = GIStoHPX(pt1);
       
		PHI.push_back(pt2.phi);
		THETA.push_back(pt2.theta);
	}

	ifp.close();

	// Now compute the minimum point separation to get proper
	// HPX resolution to use (set in output file). 
	minDist = 99999999.0;
	cout << "Now computing minimum data point separation ..." << endl;
	order = 0;
	for( i = 0; i < numRecs; i++)
	{
	  ptA.phi = PHI[i];
	  ptA.theta = THETA[i];

	  for( j = 0; j < numRecs; j++)
	  {
		   if( i != j )
		   {

				// Compute the minimum point seperaton if haven't reached maximum resolution
				if( order < MAXHPXORDER ) {
					ptB.phi = PHI[j];
					ptB.theta = THETA[j];
					curDist = RadialDist(ptA,ptB);
					if(curDist < minDist) {
						minDist = curDist;
						cout << "New MinDist: " << minDist << " Cur Order: " << order << endl;
					}
					// Compute required HEALPix resolution
					for( k = 0; k <= 29; k++ ) 
					{
						if( HpxResTable[k] < minDist )
						{
							order = k;
							break;
						}
					}
				} else {
					order = MAXHPXORDER;
					j = numRecs;
					i = numRecs;
				}
		   }
	  }
	}

	// Open the output file
	ofp.open(outputFile.c_str());

	if(ofp.fail()){
	  cout << "Unable to open " << outputFile.c_str() << "!" << endl;
	  exit(1);     
	}       
       
	// Write out "Measurement" file header
	ofp << numRecs << "\t" << order << endl;
	ofp << "REC\tPHI\tTHETA\tD1\tD2\tD3\tD4\tV1\tV2\tV3\tV4\t\n";

	// Write out "Measurement" record to output file
	cout << "Lastly write out converted Moon Crater Position records in Measurement file format..." << endl;
	for( i = 0; i < numRecs; i++ ) {
		data1 = DRS[i];
		data2 = PS[i];
		data3 = 0;
		data4 = 0;
		val1 = RS[i];
		val2 = DS[i];
		val3 = 0.0;
		val4 = 0.0;

		ofp << i << "\t"
			<< std::setprecision(21) << PHI[i] << "\t" << THETA[i] << "\t"
			<< std::setprecision(1) << data1 << "\t" << data2 << "\t" << data3 << "\t" << data4 << "\t"
			<< std::setprecision(10) << val1 << "\t" << val2 << "\t" << val3 << "\t" << val4 << "\n";
	}
	ofp.close();
	cout << "Moon Crater Position file processing complete!" << endl;
}



//######### WRITE & READ RANDOM MEASUREMENT DATA TO FILE #########

void WriteRandomMeasurementsHeader(ofstream& fp,int64 numRec,int64 hpxOrder)
{     
	fp << numRec << "\t" << hpxOrder << '\n';
	fp << "REC\t";
	fp << "PHI\t";
	fp << "THETA\t";
	fp << "D1\t";
	fp << "D2\t";
	fp << "D3\t";
	fp << "D4\t";
	fp << "V1\t";
	fp << "V2\t";
	fp << "V3\t";
	fp << "V4\t\n";
}

void WriteRandomMeasurements(ofstream& fp, std::vector<Measurement> m)
{
	unsigned int i;
	for( i = 0; i < m.size(); i++)
	{
		fp << std::setprecision(32) << m[i].rec << "\t";
		fp << std::setprecision(32) << m[i].pt.phi << "\t";
		fp << std::setprecision(32) << m[i].pt.theta << "\t";
		fp << m[i].data1 << "\t";
		fp << m[i].data2 << "\t";
		fp << m[i].data3 << "\t";
		fp << m[i].data4 << "\t";
		fp << m[i].val1 << "\t";
		fp << m[i].val2 << "\t";
		fp << m[i].val3 << "\t";
		fp << m[i].val4 << "\t\n";
	}
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


//######### WRITE RANDOM DISC DEFINITIONS TO FILE #########

void WriteRandomDiscHeader(ofstream& fp,int numRec)
{     
	fp << numRec << '\n';
	fp << "PHI\t";
	fp << "THETA\t";
	fp << "RADIUS\t\n";
}

void WriteRandomDiscs(ofstream& fp, std::vector<DiscType> d)
{
	unsigned int i;
	for( i = 0; i < d.size(); i++)
	{
		fp << d[i].pt.phi << "\t";
		fp << d[i].pt.theta << "\t";
		fp << d[i].radius << "\t\n";
	}
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

//######### WRITE RANDOM POLYGON DEFINITIONS TO FILE #########

void WriteRandomPolysHeader(ofstream& fp,int numRec)
{     
	fp << numRec << '\n';
	fp << "NUMPTS" << "\t";
	fp << "PHI\t";
	fp << "THETA\t";
	fp << "ETC...\t\n";
}

void WriteRandomPolys(ofstream& fp, std::vector<PolyType> p)
{
	unsigned int i,j;
	for( i = 0; i < p.size(); i++)
	{
		fp << p[i].pts.size() << "\t";
		for( j = 0; j < p[i].pts.size(); j++) {
			fp << p[i].pts[j].phi << "\t";
			fp << p[i].pts[j].theta << "\t";
		}
		fp << "\n";
	}
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

//######### WRITE RANDOM THETA STRIP DEFINITIONS TO FILE #########

void WriteRandomStripsHeader(ofstream& fp,int numRec)
{     
	fp << numRec << '\n';
	fp << "THETA1\t";
	fp << "THETA2\t\n";
}

void WriteRandomStrips(ofstream& fp, std::vector<StripType> s)
{
	unsigned int i;
	for( i = 0; i < s.size(); i++)
	{
		fp << s[i].theta1 << "\t";
		fp << s[i].theta2 << "\t\n";
	}
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


//######### WRITE RANDOM NEIGHBOR INDICES TO FILE #########

void WriteRandomNeighborsHeader(ofstream& fp,int numRec)
{     
	fp << numRec << '\n';
	fp << "PHI\tTHETA\t\n";
}

void WriteRandomNeighbors(ofstream& fp, std::vector<NeighborType> n)
{
	unsigned int i;
	for( i = 0; i < n.size(); i++)
	{
		fp << n[i].pt.phi << "\t" << n[i].pt.theta << "\t\n";
	}
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


//######### GENERATE RANDOM POINT LOCATIONS #########

void GenRandomPhiTheta
(
 int numPoints,
 std::vector<pointing>& points,
 int64& order,
 bool verbose,
 double minPhi,
 double maxPhi,
 double minTheta,
 double maxTheta
)
{
	// GIS Longitude (Phi) Range: [-180.0,+180.0] degrees, relative to Prime Meridian.
	// GIS Latitude (Theta) Range: [-90.0,+90.0] degrees, relative to Equator.
	// HPX Longitude (Phi) Range: [0.0,360.0] degrees, relative to Prime Meridian and progressing East.
	// HPX Colatitude (Theta) Range: [0.0,180.0] degrees, relative to North Pole and progressing South to
	// South Pole.
	pointing p1,p2;
	double minRad,curRad;
	unsigned int i,j;
	int64 nside;
	double pixres;
	int count = 0;
	std::vector<double> HpxOrderResTable;
	points.clear();
	order = -1;
	minRad = 9999999999999.0;


	// Generate HEALPix Stats Tables
	for( i = 0; i <= 29; i++ ) 
	{
	  order = i;
      nside  = int64(1)<<order;
	  pixres = sqrt(pi/(3.0*nside*nside));
	  HpxOrderResTable.push_back(pixres);
	  if(verbose) {
		printf("HPX order = %d HPX cell res = %2.15f\n",order,pixres);
	  }
	}

	// Generate random HPX longitude, latitude pairs
	if(verbose) {
		cout << "#,PHI,THETA,P,Z" << endl;
	}

	order = -1;
	while( count < numPoints )
	{
		p1.phi = minPhi + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(maxPhi-minPhi)));
		p1.theta = minTheta + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(maxTheta-minTheta)));

		// Compute the minimum point seperaton if haven't reached maximum resolution
		if( order != MAXHPXORDER ) {
			for( j = 0; j < points.size(); j++ ) {
				p2 = points[j];
				curRad = acos(fabs(cosdist_zphi(cos(p1.theta),p1.phi,cos(p2.theta),p2.phi)));
				if(curRad < minRad) {
					minRad = curRad;
				}
			}
			// Compute required HEALPix resolution
			for( i = 0; i <= 29; i++ ) 
			{
				if( HpxOrderResTable[i] < minRad )
				{
					order = i;
					break;
				}
			}
		}
		points.push_back(p1);
		if(verbose) {
			cout << count+1 << "," 
				 << p1.phi << "," 
				 << p1.theta << "," 
				 << p1.phi/pi << "," 
				 << cos(p1.theta) << endl;
		}
		count++;
		if(count%1000 == 0) {
			cout << "Generated point number " << count << " of " << numPoints << endl;
			if(order == MAXHPXORDER) {
				cout << "	Reached maximum HEALPix resolution of order: " << order << endl;
			} else {
				cout << "	Current required HEALPix resolution of order: " << order << endl;
			}
		}
	}

	if(verbose) {
		printf("For random dataset computed: HPX order = %d HPX cell res = %2.15f\n\n",order,minRad);
	}

}

//######### GENERATE RANDOM MEASUREMENTS #########

void GenRandomData
(
 int numPoints,
 std::vector<Measurement>& ms,
 double minPhi,
 double maxPhi,
 double minTheta,
 double maxTheta
)
{
	int MIN_INT = -999; int MAX_INT = 999;
	double MIN_DOUBLE = -9.0; double MAX_DOUBLE = 9.0;
	Measurement m;
	std::vector<pointing> pts;
	int64 order;
	ms.clear();
	int64 recNum = 0;
	for(int i = 0; i < numPoints; i++) {
		GenRandomPhiTheta(1,pts,order,false,minPhi,maxPhi,minTheta,maxTheta);
		m.rec = recNum;
		m.pt.phi = pts[0].phi;
		m.pt.theta = pts[0].theta;
		m.data1 = MIN_INT + static_cast <int> (rand()) /( static_cast <double> (RAND_MAX/(MAX_INT-MIN_INT)));
		m.data2 = MIN_INT + static_cast <int> (rand()) /( static_cast <double> (RAND_MAX/(MAX_INT-MIN_INT)));
		m.data3 = MIN_INT + static_cast <int> (rand()) /( static_cast <double> (RAND_MAX/(MAX_INT-MIN_INT)));
		m.data4 = MIN_INT + static_cast <int> (rand()) /( static_cast <double> (RAND_MAX/(MAX_INT-MIN_INT)));

		m.val1 = MIN_DOUBLE + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(MAX_DOUBLE-MIN_DOUBLE)));
		m.val2 = MIN_DOUBLE + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(MAX_DOUBLE-MIN_DOUBLE)));
		m.val3 = MIN_DOUBLE + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(MAX_DOUBLE-MIN_DOUBLE)));
		m.val4 = MIN_DOUBLE + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(MAX_DOUBLE-MIN_DOUBLE)));
		
		ms.push_back(m);		
		recNum += 1;
	}
}


//######### GENERATE RANDOM QUERY DISCS #########

void RandomDisc
(
 double& phi, 
 double& theta, 
 double& radius,
 double minPhi, 
 double maxPhi,
 double minTheta,
 double maxTheta,
 double minRad,
 double maxRad
 ){
	// GIS Longitude (Phi) Range: [-180.0,+180.0] degrees, relative to Prime Meridian.
	// GIS Latitude (Theta) Range: [-90.0,+90.0] degrees, relative to Equator.
	// HPX Longitude (Phi) Range: [0.0,360.0] degrees, relative to Prime Meridian and progressing East.
	// HPX Colatitude (Theta) Range: [0.0,180.0] degrees, relative to North Pole and progressing South to
	// South Pole.
	pointing pt;
	bool validDisc = false;

	while( validDisc == false ) {
		// Generate random query center location in HPX Coordinate System
		phi = minPhi + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(maxPhi-minPhi)));
		theta = minTheta + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(maxTheta-minTheta)));

		// Generate random query radius, making sure radius doesn't sweep
		// off of allowed latitude,longitude limits!
		radius = minRad + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(maxRad-minRad)));

		if( phi+radius < maxPhi &&
			phi-radius > minPhi &&
			theta+radius < maxTheta &&
			theta-radius > minTheta ) {
			validDisc = true;
		}
	}
}

//######### GENERATE RANDOM QUERY POLYGONS #########

vector<pointing>  RandomConvexPoly
(
 bool RANDPOLYVERBOSE,
 double minPhi, 
 double maxPhi,
 double minTheta,
 double maxTheta,
 double minRad,
 double maxRad
) 
{
	// GIS Longitude (Phi) Range: [-180.0,+180.0] degrees, relative to Prime Meridian.
	// GIS Latitude (Theta) Range: [-90.0,+90.0] degrees, relative to Equator.
	// HPX Longitude (Phi) Range: [0.0,360.0] degrees, relative to Prime Meridian and progressing East.
	// HPX Colatitude (Theta) Range: [0.0,180.0] degrees, relative to North Pole and progressing South to
	// South Pole.
	int MIN_PTS = 3; 
	double MIN_ARC = 5.0*D2R; double MAX_ARC,nArc;
	double ttlArc;
	double qPhi,qTheta,qRadius;
	pointing pt;
	vector<pointing> poly;

	// First get a valid query disc (fits in the bounds of HPX coordinate system).
	RandomDisc(qPhi,qTheta,qRadius,minPhi,maxPhi,minTheta,maxTheta,minRad,maxRad);
	
	if(RANDPOLYVERBOSE) {
		cout << qPhi << "," << qTheta << "," << qRadius << endl;
	}

	// Generate random maximum arc
	MAX_ARC = MIN_ARC + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(180.0*D2R-MIN_ARC)));

	// Generate an N sided, clockwise, around query disc
	// convex polyton. Basically keep adding sides to the polyton until
	// arc excedes 2PI (full circle). Coordinates of poly are in HPX.
	ttlArc = 0.0;
	poly.clear();
	while( ttlArc < TWOPI ) {
		// Compute points on disc
		pt.phi   = qPhi + qRadius * cos(ttlArc);
		pt.theta = qTheta + qRadius * sin(ttlArc);
		poly.push_back(pt);

		// Compute next arc length
		nArc = MIN_ARC + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(MAX_ARC-MIN_ARC)));
		ttlArc += nArc;
	}
	if(RANDPOLYVERBOSE) {
		cout << endl << "### CW, Concave Polygon Definition ###\n";
		cout << "POINT,PHI,THETA,P,Z\n";
		for(unsigned int j = 0; j < poly.size(); j++ ) {
			cout << j << "," << poly[j].phi << "," << poly[j].theta << "," 
				<< poly[j].phi/pi << "," << cos(poly[j].theta) << endl;		
		}	
		cout << 0 << "," << poly[0].phi << "," << poly[0].theta << "," 
			<< poly[0].phi/pi << "," << cos(poly[0].theta) << endl;
	}

	return poly;
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


//####################################################################################
//####################################################################################
//######																		######
//######   MRH AND HPX DATA POINT INSERTION AND MEMORY FOOT PRINT BENCHMARKING	######
//######																		######
//####################################################################################
//####################################################################################


// Test to make sure that MRH will over-write previously loaded records and spatial 
// references in MortonLQT for Benchmarking purposes.
void MRHDuplicatePointInsertionTest
(
 int maxDepth,
 int numPoints,
 int duplicateYN,
 double minPhi,
 double maxPhi,
 double minTheta,
 double maxTheta
)
{
	vector<pointing> points;
	vector<Measurement> measurements;
	int64 order;
	int i;
	points.clear();
	measurements.clear();
	MultiResHpx_Map<Measurement> mMRH(maxDepth,NEST);

	// Generate ONE random GIS longitude, latitude pairs
	if(duplicateYN) {
		GenRandomPhiTheta(1,points,order,true,minPhi,maxPhi,minTheta,maxTheta);
	} else {
		GenRandomPhiTheta(numPoints,points,order,true,minPhi,maxPhi,minTheta,maxTheta);
	}

	// Generate random measurement data
	GenRandomData(numPoints,measurements,minPhi,maxPhi,minTheta,maxTheta);

	for( i = 0; i < numPoints; i++)
	{
		// Insert next GIS longitude,latitude, data index tuple into MRH
		if(duplicateYN) {
			measurements[i].pt = points[0]; // Set DUPLICATE spatial location of measurement
			mMRH.AddRecord(measurements[i],points[0]);
		} else {
			measurements[i].pt = points[i]; 
			cout << "Next Measurement To Insert:" << endl;
			cout << "	Phi: " << measurements[i].pt.phi << endl;
			cout << "	Theta: " << measurements[i].pt.theta << endl;
			cout << "	Rec: " << measurements[i].rec << endl;
			cout << "	Data 1: " << measurements[i].data1 << endl;
			cout << "	Data 2: " << measurements[i].data2 << endl;
			cout << "	Data 3: " << measurements[i].data3 << endl;
			cout << "	Data 4: " << measurements[i].data4 << endl;
			cout << "	Val 1: " << measurements[i].val1 << endl;
			cout << "	Val 2: " << measurements[i].val2 << endl;
			cout << "	Val 3: " << measurements[i].val3 << endl;
			cout << "	Val 4: " << measurements[i].val4 << endl;
			cout << endl;
			mMRH.AddRecord(measurements[i],points[i]);
		}
	}

	//Print all trees
	for( i = 0; i < 12; i++ ) 
	{
		cout << "\n\n######################\n";
		cout << "#### BASE CELL " << i << " ####\n";
		cout << "######################" << endl;
		mMRH.PrintTreeAtIndex(i); 
	}

	//Dump the MRH Map
	vector<Measurement> mMap = mMRH.Map();
	cout << "\n\n#################\n";
	cout << "#### MRH MAP ####\n";
	cout << "#################" << endl;
	for( i = 0; i < mMap.size(); i++ ) {
		cout << "Measurement at Map Location: " << i << endl;
		cout << "	Phi: " << mMap[i].pt.phi << endl;
		cout << "	Theta: " << mMap[i].pt.theta << endl;
		cout << "	Rec: " << mMap[i].rec << endl;
		cout << "	Data 1: " << mMap[i].data1 << endl;
		cout << "	Data 2: " << mMap[i].data2 << endl;
		cout << "	Data 3: " << mMap[i].data3 << endl;
		cout << "	Data 4: " << mMap[i].data4 << endl;
		cout << "	Val 1: " << mMap[i].val1 << endl;
		cout << "	Val 2: " << mMap[i].val2 << endl;
		cout << "	Val 3: " << mMap[i].val3 << endl;
		cout << "	Val 4: " << mMap[i].val4 << endl;
		cout << endl;
	}
}


// Test to make sure the following methods are working properly:
// 1. Random Measurement Data Generator
// 2. Random Query Disc Generator
// 3. Random Query Poly Generator
// 4. Random Query Strip Generator
//
// Also test to make sure can write the above to file and parse
// the resultant files properly.
void TestFileWritersReaders
(
 int numPoints,
 int numQueries,
 double minPhi,
 double maxPhi,
 double minTheta,
 double maxTheta,
 double minRad,
 double maxRad
)
{
	int i,j;
	ofstream fp;
	ifstream ip;
	std::string filename;

	// Generate Random Measurement Data
	std::vector<pointing> points;
	std::vector<Measurement> ms;
	int64 order;
	GenRandomPhiTheta(numPoints,points,order,false,minPhi,maxPhi,minTheta,maxTheta);
	GenRandomData(numPoints,ms,minPhi,maxPhi,minTheta,maxTheta);

	// Generate Random Query Discs
	double phi,theta,radius;
	std::vector<DiscType> discs;
	DiscType d;
	for(i = 0; i < numQueries; i++) {
		RandomDisc(phi,theta,radius,minPhi,maxPhi,minTheta,maxTheta,minRad,maxRad);
		d.pt.phi = phi; d.pt.theta = theta; d.radius = radius;
		discs.push_back(d);
	}

	// Generate Random Query Polys
	std::vector<pointing> pts;
	std::vector<PolyType> polys;
	PolyType p;
	for(i = 0; i < numQueries; i++) {
		pts = RandomConvexPoly(false,minPhi,maxPhi,minTheta,maxTheta,minRad,maxRad);
		p.pts.clear(); 
		for(j = 0; j < pts.size(); j++) {
			p.pts.push_back(pts[j]);
		}
		polys.push_back(p);
	}

	// Generate Random Query Strips
	double LO_THETA = 0.01*D2R; 
	double HI_THETA = 179.99*D2R; 
	std::vector<StripType> strips;
	StripType s;
	for( i = 0; i < numQueries; i++ ) {
		s.theta1 = LO_THETA + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(HI_THETA-LO_THETA)));
		s.theta2 = LO_THETA + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(HI_THETA-LO_THETA)));
		strips.push_back(s);
	}

	// Write Measurement Data To File
	filename = "rand_measurements.txt";
	fp.open(filename.c_str());
	if(fp.fail()){
	  cout << "Unable to open " << filename.c_str() << "!" << endl;
	  exit(1);     
	}
	WriteRandomMeasurementsHeader(fp,numPoints,order);
	WriteRandomMeasurements(fp,ms);
	fp.close();

	// Write Discs To File
	filename = "rand_discs.txt";
	fp.open(filename.c_str());
	if(fp.fail()){
	  cout << "Unable to open " << filename.c_str() << "!" << endl;
	  exit(1);     
	}
	WriteRandomDiscHeader(fp,numQueries);
	WriteRandomDiscs(fp,discs);
	fp.close();

	// Write Polys To File
	filename = "rand_polys.txt";
	fp.open(filename.c_str());
	if(fp.fail()){
	  cout << "Unable to open " << filename.c_str() << "!" << endl;
	  exit(1);     
	}
	WriteRandomPolysHeader(fp,numQueries);
	WriteRandomPolys(fp,polys);
	fp.close();

	// Write Strips To File
	filename = "rand_strips.txt";
	fp.open(filename.c_str());
	if(fp.fail()){
	  cout << "Unable to open " << filename.c_str() << "!" << endl;
	  exit(1);     
	}
	WriteRandomStripsHeader(fp,numQueries);
	WriteRandomStrips(fp,strips);
	fp.close();

	// Read & STDOUT Measurements From File
	filename = "rand_measurements.txt";
	ip.open(filename.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << filename.c_str() << "!" << endl;
	  exit(1);     
	}
	ms = ReadRandomMeasurements(ip,order);
	ip.close();
	cout << "#### RANDOM MEASUREMENTS ####\n";
	cout << "HPX ORDER: " << order << endl;
	cout << "REC,PHI,THETA,D1,D2,D3,D4,V1,V2,V3,V4\n";
	for(i = 0; i < ms.size(); i++) {
		cout << ms[i].rec << " " <<  ms[i].pt.phi << " " << ms[i].pt.theta << " "
			<<  ms[i].data1 << " " <<  ms[i].data2 << " " <<  ms[i].data3 << " " <<  ms[i].data4 << " "
			<<  ms[i].val1 << " " <<  ms[i].val2 << " " <<  ms[i].val3 << " " <<  ms[i].val4 << endl;
	}

	// Read Discs From File
	filename = "rand_discs.txt";
	ip.open(filename.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << filename.c_str() << "!" << endl;
	  exit(1);     
	}
	discs = ReadRandomDiscs(ip);
	ip.close();
	cout << "\n#### RANDOM DISCS ####\n";
	cout << "REC,PHI,THETA,RADIUS\n";
	for(i = 0; i < discs.size(); i++) {
		cout << i << " " <<  discs[i].pt.phi << " " << discs[i].pt.theta << " " << discs[i].radius << endl;
	}

	// Read Polys From File
	filename = "rand_polys.txt";
	ip.open(filename.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << filename.c_str() << "!" << endl;
	  exit(1);     
	}
	polys = ReadRandomPolys(ip);
	ip.close();
	cout << "\n#### RANDOM POLYS ####\n";
	cout << "REC,NUMPTS,PHI,THETA,...\n";
	for(i = 0; i < polys.size(); i++) {
		cout << i << " " << polys[i].pts.size() << " ";
		for(j = 0; j < polys[i].pts.size(); j++) {
			cout << polys[i].pts[j].phi << " " << polys[i].pts[j].theta << " ";
		}
		cout << endl;
	}

	// Read Strips From File
	filename = "rand_strips.txt";
	ip.open(filename.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << filename.c_str() << "!" << endl;
	  exit(1);     
	}
	strips = ReadRandomStrips(ip);
	ip.close();
	cout << "\n#### RANDOM STRIPS ####\n";
	cout << "REC,THETA1,THETA2\n";
	for(i = 0; i < strips.size(); i++) {
		cout << i << " " <<  strips[i].theta1 << " " << strips[i].theta2 << endl;
	}
}

void CreateRandomMeasurementsFile
(
 int numPoints,
 std::string filename,
 double minPhi,
 double maxPhi,
 double minTheta,
 double maxTheta 
)
{
	ofstream fp;

	// Generate Random Measurement Data
	std::vector<pointing> points;
	std::vector<Measurement> ms;
	int64 order;
	GenRandomPhiTheta(numPoints,points,order,false,minPhi,maxPhi,minTheta,maxTheta);
	GenRandomData(numPoints,ms,minPhi,maxPhi,minTheta,maxTheta);

	// Write Measurement Data To File
	fp.open(filename.c_str());
	if(fp.fail()){
	  cout << "Unable to open " << filename.c_str() << "!" << endl;
	  exit(1);     
	}
	WriteRandomMeasurementsHeader(fp,numPoints,order);
	WriteRandomMeasurements(fp,ms);
	fp.close();
}

void CreateRandomMeasurementsFileOfDensity(int order,float density,std::string filename)
{
	// Given HPX order compute the NPIX (total number of cells in HPX at order). 
	// Next using given density, draw density*NPIX (where density is 0.0-1.0)
	// random HPX indices between 0-NPIX and random Measurement data and write 
	// to file.
	int MIN_INT = -999; int MAX_INT = 999;
	double MIN_DOUBLE = -9.0; double MAX_DOUBLE = 9.0;
	ofstream fp;
	Healpix_Custom mHPX;
	pointing pt;
	int64 HPXid = 0;
	int64 MIN_HPX = 0;
	int64 MAX_HPX = 0;
	int64 nSide  = int64(1)<<order;
	int64 npFace = nSide<<order;
	int64 nPix   = 12*npFace;

	// Compute number of cells required for density specification
	int64 nCells = int64(float(nPix)*density);
	MAX_HPX = nPix-1;

	std::vector<Measurement> ms(nCells);

	mHPX.Set(order,NEST);

	cout << "Constructing Random Measurement Data File by HPX File Density Specification" << endl;
	cout << "	Required order: " << order << endl;
	cout << "	npFace: " << npFace << endl;
	cout << "	nPix: " << nPix << endl;
	cout << "	Required density: " << density << endl;
	cout << "	Num Cells: " << nCells << endl;
	
	// Open Measurement Data File
	fp.open(filename.c_str());
	if(fp.fail()){
	  cout << "Unable to open " << filename.c_str() << "!" << endl;
	  exit(1);     
	}
	WriteRandomMeasurementsHeader(fp,nCells,order);

	// Draw "nCells" random HPXids in range MIN_HPX-MAX_HPX
	for(int i = 0; i < nCells; i++) {
		// Compute random HPX index
		HPXid = MIN_HPX + static_cast <int64> (rand()) /( static_cast <double> (RAND_MAX/double(MAX_HPX-MIN_HPX)));
		
		// Convert HPX index into pointing angle (phi,theta)
		pt = mHPX.pix2ang(HPXid);

		// Update Measurement data location with pointing angle
		ms[i].rec = i;
		ms[i].pt = pt;
		ms[i].data1 = MIN_INT + static_cast <int> (rand()) /( static_cast <double> (RAND_MAX/(MAX_INT-MIN_INT)));
		ms[i].data2 = MIN_INT + static_cast <int> (rand()) /( static_cast <double> (RAND_MAX/(MAX_INT-MIN_INT)));
		ms[i].data3 = MIN_INT + static_cast <int> (rand()) /( static_cast <double> (RAND_MAX/(MAX_INT-MIN_INT)));
		ms[i].data4 = MIN_INT + static_cast <int> (rand()) /( static_cast <double> (RAND_MAX/(MAX_INT-MIN_INT)));

		ms[i].val1 = MIN_DOUBLE + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(MAX_DOUBLE-MIN_DOUBLE)));
		ms[i].val2 = MIN_DOUBLE + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(MAX_DOUBLE-MIN_DOUBLE)));
		ms[i].val3 = MIN_DOUBLE + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(MAX_DOUBLE-MIN_DOUBLE)));
		ms[i].val4 = MIN_DOUBLE + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(MAX_DOUBLE-MIN_DOUBLE)));

	}
	WriteRandomMeasurements(fp,ms);
	fp.close();

}

void CreateRandomDiscsFile
(
 int numQueries,
 std::string filename,
 double minPhi,
 double maxPhi,
 double minTheta,
 double maxTheta,
 double minRad,
 double maxRad
)
{
	// Generate Random Query Discs
	int i;
	ofstream fp;
	double phi,theta,radius;
	std::vector<DiscType> discs;
	DiscType d;
	for(i = 0; i < numQueries; i++) {
		RandomDisc(phi,theta,radius,minPhi,maxPhi,minTheta,maxTheta,minRad,maxRad);
		d.pt.phi = phi; d.pt.theta = theta; d.radius = radius;
		discs.push_back(d);
	}

	// Write Discs To File
	fp.open(filename.c_str());
	if(fp.fail()){
	  cout << "Unable to open " << filename.c_str() << "!" << endl;
	  exit(1);     
	}
	WriteRandomDiscHeader(fp,numQueries);
	WriteRandomDiscs(fp,discs);
	fp.close();
}

void CreateRandomPolysFile
(
 int numQueries,
 std::string filename,
 double minPhi,
 double maxPhi,
 double minTheta,
 double maxTheta,
 double minRad,
 double maxRad
 )
{
	// Generate Random Query Polys
	int i,j;
	ofstream fp;
	std::vector<pointing> pts;
	std::vector<PolyType> polys;
	PolyType p;
	for(i = 0; i < numQueries; i++) {
		pts = RandomConvexPoly(false,minPhi,maxPhi,minTheta,maxTheta,minRad,maxRad);
		p.pts.clear();
		for(j = 0; j < pts.size(); j++) {
			p.pts.push_back(pts[j]);
		}
		polys.push_back(p);
	}

	// Write Polys To File
	fp.open(filename.c_str());
	if(fp.fail()){
	  cout << "Unable to open " << filename.c_str() << "!" << endl;
	  exit(1);     
	}
	WriteRandomPolysHeader(fp,numQueries);
	WriteRandomPolys(fp,polys);
	fp.close();
}

void CreateRandomStripsFile
(
 int numQueries,
 std::string filename,
 double minTheta,
 double maxTheta
 )
{
	// Generate Random Query Strips
	int i;
	ofstream fp;
	std::vector<StripType> strips;
	StripType s;
	for( i = 0; i < numQueries; i++ ) {
		s.theta1 = minTheta + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(maxTheta-minTheta)));
		s.theta2 = minTheta + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(maxTheta-minTheta)));
		strips.push_back(s);
	}

	// Write Strips To File
	fp.open(filename.c_str());
	if(fp.fail()){
	  cout << "Unable to open " << filename.c_str() << "!" << endl;
	  exit(1);     
	}
	WriteRandomStripsHeader(fp,numQueries);
	WriteRandomStrips(fp,strips);
	fp.close();
}

void CreateRandomNeighborsFile
(
 int numQueries,
 std::string dataFile,
 std::string outFile
)
{
	int i,index;
	int64 order;
	pointing nextPt;
	ifstream ip;
	ofstream fp;
	std::vector<NeighborType> neighbors;
	std::vector<Measurement> ms;
	NeighborType n;
	std::vector<pointing> points;
	int MIN_INT,MAX_INT;

	// Parse the point data file 
	ip.open(dataFile.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << dataFile.c_str() << "!" << endl;
	  exit(1);     
	}
	ms = ReadRandomMeasurements(ip,order);
	ip.close();

	// Now randomly select numQueries from point data list to use
	// for neighbor queries
	MIN_INT = 0;
	MAX_INT = ms.size()-1;
	for( i = 0; i < numQueries; i++ ) 
	{
		// Random index draw
		index = MIN_INT + rand() % (MAX_INT-MIN_INT+1);
		n.pt = ms[index].pt;
		neighbors.push_back(n);
	}

	// Write Neighbor Indices To File
	fp.open(outFile.c_str());
	if(fp.fail()){
	  cout << "Unable to open " << outFile.c_str() << "!" << endl;
	  exit(1);     
	}
	WriteRandomNeighborsHeader(fp,numQueries);
	WriteRandomNeighbors(fp,neighbors);
	fp.close();
}


void MRHPointInsertionAndSave(std::string dataInputFile, std::string outFile)
{
	int i;
	int64 order;
	ifstream ip;
	ofstream ofs("test_mrh_map.txt");
	std::vector<Measurement> ms;
	// Open and Parse Data Input File
	ip.open(dataInputFile.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << dataInputFile.c_str() << "!" << endl;
	  exit(1);     
	}
	ms = ReadRandomMeasurements(ip,order);
	ip.close();
	MultiResHpx_Map<Measurement> mMRH(order,NEST);
	cout << "First load data points into MRH data structure...\n";
	for( i = 0; i < ms.size(); i++)
	{
		if( i%1000 == 0 ) {
			cout << "	Inserting Point " << i << " of " << ms.size() << endl;
		}
		mMRH.AddRecord(ms[i],ms[i].pt);
	}
	cout << endl << endl;

	// Now save the MRH data structure to file
	cout << "Now save MRH data structure to file...\n\n";
	mMRH.SaveToFile(outFile);
	cout << "MRH data structure saved to: " << outFile.c_str() << endl;

	// Write out the MRH Map to file
	//ofs << mMRH;
	//ofs.close();

	//// TEST: Read MRH Map back in.
	//ifstream ifs("test_mrh_map.txt");
	//if(ifs >> mMRH)
	//{
	//	cout << "Read MRH Map back into memory from file!\n";
	//}



}

void MRHPointInsertionBenchMark(std::string dataInputFile, std::string outFile, int numTrials)
{
	ofstream fp;
	ifstream ip;
	uint64 start_s,stop_s;
	int i,j;
	int64 order,hpxid;
	std::vector<Measurement> ms;
	std::vector<int64> foundV;
	double avgPtInsertTime = 0.0;
	double avgPtInsertTrials = 0.0;

	// Open and Parse Data Input File
	ip.open(dataInputFile.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << dataInputFile.c_str() << "!" << endl;
	  exit(1);     
	}
	ms = ReadRandomMeasurements(ip,order);
	ip.close();


	// Open Results Output File
	fp.open(outFile.c_str());
	if(fp.fail()){
	  cout << "Unable to open " << outFile.c_str() << "!" << endl;
	  exit(1);     
	}

	cout << "#############################################" << endl;
	cout << "#### MRH POINT INSERTION BENCHMARK ####" << endl;
	cout << "#############################################" << endl << endl;

	fp << "#############################################" << "\n";
	fp << "#### MRH POINT INSERTION BENCHMARK ####" << "\n";
	fp << "#############################################" << "\n" << "\n";

	// Build MultiResHpx Data Structure
	// START Benchmark
	for( j = 0; j < numTrials; j++ )
	{
		cout << "Begin " << ms.size() << " MRH Point Insertions...\n";
		fp << "Begin " << ms.size() << " MRH Point Insertions...\n";
		MultiResHpx_Map<Measurement> mMRH(order,NEST);
		start_s=GetTimeMs64();
		for( i = 0; i < ms.size(); i++)
		{
			mMRH.AddRecord(ms[i],ms[i].pt);
		}
		// STOP Benchmark
		stop_s=GetTimeMs64();
		avgPtInsertTime = double(stop_s-start_s)/double(ms.size());
		avgPtInsertTrials += avgPtInsertTime;

		// Report Trial Results
		cout << "	Completed " << ms.size() << " MRH Point Insertions in " << double(stop_s-start_s) 
			<< " ms! Average MRH Point Insertion: " << avgPtInsertTime << " ms\n\n";

		// Write out Trial Result to output file
		fp << "	Completed " << ms.size() << " MRH Point Insertions in " << double(stop_s-start_s) 
			<< " ms! Average MRH Point Insertion: " << avgPtInsertTime << " ms\n\n";
	}
	fp << "\n\n";
	cout << "\n\n";
	
	// Report Results of All Trials
	cout << "\nCompleted ALL Trials of " << numTrials*ms.size() << " Point Insertion in " << avgPtInsertTrials 
		 << " ms! Average MRH Point Insertion: " << avgPtInsertTrials/double(numTrials) << " ms\n";

	// Write out results of All Trials to output file
	fp << "\nCompleted ALL Trials of " << numTrials*ms.size() << " Point Insertion in " << avgPtInsertTrials 
		 << " ms! Average MRH Point Insertion: " << avgPtInsertTrials/double(numTrials) << " ms\n";

	fp.close();

}


void HPXPointInsertionBenchMark(std::string dataInputFile, std::string outFile, int numTrials)
{
	ofstream fp;
	ifstream ip;
	uint64 start_s,stop_s;
	int i,j;
	int64 order,hpxid;
	std::vector<Measurement> ms;
	std::vector<int64> foundV;
	double avgPtInsertTime = 0.0;
	double avgPtInsertTrials = 0.0;

	// Open and Parse Data Input File
	ip.open(dataInputFile.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << dataInputFile.c_str() << "!" << endl;
	  exit(1);     
	}
	ms = ReadRandomMeasurements(ip,order);
	ip.close();

	// Open Results Output File
	fp.open(outFile.c_str());
	if(fp.fail()){
	  cout << "Unable to open " << outFile.c_str() << "!" << endl;
	  exit(1);     
	}

	cout << "#############################################" << endl;
	cout << "#### HPX POINT INSERTION BENCHMARK ####" << endl;
	cout << "#############################################" << endl << endl;

	fp << "#############################################" << "\n";
	fp << "#### HPX POINT INSERTION BENCHMARK ####" << "\n";
	fp << "#############################################" << "\n" << "\n";


	// Build HPX Data Structure
	// START Benchmark
	for( j = 0; j < numTrials; j++ )
	{
		cout << "Begin " << ms.size() << " HPX Point Insertions...\n";
		Healpix_Map<Measurement> mHPX(order,NEST);
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
		avgPtInsertTrials += avgPtInsertTime;

		// Report Trial Results
		cout << "	Completed " << ms.size() << " HPX Point Insertions in " << std::setprecision(20) << double(stop_s-start_s) 
			<< " ms! Average HPX Point Insertion: " << avgPtInsertTime << " ms\n\n";
		
		fp << "	Completed " << ms.size() << " HPX Point Insertions in " << std::setprecision(20) << double(stop_s-start_s) 
			<< " ms! Average HPX Point Insertion: " << avgPtInsertTime << " ms\n\n";
	}
	cout << "\n\n";
	fp << "\n\n";

	// Report Results of All Trials
	cout << "\nCompleted ALL Trials of " << std::setprecision(20) << numTrials*ms.size() << " Point Insertion in " << avgPtInsertTrials 
		 << " ms! Average HPX Point Insertion: " << avgPtInsertTrials/double(numTrials) << " ms\n";
	fp << "\nCompleted ALL Trials of " << std::setprecision(20) << numTrials*ms.size() << " Point Insertion in " << avgPtInsertTrials 
		 << " ms! Average HPX Point Insertion: " << avgPtInsertTrials/double(numTrials) << " ms\n";

	fp.close();
}


// #### OFFICIAL MRH DISC QUERY BENCHMARK ####
// Reads in file that contains "Measurement" data that is inserted into the MRH data structure.
// In addition, reads in a query input file that contains all the queries that will be used for
// the bench mark test. Last parameter is the number of trials. Basically how many repeats of the
// query bench mark to perform. This helps smooth out the average query estimation by repeating the
// same batch of queries multiple times.
void MRHDiscQueryBenchMark(std::string dataInputFile, std::string queryInputFile, std::string outFile, int numTrials)
{
	ofstream fp;
	ifstream ip;
	int i,j;
	int64 order;
	uint64 start_s,stop_s;
	std::vector<Measurement> ms;
	std::vector<Measurement> found;
	std::vector<DiscType> discs;
	pointing pt;
	double avgQueryTime = 0.0;
	double avgQueryTrials = 0.0;

	// Open and Parse Data Input File
	ip.open(dataInputFile.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << dataInputFile.c_str() << "!" << endl;
	  exit(1);     
	}
	ms = ReadRandomMeasurements(ip,order);
	ip.close();

	// Open Results Output File
	fp.open(outFile.c_str());
	if(fp.fail()){
	  cout << "Unable to open " << outFile.c_str() << "!" << endl;
	  exit(1);     
	}

	cout << "##################################" << endl;
	cout << "#### MRH DISC QUERY BENCHMARK ####" << endl;
	cout << "##################################" << endl << endl;
	
	fp << "##################################" << "\n";
	fp << "#### MRH DISC QUERY BENCHMARK ####" << "\n";
	fp << "##################################" << "\n" << "\n";

	// Build MultiResHpx Data Structure
	cout << "Begin adding " << ms.size() << " records to the Data Structure...\n";
	fp << "Begin adding " << ms.size() << " records to the Data Structure...\n";
	MultiResHpx_Map<Measurement> mMRH(order,NEST);
	for( i = 0; i < ms.size(); i++)
	{
		if( i%1000 == 0 ) {
			cout << "	Inserting Point " << i+1 << " of " << ms.size() << endl;
			fp << "	Inserting Point " << i+1 << " of " << ms.size() << endl;
		}
		mMRH.AddRecord(ms[i],ms[i].pt);
	}
	cout << "Completed adding " << ms.size() << " records to the Data Structure!\n";
	fp << "Completed adding " << ms.size() << " records to the Data Structure!\n";

	// Open and Parse Query Input FIle
	ip.open(queryInputFile.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << queryInputFile.c_str() << "!" << endl;
	  exit(1);     
	}
	discs = ReadRandomDiscs(ip);
	ip.close();

	// START Benchmark
	for( j = 0; j < numTrials; j++ )
	{
		cout << "Begin " << discs.size() << " Disc Queries..." << endl;
		fp << "Begin " << discs.size() << " Disc Queries..." << endl;
		start_s=GetTimeMs64();
		for( i = 0; i < discs.size(); i++)
		{          
			found = mMRH.QueryDisc(discs[i].pt,discs[i].radius);
		}

		// STOP Benchmark
		stop_s=GetTimeMs64();
		avgQueryTime = double(stop_s-start_s)/double(discs.size());
		avgQueryTrials += double(stop_s-start_s);
		cout << "	End " << discs.size() << " Disc Queries!" << endl;
		fp << "	End " << discs.size() << " Disc Queries!" << endl;


		// Report Trial Results
		cout << "	Completed " << discs.size() << " Disc Queries in " << double(stop_s-start_s) 
			<< " ms! Average Query: " << avgQueryTime << " ms\n\n";
		fp << "	Completed " << discs.size() << " Disc Queries in " << double(stop_s-start_s) 
			<< " ms! Average Query: " << avgQueryTime << " ms\n\n";
	}
	cout << endl << endl;
	fp << "\n\n";

	// Report Results of All Trials
	cout << "\nCompleted ALL Trials of " << numTrials*discs.size() << " Disc Queries in " << avgQueryTrials 
		<< " ms! Average Query: " << avgQueryTrials/double(numTrials*discs.size()) << " ms\n";
	fp << "\nCompleted ALL Trials of " << numTrials*discs.size() << " Disc Queries in " << avgQueryTrials 
		<< " ms! Average Query: " << avgQueryTrials/double(numTrials*discs.size()) << " ms\n";
	fp.close();
}


// #### OFFICIAL HPX DISC QUERY BENCHMARK ####
// Reads in file that contains "Measurement" data that is inserted into the HPX data structure.
// In addition, reads in a query input file that contains all the queries that will be used for
// the bench mark test. Last parameter is the number of trials. Basically how many repeats of the
// query bench mark to perform. This helps smooth out the average query estimation by repeating the
// same batch of queries multiple times.
void HPXDiscQueryBenchMark(std::string dataInputFile, std::string queryInputFile, std::string outFile, int numTrials)
{
	ofstream fp;
	ifstream ip;
	int i,j,k;
	int64 order,hpxid;
	uint64 start_s,stop_s;
	std::vector<Measurement> ms;
	std::vector<Measurement> found;
	std::vector<DiscType> discs;
	pointing pt;
	rangeset<int64> pixset;
	std::vector<int64> foundV;
	double avgQueryTime = 0.0;
	double avgQueryTrials = 0.0;

	// Open and Parse Data Input File
	ip.open(dataInputFile.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << dataInputFile.c_str() << "!" << endl;
	  exit(1);     
	}
	ms = ReadRandomMeasurements(ip,order);
	ip.close();

	// Open Results Output File
	fp.open(outFile.c_str());
	if(fp.fail()){
	  cout << "Unable to open " << outFile.c_str() << "!" << endl;
	  exit(1);     
	}

	cout << "##################################" << endl;
	cout << "#### HPX DISC QUERY BENCHMARK ####" << endl;
	cout << "##################################" << endl << endl;

	fp << "##################################" << "\n";
	fp << "#### HPX DISC QUERY BENCHMARK ####" << "\n";
	fp << "##################################" << "\n" << "\n";

	// Build HPX Data Structure
	cout << "Begin adding " << ms.size() << " records to the Data Structure...\n";
	fp << "Begin adding " << ms.size() << " records to the Data Structure...\n";
	Healpix_Map<Measurement> mHPX(order,NEST);
	for( i = 0; i < ms.size(); i++)
	{
		if( i%1000 == 0 ) {
			cout << "	Inserting Point " << i+1 << " of " << ms.size() << endl;
			fp << "	Inserting Point " << i+1 << " of " << ms.size() << endl;
		}
		// Compute HEALPix index of next GIS longitude,latitude,data index tuple
		hpxid = mHPX.ang2pix(ms[i].pt);

		// Insert data index associated with
		mHPX[hpxid] = ms[i]; 
	}
	cout << "Completed adding " << ms.size() << " records to the Data Structure!\n";
	fp << "Completed adding " << ms.size() << " records to the Data Structure!\n";

	// Open and Parse Query Input FIle
	ip.open(queryInputFile.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << queryInputFile.c_str() << "!" << endl;
	  exit(1);     
	}
	discs = ReadRandomDiscs(ip);
	ip.close();

	// START Benchmark
	for( j = 0; j < numTrials; j++ )
	{
		cout << "Begin " << discs.size() << " Disc Queries..." << endl;
		fp << "Begin " << discs.size() << " Disc Queries..." << endl;
		start_s=GetTimeMs64();
		for( i = 0; i < discs.size(); i++)
		{          
			mHPX.query_disc(discs[i].pt,discs[i].radius,pixset);
			
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
		avgQueryTime = double(stop_s-start_s)/double(discs.size());
		avgQueryTrials += double(stop_s-start_s);
		cout << "	End " << discs.size() << " Disc Queries!" << endl;
		fp << "	End " << discs.size() << " Disc Queries!" << endl;


		// Report Trial Results
		cout << "	Completed " << discs.size() << " Disc Queries in " << double(stop_s-start_s) 
			<< " ms! Average Query: " << avgQueryTime << " ms\n\n";
		fp << "	Completed " << discs.size() << " Disc Queries in " << double(stop_s-start_s) 
			<< " ms! Average Query: " << avgQueryTime << " ms\n\n";
	}
	cout << endl << endl;
	fp << "\n\n";

	// Report Results of All Trials
	cout << "\nCompleted ALL Trials of " << numTrials*discs.size() << " Disc Queries in " << avgQueryTrials 
		<< " ms! Average Query: " << avgQueryTrials/double(numTrials*discs.size()) << " ms\n";
	fp << "\nCompleted ALL Trials of " << numTrials*discs.size() << " Disc Queries in " << avgQueryTrials 
		<< " ms! Average Query: " << avgQueryTrials/double(numTrials*discs.size()) << " ms\n";
	fp.close();
}


// Want to quantify any possible differences with output of respective MRH and HPX disc queries
// given same data sets and same disc queries.
void CompareDiscQueryOutput(std::string dataInputFile, std::string queryInputFile, std::string outFile)
{
	ofstream fp;
	ifstream ip;
	int i,j;
	int64 order,hpxid;
	int	totalMrhFound,totalHpxFound,totalMatches,totalMrhUnique,totalHpxUnique;
	int nMRHnodes,nHPXnodes;
	std::vector<Measurement> ms;
	std::vector<Measurement> foundMRH,foundHPX;
	std::vector<DiscType> discs;
	pointing pt;
	rangeset<int64> pixset;
	std::vector<int64> foundV;
	double avgQueryTime = 0.0;
	double avgQueryTrials = 0.0;
	int Matches = 0;
	int MRHUnique = 0;
	int HPXUnique = 0;

    // Open Query Results Output File
	fp.open(outFile.c_str());
	if(fp.fail()){
	  cout << "Unable to open " << outFile.c_str() << "!" << endl;
	  exit(1);     
	}

	// Open and Parse Data Input File
	ip.open(dataInputFile.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << dataInputFile.c_str() << "!" << endl;
	  exit(1);     
	}
	ms = ReadRandomMeasurements(ip,order);
	ip.close();

	cout << "###############################################" << endl;
	cout << "#### COMPARE MRH VS DISC DISC QUERY OUTPUT ####" << endl;
	cout << "###############################################" << endl << endl;

	fp << "###############################################" << "\n";
	fp << "#### COMPARE MRH VS DISC DISC QUERY OUTPUT ####" << "\n";
	fp << "###############################################" << "\n" << "\n";

	// Build MultiResHpx Data Structure
	cout << "Begin adding " << ms.size() << " records to the MRH Data Structure...\n";
	fp << "Begin adding " << ms.size() << " records to the MRH Data Structure...\n";
	MultiResHpx_Map<Measurement> mMRH(order,NEST);
	fp << "Data Set:\n";
	fp << "Phi,Theta,P,Z\n";
	for( i = 0; i < ms.size(); i++)
	{
		if( i%1000 == 0 ) {
			cout << "	Inserting Point " << i+1 << " of " << ms.size() << endl;
			fp << "	Inserting Point " << i+1 << " of " << ms.size() << "\n";
		}
		mMRH.AddRecord(ms[i],ms[i].pt);

		// Write out data location to query comparision file
		fp << ms[i].pt.phi << "," << ms[i].pt.theta << ","
		   << ms[i].pt.phi/pi << "," << cos(ms[i].pt.theta) << "\n";
	}
	fp << "\n\n";
	cout << "Completed adding " << ms.size() << " records to the MRH Data Structure!\n";
	fp << "Completed adding " << ms.size() << " records to the MRH Data Structure!\n";

	// Build HPX Data Structure
	cout << "Begin adding " << ms.size() << " records to the HPX Data Structure...\n";
	fp << "Begin adding " << ms.size() << " records to the HPX Data Structure...\n";
	Healpix_Map<Measurement> mHPX(order,NEST);
	for( i = 0; i < ms.size(); i++)
	{
		if( i%1000 == 0 ) {
			cout << "	Inserting Point " << i+1 << " of " << ms.size() << endl;
			fp << "	Inserting Point " << i+1 << " of " << ms.size() << "\n";
		}
		// Compute HEALPix index of next GIS longitude,latitude,data index tuple
		hpxid = mHPX.ang2pix(ms[i].pt);

		// Insert data index associated with
		mHPX[hpxid] = ms[i]; 
	}
	cout << "Completed adding " << ms.size() << " records to the HPX Data Structure!\n";
	fp << "Completed adding " << ms.size() << " records to the HPX Data Structure!\n";


	// Open and Parse Query Input FIle
	ip.open(queryInputFile.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << queryInputFile.c_str() << "!" << endl;
	  exit(1);     
	}
	discs = ReadRandomDiscs(ip);
	ip.close();

	cout << "Begin " << discs.size() << " Disc Queries..." << endl;
	fp << "Begin " << discs.size() << " Disc Queries..." << "\n";
	totalMrhFound = 0;totalHpxFound=0;totalMatches=0;totalMrhUnique=0;totalHpxUnique=0;
	for( i = 0; i < discs.size(); i++)
	{          
		Matches = 0; MRHUnique = 0; HPXUnique = 0;
		foundMRH.clear(); foundHPX.clear(); foundV.clear();

		// Do MRH Disc Query
		foundMRH = mMRH.QueryDisc(discs[i].pt,discs[i].radius);

		// Do HPX Disc Query
		mHPX.query_disc(discs[i].pt,discs[i].radius,pixset);
		
		// Convert cell rangeset to list of individual
		// cells.
		pixset.toVector(foundV);

		// Spin through HPX's found indices keep only the
		// non empty indices that pass the point-in-disc test
		for( j=0; j < foundV.size(); j++ ) {
			if( mHPX[foundV[j]].rec != EMPTY ) {
				if( IsPointInDisc(discs[i].pt,discs[i].radius,mHPX[foundV[j]].pt) ){
					foundHPX.push_back(mHPX[foundV[j]]);
				}
			} 
		}
		AnalyzeResults(foundHPX,foundMRH,Matches,HPXUnique,MRHUnique,false);
		totalMatches += Matches;
		totalMrhUnique += MRHUnique;
		totalHpxUnique += HPXUnique;
		totalMrhFound += foundMRH.size();
		totalHpxFound += foundHPX.size();
		cout << "Query " << i+1 << " Matches: " << Matches << " MRH Only: " << MRHUnique << " HPX Only: " << HPXUnique << endl;
		fp << "Query " << i+1 << " Matches: " << Matches << " MRH Only: " << MRHUnique << " HPX Only: " << HPXUnique << "\n";
	
		// Output Query and Query results to output file
		// Write Out Query # and the Query in phi,theta,p,z
		fp << "Query #" << i+1 << "\n";
		fp << "Phi,Theta,P,Z,Radius\n";
		fp << discs[i].pt.phi << "," << discs[i].pt.theta << ","
		   << discs[i].pt.phi/pi << "," << cos(discs[i].pt.theta) << ","
		   << discs[i].radius << "\n";

		// Now write out the MRH found results
		fp << "Query #" << i+1 << " MRH Found Results\n";
		fp << "Phi,Theta,P,Z\n";
		for( j=0; j < foundMRH.size(); j++) {
			fp << foundMRH[j].pt.phi << "," << foundMRH[j].pt.theta << ","
			   << foundMRH[j].pt.phi/pi << "," << cos(foundMRH[j].pt.theta) << "\n";
		}
		fp << "\n";

		// Now write out the HPX found results
		fp << "Query #" << i+1 << " HPX Found Results\n";
		fp << "Phi,Theta,P,Z\n";
		for( j=0; j < foundHPX.size(); j++) {
			fp << foundHPX[j].pt.phi << "," << foundHPX[j].pt.theta << ","
			   << foundHPX[j].pt.phi/pi << "," << cos(foundHPX[j].pt.theta) << "\n";	
		}
		fp << "\n\n";
	}
	cout << "	Finished Computing MRH vs HPX Disc Query Accuracy Test!\n";
	fp << "	Finished Computing MRH vs HPX Disc Query Accuracy Test!\n";

	// Compute/get statistics of the final data structures 
	nMRHnodes = mMRH.NumNodes();
	nHPXnodes = mHPX.Npix();
    // Get Min,Max,Average Depth of QuadTrees of each data stucture
	int maxDepthMRH = mMRH.MaxDepth();
	int minDepthMRH = mMRH.MinDepth();
	int avgDepthMRH = mMRH.AvgDepth();

	cout << "Number Points: " << ms.size() << endl;
	cout << "Number Disc Queries: " << discs.size() << endl;
	cout << "Total MRH Nodes: " << nMRHnodes << " Total HPX Nodes: " << nHPXnodes << endl;
	cout << "MRH Min Depth: " << minDepthMRH << " Max Depth: " << maxDepthMRH << " Avg Depth: " << avgDepthMRH << endl;
	cout << "HPX Depth: " << order << endl;
	cout << "Total Matches: " << totalMatches << endl;
	cout << "Total MRH Found: " << totalMrhFound << " Total MRH Unique: " << totalMrhUnique << endl;
	cout << "Total HPX Found: " << totalHpxFound << " Total HPX Unique: " << totalHpxUnique << endl;

	fp << "Number Points: " << ms.size() << "\n";
	fp << "Number Disc Queries: " << discs.size() << "\n";
	fp << "Total MRH Nodes: " << nMRHnodes << " Total HPX Nodes: " << nHPXnodes << "\n";
	fp << "MRH Min Depth: " << minDepthMRH << " Max Depth: " << maxDepthMRH << " Avg Depth: " << avgDepthMRH << "\n";
	fp << "HPX Depth: " << order << "\n";
	fp << "Total Matches: " << totalMatches << endl;
	fp << "Total MRH Found: " << totalMrhFound << " Total MRH Unique: " << totalMrhUnique << "\n";
	fp << "Total HPX Found: " << totalHpxFound << " Total HPX Unique: " << totalHpxUnique << "\n";

	fp.close();

}



// #### OFFICIAL MRH POLY QUERY BENCHMARK ####
// Reads in file that contains "Measurement" data that is inserted into the MRH data structure.
// In addition, reads in a query input file that contains all the queries that will be used for
// the bench mark test. Last parameter is the number of trials. Basically how many repeats of the
// query bench mark to perform. This helps smooth out the average query estimation by repeating the
// same batch of queries multiple times.
void MRHPolyQueryBenchMark(std::string dataInputFile, std::string queryInputFile, std::string outFile, int numTrials)
{
	ofstream fp;
	ifstream ip;
	int i,j;
	int64 order;
	uint64 start_s,stop_s;
	std::vector<Measurement> ms;
	std::vector<Measurement> found;
	std::vector<PolyType> polys;
	pointing pt;
	double avgQueryTime = 0.0;
	double avgQueryTrials = 0.0;

    // Open Query Results Output File
	fp.open(outFile.c_str());
	if(fp.fail()){
	  cout << "Unable to open " << outFile.c_str() << "!" << endl;
	  exit(1);     
	}

	// Open and Parse Data Input File
	ip.open(dataInputFile.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << dataInputFile.c_str() << "!" << endl;
	  exit(1);     
	}
	ms = ReadRandomMeasurements(ip,order);
	ip.close();

	cout << "##################################" << endl;
	cout << "#### MRH POLY QUERY BENCHMARK ####" << endl;
	cout << "##################################" << endl << endl;

	fp << "##################################" << "\n";
	fp << "#### MRH POLY QUERY BENCHMARK ####" << "\n";
	fp << "##################################" << "\n" << "\n";


	// Build MultiResHpx Data Structure
	cout << "Begin adding " << ms.size() << " records to the Data Structure...\n";
	fp << "Begin adding " << ms.size() << " records to the Data Structure...\n";
	MultiResHpx_Map<Measurement> mMRH(order,NEST);
	for( i = 0; i < ms.size(); i++)
	{
		if( i%1000 == 0 ) {
			cout << "	Inserting Point " << i+1 << " of " << ms.size() << endl;
			fp << "	Inserting Point " << i+1 << " of " << ms.size() << "\n";
		}
		mMRH.AddRecord(ms[i],ms[i].pt);
	}
	cout << "Completed adding " << ms.size() << " records to the Data Structure!\n";
	fp << "Completed adding " << ms.size() << " records to the Data Structure!\n";

	// Open and Parse Query Input FIle
	ip.open(queryInputFile.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << queryInputFile.c_str() << "!" << endl;
	  exit(1);     
	}
	polys = ReadRandomPolys(ip);
	ip.close();

	// START Benchmark
	for( j = 0; j < numTrials; j++ )
	{
		cout << "Begin " << polys.size() << " Poly Queries..." << endl;
		fp << "Begin " << polys.size() << " Poly Queries..." << "\n";
		start_s=GetTimeMs64();
		for( i = 0; i < polys.size(); i++)
		{          
			found = mMRH.QueryPolygon(polys[i].pts);
		}

		// STOP Benchmark
		stop_s=GetTimeMs64();
		avgQueryTime = double(stop_s-start_s)/double(polys.size());
		avgQueryTrials += double(stop_s-start_s);
		cout << "	End " << polys.size() << " Poly Queries!" << endl;
		fp << "	End " << polys.size() << " Poly Queries!" << "\n";

		// Report Trial Results
		cout << "	Completed " << polys.size() << " Poly Queries in " << double(stop_s-start_s) 
			<< " ms! Average Query: " << avgQueryTime << " ms\n\n";
		fp << "	Completed " << polys.size() << " Poly Queries in " << double(stop_s-start_s) 
			<< " ms! Average Query: " << avgQueryTime << " ms\n\n";
	}
	cout << endl << endl;
	fp << "\n\n";

	// Report Results of All Trials
	cout << "\nCompleted ALL Trials of " << numTrials*polys.size() << " Poly Queries in " << avgQueryTrials 
		<< " ms! Average Query: " << avgQueryTrials/double(numTrials*polys.size()) << " ms\n";
	fp << "\nCompleted ALL Trials of " << numTrials*polys.size() << " Poly Queries in " << avgQueryTrials 
		<< " ms! Average Query: " << avgQueryTrials/double(numTrials*polys.size()) << " ms\n";
	fp.close();
}


// #### OFFICIAL HPX POLY QUERY BENCHMARK ####
// Reads in file that contains "Measurement" data that is inserted into the HPX data structure.
// In addition, reads in a query input file that contains all the queries that will be used for
// the bench mark test. Last parameter is the number of trials. Basically how many repeats of the
// query bench mark to perform. This helps smooth out the average query estimation by repeating the
// same batch of queries multiple times.
void HPXPolyQueryBenchMark(std::string dataInputFile, std::string queryInputFile, std::string outFile, int numTrials)
{
	ofstream fp;
	ifstream ip;
	int i,j,k;
	int64 order,hpxid;
	uint64 start_s,stop_s;
	std::vector<Measurement> ms;
	std::vector<Measurement> found;
	std::vector<PolyType> polys;
	pointing pt;
	rangeset<int64> pixset;
	std::vector<int64> foundV;
	double avgQueryTime = 0.0;
	double avgQueryTrials = 0.0;

    // Open Query Results Output File
	fp.open(outFile.c_str());
	if(fp.fail()){
	  cout << "Unable to open " << outFile.c_str() << "!" << endl;
	  exit(1);     
	}

	// Open and Parse Data Input File
	ip.open(dataInputFile.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << dataInputFile.c_str() << "!" << endl;
	  exit(1);     
	}
	ms = ReadRandomMeasurements(ip,order);
	ip.close();

	cout << "##################################" << endl;
	cout << "#### HPX POLY QUERY BENCHMARK ####" << endl;
	cout << "##################################" << endl << endl;

	fp << "##################################" << "\n";
	fp << "#### HPX POLY QUERY BENCHMARK ####" << "\n";
	fp << "##################################" << "\n" << "\n";

	// Build HPX Data Structure
	cout << "Begin adding " << ms.size() << " records to the Data Structure...\n";
	fp << "Begin adding " << ms.size() << " records to the Data Structure...\n";
	Healpix_Map<Measurement> mHPX(order,NEST);
	for( i = 0; i < ms.size(); i++)
	{
		if( i%1000 == 0 ) {
			cout << "	Inserting Point " << i+1 << " of " << ms.size() << endl;
			fp << "	Inserting Point " << i+1 << " of " << ms.size() << "\n";
		}
		// Compute HEALPix index of next GIS longitude,latitude,data index tuple
		hpxid = mHPX.ang2pix(ms[i].pt);

		// Insert data index associated with
		mHPX[hpxid] = ms[i]; 
	}
	cout << "Completed adding " << ms.size() << " records to the Data Structure!\n";
	fp << "Completed adding " << ms.size() << " records to the Data Structure!\n";

	// Open and Parse Query Input FIle
	ip.open(queryInputFile.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << queryInputFile.c_str() << "!" << endl;
	  exit(1);     
	}
	polys = ReadRandomPolys(ip);
	ip.close();

	// START Benchmark
	for( j = 0; j < numTrials; j++ )
	{
		cout << "Begin " << polys.size() << " Poly Queries..." << endl;
		fp << "Begin " << polys.size() << " Poly Queries..." << "\n";
		start_s=GetTimeMs64();
		for( i = 0; i < polys.size(); i++)
		{          
			mHPX.query_polygon_inclusive(polys[i].pts,pixset,1);
			
			// Convert cell rangeset to list of individual
			// cells.
			pixset.toVector(foundV);

			// Spin through HPX's found indices keep only the
			// non empty indices that pass the point-in-disc test
			for( k=0; k < foundV.size(); k++ ) {
				if( mHPX[foundV[k]].rec != EMPTY ) {
					if( IsPointInPoly(mHPX[foundV[k]].pt,polys[i].pts) ){
						found.push_back(mHPX[foundV[k]]);
					}
				} 
			}
		}

		// STOP Benchmark
		stop_s=GetTimeMs64();
		avgQueryTime = double(stop_s-start_s)/double(polys.size());
		avgQueryTrials += double(stop_s-start_s);
		cout << "	End " << polys.size() << " Poly Queries!" << endl;
		fp << "	End " << polys.size() << " Poly Queries!" << "\n";


		// Report Trial Results
		cout << "	Completed " << polys.size() << " Poly Queries in " << double(stop_s-start_s) 
			<< " ms! Average Query: " << avgQueryTime << " ms\n\n";
		fp << "	Completed " << polys.size() << " Poly Queries in " << double(stop_s-start_s) 
			<< " ms! Average Query: " << avgQueryTime << " ms\n\n";
	}
	cout << endl << endl;
	fp << "\n\n";

	// Report Results of All Trials
	cout << "\nCompleted ALL Trials of " << numTrials*polys.size() << " Disc Queries in " << avgQueryTrials 
		<< " ms! Average Query: " << avgQueryTrials/double(numTrials*polys.size()) << " ms\n";
	fp << "\nCompleted ALL Trials of " << numTrials*polys.size() << " Disc Queries in " << avgQueryTrials 
		<< " ms! Average Query: " << avgQueryTrials/double(numTrials*polys.size()) << " ms\n";
	fp.close();
}


// Want to quantify any possible differences with output of respective MRH and HPX poly queries
// given same data sets and same poly queries.
void ComparePolyQueryOutput(std::string dataInputFile, std::string queryInputFile, std::string outFile)
{
	ofstream fp;
	ifstream ip;
	int i,j;
	int64 order,hpxid;
	int	totalMrhFound,totalHpxFound,totalMatches,totalMrhUnique,totalHpxUnique;
	int nMRHnodes,nHPXnodes;
	std::vector<Measurement> ms;
	std::vector<Measurement> foundMRH,foundHPX;
	std::vector<PolyType> polys;
	pointing pt;
	rangeset<int64> pixset;
	std::vector<int64> foundV;
	double avgQueryTime = 0.0;
	double avgQueryTrials = 0.0;
	int Matches = 0;
	int MRHUnique = 0;
	int HPXUnique = 0;

    // Open Query Results Output File
	fp.open(outFile.c_str());
	if(fp.fail()){
	  cout << "Unable to open " << outFile.c_str() << "!" << endl;
	  exit(1);     
	}

	// Open and Parse Data Input File
	ip.open(dataInputFile.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << dataInputFile.c_str() << "!" << endl;
	  exit(1);     
	}
	ms = ReadRandomMeasurements(ip,order);
	ip.close();

	cout << "##############################################" << endl;
	cout << "#### COMPARE MRH VS HPX POLY QUERY OUTPUT ####" << endl;
	cout << "##############################################" << endl << endl;

	fp << "##############################################" << "\n";
	fp << "#### COMPARE MRH VS HPX POLY QUERY OUTPUT ####" << "\n";
	fp << "##############################################" << "\n" << "\n";

	// Build MultiResHpx Data Structure
	cout << "Begin adding " << ms.size() << " records to the MRH Data Structure...\n";
	fp << "Begin adding " << ms.size() << " records to the MRH Data Structure...\n";
	MultiResHpx_Map<Measurement> mMRH(order,NEST);
	fp << "Data Set:\n";
	fp << "Phi,Theta,P,Z\n";
	for( i = 0; i < ms.size(); i++)
	{
		if( i%1000 == 0 ) {
			cout << "	Inserting Point " << i+1 << " of " << ms.size() << endl;
			fp << "	Inserting Point " << i+1 << " of " << ms.size() << "\n";
		}
		mMRH.AddRecord(ms[i],ms[i].pt);
 		// Write out data location to query comparision file
		fp << ms[i].pt.phi << "," << ms[i].pt.theta << ","
		   << ms[i].pt.phi/pi << "," << cos(ms[i].pt.theta) << "\n";
	}
	fp << "\n\n";
	cout << "Completed adding " << ms.size() << " records to the MRH Data Structure!\n";
	fp << "Completed adding " << ms.size() << " records to the MRH Data Structure!\n";

	// Build HPX Data Structure
	cout << "Begin adding " << ms.size() << " records to the HPX Data Structure...\n";
	fp << "Begin adding " << ms.size() << " records to the HPX Data Structure...\n";
	Healpix_Map<Measurement> mHPX(order,NEST);
	for( i = 0; i < ms.size(); i++)
	{
		if( i%1000 == 0 ) {
			cout << "	Inserting Point " << i+1 << " of " << ms.size() << endl;
			fp << "	Inserting Point " << i+1 << " of " << ms.size() << endl;
		}
		// Compute HEALPix index of next GIS longitude,latitude,data index tuple
		hpxid = mHPX.ang2pix(ms[i].pt);

		// Insert data index associated with
		mHPX[hpxid] = ms[i]; 
	}
	cout << "Completed adding " << ms.size() << " records to the HPX Data Structure!\n";
	fp << "Completed adding " << ms.size() << " records to the HPX Data Structure!\n";


	// Open and Parse Query Input FIle
	ip.open(queryInputFile.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << queryInputFile.c_str() << "!" << endl;
	  exit(1);     
	}
	polys = ReadRandomPolys(ip);
	ip.close();

	cout << "Begin " << polys.size() << " Poly Queries..." << endl;
	fp << "Begin " << polys.size() << " Poly Queries..." << "\n";
	totalMrhFound = 0;totalHpxFound=0;totalMatches=0;totalMrhUnique=0;totalHpxUnique=0;
	for( i = 0; i < polys.size(); i++)
	{          
		Matches = 0; MRHUnique = 0; HPXUnique = 0;
		foundMRH.clear(); foundHPX.clear(); foundV.clear();

		// Do MRH Disc Query
		foundMRH = mMRH.QueryPolygon(polys[i].pts);

		// Do HPX Disc Query
		mHPX.query_polygon_inclusive(polys[i].pts,pixset,1);
		
		// Convert cell rangeset to list of individual
		// cells.
		pixset.toVector(foundV);

		// Spin through HPX's found indices keep only the
		// non empty indices that pass the point-in-disc test
		for( j=0; j < foundV.size(); j++ ) {
			if( mHPX[foundV[j]].rec != EMPTY ) {
				if( IsPointInPoly(mHPX[foundV[j]].pt,polys[i].pts) ){
					foundHPX.push_back(mHPX[foundV[j]]);
				}
			} 
		}
		AnalyzeResults(foundHPX,foundMRH,Matches,HPXUnique,MRHUnique,false);
		totalMatches += Matches;
		totalMrhUnique += MRHUnique;
		totalHpxUnique += HPXUnique;
		totalMrhFound += foundMRH.size();
		totalHpxFound += foundHPX.size();
		cout << "Query " << i+1 << " Matches: " << Matches << " MRH Only: " << MRHUnique << " HPX Only: " << HPXUnique << endl;
		fp << "Query " << i+1 << " Matches: " << Matches << " MRH Only: " << MRHUnique << " HPX Only: " << HPXUnique << "\n";

		// Output Query and Query results to output file
		// Write Out Query # and the Query in phi,theta,p,z
		fp << "Query #" << i+1 << "\n";
		fp << "Phi,Theta,P,Z,...\n";
		for( j=0; j < polys[i].pts.size(); j++) {
		   fp << polys[i].pts[j].phi << "," << polys[i].pts[j].theta << ","
		      << polys[i].pts[j].phi/pi << "," << cos(polys[i].pts[j].theta) << "\n";
		}
		fp << polys[i].pts[0].phi << "," << polys[i].pts[0].theta << ","
		      << polys[i].pts[0].phi/pi << "," << cos(polys[i].pts[0].theta) << "\n";
		fp << "\n";

		// Now write out the MRH found results
		fp << "Query # " << i+1 << " MRH Found Results\n";
		fp << "Phi,Theta,P,Z\n";
		for( j=0; j < foundMRH.size(); j++) {
			fp << foundMRH[j].pt.phi << "," << foundMRH[j].pt.theta << ","
			   << foundMRH[j].pt.phi/pi << "," << cos(foundMRH[j].pt.theta) << "\n";
		}
		fp << "\n";
		// Now write out the HPX found results
		fp << "Query # " << i+1 << " HPX Found Results\n";
		fp << "Phi,Theta,P,Z\n";
		for( j=0; j < foundHPX.size(); j++) {
			fp << foundHPX[j].pt.phi << "," << foundHPX[j].pt.theta << ","
			   << foundHPX[j].pt.phi/pi << "," << cos(foundHPX[j].pt.theta) << "\n";	
		}
		fp << "\n\n";
	}
	cout << "	Finished Computing MRH vs HPX Poly Query Accuracy Test!\n";
	fp << "	Finished Computing MRH vs HPX Poly Query Accuracy Test!\n";

	// Compute/get statistics of the final data structures 
	nMRHnodes = mMRH.NumNodes();
	nHPXnodes = mHPX.Npix();
    // Get Min,Max,Average Depth of QuadTrees of each data stucture
	int maxDepthMRH = mMRH.MaxDepth();
	int minDepthMRH = mMRH.MinDepth();
	int avgDepthMRH = mMRH.AvgDepth();

	cout << "Number Points: " << ms.size() << endl;
	cout << "Number Poly Queries: " << polys.size() << endl;
	cout << "Total MRH Nodes: " << nMRHnodes << " Total HPX Nodes: " << nHPXnodes << endl;
	cout << "MRH Min Depth: " << minDepthMRH << " Max Depth: " << maxDepthMRH << " Avg Depth: " << avgDepthMRH << endl;
	cout << "HPX Depth: " << order << endl;
	cout << "Total Matches: " << totalMatches << endl;
	cout << "Total MRH Found: " << totalMrhFound << " Total MRH Unique: " << totalMrhUnique << endl;
	cout << "Total HPX Found: " << totalHpxFound << " Total HPX Unique: " << totalHpxUnique << endl;

	fp << "Number Points: " << ms.size() << "\n";
	fp << "Number Poly Queries: " << polys.size() << "\n";
	fp << "Total MRH Nodes: " << nMRHnodes << " Total HPX Nodes: " << nHPXnodes << "\n";
	fp << "MRH Min Depth: " << minDepthMRH << " Max Depth: " << maxDepthMRH << " Avg Depth: " << avgDepthMRH << "\n";
	fp << "HPX Depth: " << order << "\n";
	fp << "Total Matches: " << totalMatches << "\n";
	fp << "Total MRH Found: " << totalMrhFound << " Total MRH Unique: " << totalMrhUnique << "\n";
	fp << "Total HPX Found: " << totalHpxFound << " Total HPX Unique: " << totalHpxUnique << "\n";

	fp.close();
}


// #### OFFICIAL MRH STRIP QUERY BENCHMARK ####
// Reads in file that contains "Measurement" data that is inserted into the MRH data structure.
// In addition, reads in a query input file that contains all the queries that will be used for
// the bench mark test. Last parameter is the number of trials. Basically how many repeats of the
// query bench mark to perform. This helps smooth out the average query estimation by repeating the
// same batch of queries multiple times.
void MRHStripQueryBenchMark(std::string dataInputFile, std::string queryInputFile, std::string outFile, int numTrials)
{
	ofstream fp;
	ifstream ip;
	int i,j;
	int64 order;
	uint64 start_s,stop_s;
	std::vector<Measurement> ms;
	std::vector<Measurement> found;
	std::vector<StripType> strips;
	pointing pt;
	double avgQueryTime = 0.0;
	double avgQueryTrials = 0.0;

    // Open Query Results Output File
	fp.open(outFile.c_str());
	if(fp.fail()){
	  cout << "Unable to open " << outFile.c_str() << "!" << endl;
	  exit(1);     
	}

	// Open and Parse Data Input File
	ip.open(dataInputFile.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << dataInputFile.c_str() << "!" << endl;
	  exit(1);     
	}
	ms = ReadRandomMeasurements(ip,order);
	ip.close();

	cout << "###################################" << endl;
	cout << "#### MRH STRIP QUERY BENCHMARK ####" << endl;
	cout << "###################################" << endl << endl;

	fp << "###################################" << "\n";
	fp << "#### MRH STRIP QUERY BENCHMARK ####" << "\n";
	fp << "###################################" << "\n" << "\n";

	// Build MultiResHpx Data Structure
	cout << "Begin adding " << ms.size() << " records to the Data Structure...\n";
	fp << "Begin adding " << ms.size() << " records to the Data Structure...\n";
	MultiResHpx_Map<Measurement> mMRH(order,NEST);
	for( i = 0; i < ms.size(); i++)
	{
		if( i%1000 == 0 ) {
			cout << "	Inserting Point " << i+1 << " of " << ms.size() << endl;
			fp << "	Inserting Point " << i+1 << " of " << ms.size() << "\n";
		}
		mMRH.AddRecord(ms[i],ms[i].pt);
	}
	cout << "Completed adding " << ms.size() << " records to the Data Structure!\n";
	fp << "Completed adding " << ms.size() << " records to the Data Structure!\n";

	// Open and Parse Query Input FIle
	ip.open(queryInputFile.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << queryInputFile.c_str() << "!" << endl;
	  exit(1);     
	}
	strips = ReadRandomStrips(ip);
	ip.close();

	// START Benchmark
	for( j = 0; j < numTrials; j++ )
	{
		cout << "Begin " << strips.size() << " Strip Queries..." << endl;
		fp << "Begin " << strips.size() << " Strip Queries..." << "\n";
		start_s=GetTimeMs64();
		for( i = 0; i < strips.size(); i++)
		{          
			found = mMRH.QueryStrip(strips[i].theta1,strips[i].theta2);
		}

		// STOP Benchmark
		stop_s=GetTimeMs64();
		avgQueryTime = double(stop_s-start_s)/double(strips.size());
		avgQueryTrials += double(stop_s-start_s);
		cout << "	End " << strips.size() << " Strip Queries!" << endl;
		fp << "	End " << strips.size() << " Strip Queries!" << "\n";


		// Report Trial Results
		cout << "	Completed " << strips.size() << " Strip Queries in " << double(stop_s-start_s) 
			<< " ms! Average Query: " << avgQueryTime << " ms\n\n";
		fp << "	Completed " << strips.size() << " Strip Queries in " << double(stop_s-start_s) 
			<< " ms! Average Query: " << avgQueryTime << " ms\n\n";
	}
	cout << endl << endl;
	fp << "\n\n";

	// Report Results of All Trials
	cout << "\nCompleted ALL Trials of " << numTrials*strips.size() << " Strip Queries in " << avgQueryTrials 
		<< " ms! Average Query: " << avgQueryTrials/double(numTrials*strips.size()) << " ms\n";
	fp << "\nCompleted ALL Trials of " << numTrials*strips.size() << " Strip Queries in " << avgQueryTrials 
		<< " ms! Average Query: " << avgQueryTrials/double(numTrials*strips.size()) << " ms\n";
	fp.close();
}


// #### OFFICIAL HPX STRIP QUERY BENCHMARK ####
// Reads in file that contains "Measurement" data that is inserted into the HPX data structure.
// In addition, reads in a query input file that contains all the queries that will be used for
// the bench mark test. Last parameter is the number of trials. Basically how many repeats of the
// query bench mark to perform. This helps smooth out the average query estimation by repeating the
// same batch of queries multiple times.
void HPXStripQueryBenchMark(std::string dataInputFile, std::string queryInputFile, std::string outFile, int numTrials)
{
	ofstream fp;
	ifstream ip;
	int i,j,k;
	int64 order,hpxid;
	uint64 start_s,stop_s;
	std::vector<Measurement> ms;
	std::vector<Measurement> found;
	std::vector<StripType> strips;
	pointing pt;
	rangeset<int64> pixset;
	std::vector<int64> foundV;
	double avgQueryTime = 0.0;
	double avgQueryTrials = 0.0;

    // Open Query Results Output File
	fp.open(outFile.c_str());
	if(fp.fail()){
	  cout << "Unable to open " << outFile.c_str() << "!" << endl;
	  exit(1);     
	}

	// Open and Parse Data Input File
	ip.open(dataInputFile.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << dataInputFile.c_str() << "!" << endl;
	  exit(1);     
	}
	ms = ReadRandomMeasurements(ip,order);
	ip.close();

	cout << "###################################" << endl;
	cout << "#### HPX STRIP QUERY BENCHMARK ####" << endl;
	cout << "###################################" << endl << endl;

	fp << "###################################" << "\n";
	fp << "#### HPX STRIP QUERY BENCHMARK ####" << "\n";
	fp << "###################################" << "\n" << "\n";

	// Build HPX Data Structure
	cout << "Begin adding " << ms.size() << " records to the Data Structure...\n";
	fp << "Begin adding " << ms.size() << " records to the Data Structure...\n";
	Healpix_Map<Measurement> mHPX(order,RING);
	for( i = 0; i < ms.size(); i++)
	{
		if( i%1000 == 0 ) {
			cout << "	Inserting Point " << i+1 << " of " << ms.size() << endl;
			fp << "	Inserting Point " << i+1 << " of " << ms.size() << "\n";
		}
		// Compute HEALPix index of next GIS longitude,latitude,data index tuple
		hpxid = mHPX.ang2pix(ms[i].pt);

		// Insert data index associated with
		mHPX[hpxid] = ms[i]; 
	}
	cout << "Completed adding " << ms.size() << " records to the Data Structure!\n";
	fp << "Completed adding " << ms.size() << " records to the Data Structure!\n";

	// Open and Parse Query Input FIle
	ip.open(queryInputFile.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << queryInputFile.c_str() << "!" << endl;
	  exit(1);     
	}
	strips = ReadRandomStrips(ip);
	ip.close();

	// START Benchmark
	for( j = 0; j < numTrials; j++ )
	{
		cout << "Begin " << strips.size() << " Strip Queries..." << endl;
		fp << "Begin " << strips.size() << " Strip Queries..." << "\n";
		start_s=GetTimeMs64();
		for( i = 0; i < strips.size(); i++)
		{          
			mHPX.query_strip(strips[i].theta1,strips[i].theta2,true,pixset);
			
			// Convert cell rangeset to list of individual
			// cells.
			pixset.toVector(foundV);

			// Spin through HPX's found indices keep only the
			// non empty indices that pass the point-in-disc test
			for( k=0; k < foundV.size(); k++ ) {
				if( mHPX[foundV[k]].rec != EMPTY ) {
					if( IsPointInStrip(strips[i].theta1,strips[i].theta2,mHPX[foundV[k]].pt) ){
						found.push_back(mHPX[foundV[k]]);
					}
				} 
			}
		}

		// STOP Benchmark
		stop_s=GetTimeMs64();
		avgQueryTime = double(stop_s-start_s)/double(strips.size());
		avgQueryTrials += double(stop_s-start_s);
		cout << "	End " << strips.size() << " Strip Queries!" << endl;
		fp << "	End " << strips.size() << " Strip Queries!" << "\n";


		// Report Trial Results
		cout << "	Completed " << strips.size() << " Strip Queries in " << double(stop_s-start_s) 
			<< " ms! Average Query: " << avgQueryTime << " ms\n\n";
		fp << "	Completed " << strips.size() << " Strip Queries in " << double(stop_s-start_s) 
			<< " ms! Average Query: " << avgQueryTime << " ms\n\n";
	}
	cout << endl << endl;
	fp << "\n\n";

	// Report Results of All Trials
	cout << "\nCompleted ALL Trials of " << numTrials*strips.size() << " Strip Queries in " << avgQueryTrials 
		<< " ms! Average Query: " << avgQueryTrials/double(numTrials*strips.size()) << " ms\n";
	fp << "\nCompleted ALL Trials of " << numTrials*strips.size() << " Strip Queries in " << avgQueryTrials 
		<< " ms! Average Query: " << avgQueryTrials/double(numTrials*strips.size()) << " ms\n";

	fp.close();
}


// Want to quantify any possible differences with output of respective MRH and HPX strip queries
// given same data sets and same strip queries.
void CompareStripQueryOutput(std::string dataInputFile, std::string queryInputFile, std::string outFile)
{
	ofstream fp;
	ifstream ip;
	int i,j;
	int64 order,hpxid;
	int	totalMrhFound,totalHpxFound,totalMatches,totalMrhUnique,totalHpxUnique;
	int nMRHnodes,nHPXnodes;
	std::vector<Measurement> ms;
	std::vector<Measurement> foundMRH,foundHPX;
	std::vector<StripType> strips;
	pointing pt;
	rangeset<int64> pixset;
	std::vector<int64> foundV;
	double avgQueryTime = 0.0;
	double avgQueryTrials = 0.0;
	int Matches = 0;
	int MRHUnique = 0;
	int HPXUnique = 0;

    // Open Query Results Output File
	fp.open(outFile.c_str());
	if(fp.fail()){
	  cout << "Unable to open " << outFile.c_str() << "!" << endl;
	  exit(1);     
	}

	// Open and Parse Data Input File
	ip.open(dataInputFile.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << dataInputFile.c_str() << "!" << endl;
	  exit(1);     
	}
	ms = ReadRandomMeasurements(ip,order);
	ip.close();

	cout << "###############################################" << endl;
	cout << "#### COMPARE MRH VS HPX STRIP QUERY OUTPUT ####" << endl;
	cout << "###############################################" << endl << endl;

	fp << "###############################################" << "\n";
	fp << "#### COMPARE MRH VS HPX STRIP QUERY OUTPUT ####" << "\n";
	fp << "###############################################" << "\n" << "\n";


	// Build MultiResHpx Data Structure
	cout << "Begin adding " << ms.size() << " records to the MRH Data Structure...\n";
	fp << "Begin adding " << ms.size() << " records to the MRH Data Structure...\n";
	MultiResHpx_Map<Measurement> mMRH(order,NEST);
	fp << "Data Set:\n";
	fp << "Phi,Theta,P,Z\n";
	for( i = 0; i < ms.size(); i++)
	{
		if( i%1000 == 0 ) {
			cout << "	Inserting Point " << i+1 << " of " << ms.size() << endl;
		}
		mMRH.AddRecord(ms[i],ms[i].pt);
 		// Write out data location to query comparision file
		fp << ms[i].pt.phi << "," << ms[i].pt.theta << ","
		   << ms[i].pt.phi/pi << "," << cos(ms[i].pt.theta) << "\n";
	}
	fp << "\n\n";
	cout << "Completed adding " << ms.size() << " records to the MRH Data Structure!\n";
	fp << "Completed adding " << ms.size() << " records to the MRH Data Structure!\n";

	// Build HPX Data Structure
	cout << "Begin adding " << ms.size() << " records to the HPX Data Structure...\n";
	fp << "Begin adding " << ms.size() << " records to the HPX Data Structure...\n";
	Healpix_Map<Measurement> mHPX(order,RING);
	for( i = 0; i < ms.size(); i++)
	{
		if( i%1000 == 0 ) {
			cout << "	Inserting Point " << i+1 << " of " << ms.size() << endl;
			fp << "	Inserting Point " << i+1 << " of " << ms.size() << "\n";
		}
		// Compute HEALPix index of next GIS longitude,latitude,data index tuple
		hpxid = mHPX.ang2pix(ms[i].pt);

		// Insert data index associated with
		mHPX[hpxid] = ms[i]; 
	}
	cout << "Completed adding " << ms.size() << " records to the HPX Data Structure!\n";
	fp << "Completed adding " << ms.size() << " records to the HPX Data Structure!\n";


	// Open and Parse Query Input FIle
	ip.open(queryInputFile.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << queryInputFile.c_str() << "!" << endl;
	  exit(1);     
	}
	strips = ReadRandomStrips(ip);
	ip.close();

	cout << "Begin " << strips.size() << " Strip Queries..." << endl;
	fp << "Begin " << strips.size() << " Strip Queries..." << "\n";
	totalMrhFound = 0;totalHpxFound=0;totalMatches=0;totalMrhUnique=0;totalHpxUnique=0;
	for( i = 0; i < strips.size(); i++)
	{          
		Matches = 0; MRHUnique = 0; HPXUnique = 0;
		foundMRH.clear(); foundHPX.clear(); foundV.clear();

		// Do MRH Strip Query
		foundMRH = mMRH.QueryStrip(strips[i].theta1,strips[i].theta2);

		// Do HPX Strip Query
		mHPX.query_strip(strips[i].theta1,strips[i].theta2,true,pixset);
		
		// Convert cell rangeset to list of individual
		// cells.
		pixset.toVector(foundV);

		// Spin through HPX's found indices keep only the
		// non empty indices that pass the point-in-disc test
		for( j=0; j < foundV.size(); j++ ) {
			if( mHPX[foundV[j]].rec != EMPTY ) {
				if( IsPointInStrip(strips[i].theta1,strips[i].theta2,mHPX[foundV[j]].pt) ){
					foundHPX.push_back(mHPX[foundV[j]]);
				}
			} 
		}
		AnalyzeResults(foundHPX,foundMRH,Matches,HPXUnique,MRHUnique,false);
		totalMatches += Matches;
		totalMrhUnique += MRHUnique;
		totalHpxUnique += HPXUnique;
		totalMrhFound += foundMRH.size();
		totalHpxFound += foundHPX.size();
		cout << "Query " << i+1 << " Matches: " << Matches << " MRH Only: " << MRHUnique << " HPX Only: " << HPXUnique << endl;
		fp << "Query " << i+1 << " Matches: " << Matches << " MRH Only: " << MRHUnique << " HPX Only: " << HPXUnique << "\n";
	
		// Output Query and Query results to output file
		// Write Out Query # and the Query in phi,theta,p,z
		fp << "Query #" << i+1 << "\n";
		fp << "Theta1,Theta2,Z1,Z2\n";
		fp << strips[i].theta1 << "," << strips[i].theta2 << ","
		   << cos(strips[i].theta1) << "," << cos(strips[i].theta2) << "\n";
		fp << "\n";

		// Now write out the MRH found results
		fp << "Query # " << i+1 << " MRH Found Results\n";
		fp << "Phi,Theta,P,Z\n";
		for( j=0; j < foundMRH.size(); j++) {
			fp << foundMRH[j].pt.phi << "," << foundMRH[j].pt.theta << ","
			   << foundMRH[j].pt.phi/pi << "," << cos(foundMRH[j].pt.theta) << "\n";
		}
		fp << "\n";
		// Now write out the HPX found results
		fp << "Query # " << i+1 << " HPX Found Results\n";
		fp << "Phi,Theta,P,Z\n";
		for( j=0; j < foundHPX.size(); j++) {
			fp << foundHPX[j].pt.phi << "," << foundHPX[j].pt.theta << ","
			   << foundHPX[j].pt.phi/pi << "," << cos(foundHPX[j].pt.theta) << "\n";	
		}
		fp << "\n\n";	
	}
	cout << "	Finished Computing MRH vs HPX Strip Query Accuracy Test!\n";
	fp << "	Finished Computing MRH vs HPX Strip Query Accuracy Test!\n";

	// Compute/get statistics of the final data structures 
	nMRHnodes = mMRH.NumNodes();
	nHPXnodes = mHPX.Npix();
    // Get Min,Max,Average Depth of QuadTrees of each data stucture
	int maxDepthMRH = mMRH.MaxDepth();
	int minDepthMRH = mMRH.MinDepth();
	int avgDepthMRH = mMRH.AvgDepth();

	cout << "Number Points: " << ms.size() << endl;
	cout << "Number Strip Queries: " << strips.size() << endl;
	cout << "Total MRH Nodes: " << nMRHnodes << " Total HPX Nodes: " << nHPXnodes << endl;
	cout << "MRH Min Depth: " << minDepthMRH << " Max Depth: " << maxDepthMRH << " Avg Depth: " << avgDepthMRH << endl;
	cout << "HPX Depth: " << order << endl;
	cout << "Total Matches: " << totalMatches << endl;
	cout << "Total MRH Found: " << totalMrhFound << " Total MRH Unique: " << totalMrhUnique << endl;
	cout << "Total HPX Found: " << totalHpxFound << " Total HPX Unique: " << totalHpxUnique << endl;

	fp << "Number Points: " << ms.size() << "\n";
	fp << "Number Strip Queries: " << strips.size() << "\n";
	fp << "Total MRH Nodes: " << nMRHnodes << " Total HPX Nodes: " << nHPXnodes << "\n";
	fp << "MRH Min Depth: " << minDepthMRH << " Max Depth: " << maxDepthMRH << " Avg Depth: " << avgDepthMRH << "\n";
	fp << "HPX Depth: " << order << "\n";
	fp << "Total Matches: " << totalMatches << "\n";
	fp << "Total MRH Found: " << totalMrhFound << " Total MRH Unique: " << totalMrhUnique << "\n";
	fp << "Total HPX Found: " << totalHpxFound << " Total HPX Unique: " << totalHpxUnique << "\n";

	fp.close();

}


// #### OFFICIAL MRH NEIGHBOR QUERY BENCHMARK ####
// Reads in file that contains "Measurement" data that is inserted into the MRH data structure.
// In addition, reads in a query input file that contains all the queries that will be used for
// the bench mark test. Last parameter is the number of trials. Basically how many repeats of the
// query bench mark to perform. This helps smooth out the average query estimation by repeating the
// same batch of queries multiple times.
void MRHNeighborQueryBenchMark(std::string dataInputFile, std::string queryInputFile, std::string outFile, int numTrials)
{
	ofstream fp;
	ifstream ip;
	int i,j;
	int64 order;
	uint64 start_s,stop_s;
	std::vector<Measurement> ms;
	std::vector<Measurement> found;
	std::vector<NeighborType> neighbors;
	pointing pt;
	double avgQueryTime = 0.0;
	double avgQueryTrials = 0.0;

    // Open Query Results Output File
	fp.open(outFile.c_str());
	if(fp.fail()){
	  cout << "Unable to open " << outFile.c_str() << "!" << endl;
	  exit(1);     
	}

	// Open and Parse Data Input File
	ip.open(dataInputFile.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << dataInputFile.c_str() << "!" << endl;
	  exit(1);     
	}
	ms = ReadRandomMeasurements(ip,order);
	ip.close();

	cout << "#######################################" << endl;
	cout << "#### MRH NEIGHBORS QUERY BENCHMARK ####" << endl;
	cout << "#######################################" << endl << endl;

	fp << "#######################################" << "\n";
	fp << "#### MRH NEIGHBORS QUERY BENCHMARK ####" << "\n";
	fp << "#######################################" << "\n" << "\n";

	// Build MultiResHpx Data Structure
	cout << "Begin adding " << ms.size() << " records to the Data Structure...\n";
	fp << "Begin adding " << ms.size() << " records to the Data Structure...\n";
	MultiResHpx_Map<Measurement> mMRH(order,NEST);
	for( i = 0; i < ms.size(); i++)
	{
		if( i%1000 == 0 ) {
			cout << "	Inserting Point " << i+1 << " of " << ms.size() << endl;
			fp << "	Inserting Point " << i+1 << " of " << ms.size() << "\n";
		}
		mMRH.AddRecord(ms[i],ms[i].pt);
	}
	cout << "Completed adding " << ms.size() << " records to the Data Structure!\n";
	fp << "Completed adding " << ms.size() << " records to the Data Structure!\n";

	ip.open(queryInputFile.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << queryInputFile.c_str() << "!" << endl;
	  exit(1);     
	}
	neighbors = ReadRandomNeighbors(ip);
	ip.close();

	// START Benchmark
	for( j = 0; j < numTrials; j++ )
	{
		cout << "Begin " << neighbors.size() << " Neighbor Queries..." << endl;
		fp << "Begin " << neighbors.size() << " Neighbor Queries..." << "\n";
		start_s=GetTimeMs64();
		for( i = 0; i < neighbors.size(); i++)
		{          
			found = mMRH.Neighbors(neighbors[i].pt,order);
		}

		// STOP Benchmark
		stop_s=GetTimeMs64();
		avgQueryTime = double(stop_s-start_s)/double(neighbors.size());
		avgQueryTrials += double(stop_s-start_s);
		cout << "	End " << neighbors.size() << " Neighbor Queries!" << endl;
		fp << "	End " << neighbors.size() << " Neighbor Queries!" << "\n";


		// Report Trial Results
		cout << "	Completed " << neighbors.size() << " Neighbor Queries in " << double(stop_s-start_s) 
			<< " ms! Average Query: " << avgQueryTime << " ms\n\n";
		fp << "	Completed " << neighbors.size() << " Neighbor Queries in " << double(stop_s-start_s) 
			<< " ms! Average Query: " << avgQueryTime << " ms\n\n";
	}
	cout << endl << endl;
	fp << "\n\n";

	// Report Results of All Trials
	cout << "\nCompleted ALL Trials of " << numTrials*neighbors.size() << " Neighbor Queries in " << avgQueryTrials 
		<< " ms! Average Query: " << avgQueryTrials/double(numTrials*neighbors.size()) << " ms\n";
	fp << "\nCompleted ALL Trials of " << numTrials*neighbors.size() << " Neighbor Queries in " << avgQueryTrials 
		<< " ms! Average Query: " << avgQueryTrials/double(numTrials*neighbors.size()) << " ms\n";
	fp.close();
}

// #### OFFICIAL HPX NEIGHBOR QUERY BENCHMARK ####
// Reads in file that contains "Measurement" data that is inserted into the HPX data structure.
// In addition, reads in a query input file that contains all the queries that will be used for
// the bench mark test. Last parameter is the number of trials. Basically how many repeats of the
// query bench mark to perform. This helps smooth out the average query estimation by repeating the
// same batch of queries multiple times.
void HPXNeighborQueryBenchMark(std::string dataInputFile, std::string queryInputFile, std::string outFile, int numTrials)
{
	ofstream fp;
	ifstream ip;
	int i,j,k;
	int64 order,hpxid;
	uint64 start_s,stop_s;
	std::vector<Measurement> ms;
	std::vector<Measurement> found;
	std::vector<NeighborType> neighbors;
	pointing pt;
	fix_arr<int64,8> pixset;
	std::vector<int64> foundV;
	double avgQueryTime = 0.0;
	double avgQueryTrials = 0.0;

    // Open Query Results Output File
	fp.open(outFile.c_str());
	if(fp.fail()){
	  cout << "Unable to open " << outFile.c_str() << "!" << endl;
	  exit(1);     
	}

	// Open and Parse Data Input File
	ip.open(dataInputFile.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << dataInputFile.c_str() << "!" << endl;
	  exit(1);     
	}
	ms = ReadRandomMeasurements(ip,order);
	ip.close();

	cout << "######################################" << endl;
	cout << "#### HPX NEIGHBOR QUERY BENCHMARK ####" << endl;
	cout << "######################################" << endl << endl;

	fp << "######################################" << "\n";
	fp << "#### HPX NEIGHBOR QUERY BENCHMARK ####" << "\n";
	fp << "######################################" << "\n" << "\n";

	// Build HPX Data Structure
	cout << "Begin adding " << ms.size() << " records to the Data Structure with HPX order " << order << "...\n";
	fp << "Begin adding " << ms.size() << " records to the Data Structure with HPX order " << order << "...\n";
	Healpix_Map<Measurement> mHPX(order,RING);
	for( i = 0; i < ms.size(); i++)
	{
		if( i%1000 == 0 ) {
			cout << "	Inserting Point " << i+1 << " of " << ms.size() << endl;
			fp << "	Inserting Point " << i+1 << " of " << ms.size() << "\n";
		}
		// Compute HEALPix index of next GIS longitude,latitude,data index tuple
		hpxid = mHPX.ang2pix(ms[i].pt);

		// Insert data index associated with
		mHPX[hpxid] = ms[i]; 
	}
	cout << "Completed adding " << ms.size() << " records to the Data Structure!\n";
	fp << "Completed adding " << ms.size() << " records to the Data Structure!\n";

	// Open and Parse Query Input FIle
	ip.open(queryInputFile.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << queryInputFile.c_str() << "!" << endl;
	  exit(1);     
	}
	neighbors = ReadRandomNeighbors(ip);
	ip.close();

	// START Benchmark
	for( j = 0; j < numTrials; j++ )
	{
		cout << "Begin " << neighbors.size() << " Neighbor Queries..." << endl;
		fp << "Begin " << neighbors.size() << " Neighbor Queries..." << "\n";
		start_s=GetTimeMs64();
		for( i = 0; i < neighbors.size(); i++)
		{          

			hpxid = mHPX.ang2pix(neighbors[i].pt);
			mHPX.neighbors(hpxid,pixset);
			
			// Spin through HPX's found indices keep only the
			// non empty indices that pass the point-in-disc test
			for( k=0; k < 8; k++ ) {
				if( mHPX[pixset[k]].rec != EMPTY ) {
					found.push_back(mHPX[pixset[k]]);
				} 
			}
		}

		// STOP Benchmark
		stop_s=GetTimeMs64();
		avgQueryTime = double(stop_s-start_s)/double(neighbors.size());
		avgQueryTrials += double(stop_s-start_s);
		cout << "	End " << neighbors.size() << " Neighbor Queries!" << endl;
		fp << "	End " << neighbors.size() << " Neighbor Queries!" << "\n";


		// Report Trial Results
		cout << "	Completed " << neighbors.size() << " Neighbor Queries in " << double(stop_s-start_s) 
			<< " ms! Average Query: " << avgQueryTime << " ms\n\n";
		fp << "	Completed " << neighbors.size() << " Neighbor Queries in " << double(stop_s-start_s) 
			<< " ms! Average Query: " << avgQueryTime << " ms\n\n";
	}
	cout << endl << endl;
	fp << "\n\n";

	// Report Results of All Trials
	cout << "\nCompleted ALL Trials of " << numTrials*neighbors.size() << " Neighbor Queries in " << avgQueryTrials 
		<< " ms! Average Query: " << avgQueryTrials/double(numTrials*neighbors.size()) << " ms\n";
	fp << "\nCompleted ALL Trials of " << numTrials*neighbors.size() << " Neighbor Queries in " << avgQueryTrials 
		<< " ms! Average Query: " << avgQueryTrials/double(numTrials*neighbors.size()) << " ms\n";

	fp.close();
}


// Want to quantify any possible differences with output of respective MRH and HPX neighbor queries
// given same data sets and same neighbor queries.
void CompareNeighborQueryOutput(std::string dataInputFile, std::string queryInputFile, std::string outFile)
{
	ofstream fp;
	ifstream ip;
	int i,j,k;
	int64 order,hpxid;
	int	totalMrhFound,totalHpxFound,totalMatches,totalMrhUnique,totalHpxUnique;
	int nMRHnodes,nHPXnodes;
	std::vector<Measurement> ms;
	std::vector<Measurement> foundMRH,foundHPX;
	std::vector<NeighborType> neighbors;
	pointing pt;
	rangeset<int64> pixset;
	fix_arr<int64,8> pixset2;
	std::vector<int64> foundV;
	double avgQueryTime = 0.0;
	double avgQueryTrials = 0.0;
	int Matches = 0;
	int MRHUnique = 0;
	int HPXUnique = 0;

    // Open Query Results Output File
	fp.open(outFile.c_str());
	if(fp.fail()){
	  cout << "Unable to open " << outFile.c_str() << "!" << endl;
	  exit(1);     
	}

	// Open and Parse Data Input File
	ip.open(dataInputFile.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << dataInputFile.c_str() << "!" << endl;
	  exit(1);     
	}
	ms = ReadRandomMeasurements(ip,order);
	ip.close();

	cout << "##################################################" << endl;
	cout << "#### COMPARE MRH VS HPX NEIGHBOR QUERY OUTPUT ####" << endl;
	cout << "##################################################" << endl << endl;

	fp << "##################################################" << "\n";
	fp << "#### COMPARE MRH VS HPX NEIGHBOR QUERY OUTPUT ####" << "\n";
	fp << "##################################################" << "\n" << "\n";

	// Build MultiResHpx Data Structure
	cout << "Begin adding " << ms.size() << " records to the MRH Data Structure...\n";
	fp << "Begin adding " << ms.size() << " records to the MRH Data Structure...\n";
	MultiResHpx_Map<Measurement> mMRH(order,NEST);
	fp << "Data Set:\n";
	fp << "Phi,Theta,P,Z\n";
	for( i = 0; i < ms.size(); i++)
	{
		if( i%1000 == 0 ) {
			cout << "	Inserting Point " << i+1 << " of " << ms.size() << endl;
			fp << "	Inserting Point " << i+1 << " of " << ms.size() << endl;
		}
		mMRH.AddRecord(ms[i],ms[i].pt);
 		// Write out data location to query comparision file
		fp << ms[i].pt.phi << "," << ms[i].pt.theta << ","
		   << ms[i].pt.phi/pi << "," << cos(ms[i].pt.theta) << "\n";
	}
	fp << "\n\n";
	cout << "Completed adding " << ms.size() << " records to the MRH Data Structure!\n";
	fp << "Completed adding " << ms.size() << " records to the MRH Data Structure!\n";

	// Build HPX Data Structure
	cout << "Begin adding " << ms.size() << " records to the HPX Data Structure...\n";
	fp << "Begin adding " << ms.size() << " records to the HPX Data Structure...\n";
	Healpix_Map<Measurement> mHPX(order,RING);
	for( i = 0; i < ms.size(); i++)
	{
		if( i%1000 == 0 ) {
			cout << "	Inserting Point " << i+1 << " of " << ms.size() << endl;
			fp << "	Inserting Point " << i+1 << " of " << ms.size() << "\n";
		}
		// Compute HEALPix index of next GIS longitude,latitude,data index tuple
		hpxid = mHPX.ang2pix(ms[i].pt);

		// Insert data index associated with
		mHPX[hpxid] = ms[i]; 
	}
	cout << "Completed adding " << ms.size() << " records to the HPX Data Structure!\n";
	fp << "Completed adding " << ms.size() << " records to the HPX Data Structure!\n";

	// Open and Parse Query Input FIle
	ip.open(queryInputFile.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << queryInputFile.c_str() << "!" << endl;
	  exit(1);     
	}
	neighbors = ReadRandomNeighbors(ip);
	ip.close();

	cout << "Begin " << neighbors.size() << " Neighbor Queries..." << endl;
	fp << "Begin " << neighbors.size() << " Neighbor Queries..." << "\n";
	totalMrhFound = 0;totalHpxFound=0;totalMatches=0;totalMrhUnique=0;totalHpxUnique=0;
	for( i = 0; i < neighbors.size(); i++)
	{          
		Matches = 0; MRHUnique = 0; HPXUnique = 0;
		foundMRH.clear(); foundHPX.clear(); foundV.clear();

		// Do MRH Neighbor Query
		foundMRH = mMRH.Neighbors(neighbors[i].pt,order);

		// Do HPX Strip Query
		hpxid = mHPX.ang2pix(neighbors[i].pt);
		mHPX.neighbors(hpxid,pixset2);
		
		// Spin through HPX's found indices keep only the
		// non empty indices that pass the point-in-disc test
		for( k=0; k < 8; k++ ) {
			if( mHPX[pixset2[k]].rec != EMPTY ) {
				foundHPX.push_back(mHPX[pixset2[k]]);
			} 
		}

		AnalyzeResults(foundHPX,foundMRH,Matches,HPXUnique,MRHUnique,false);
		totalMatches += Matches;
		totalMrhUnique += MRHUnique;
		totalHpxUnique += HPXUnique;
		totalMrhFound += foundMRH.size();
		totalHpxFound += foundHPX.size();
		cout << "Query " << i+1 << " Matches: " << Matches << " MRH Only: " << MRHUnique << " HPX Only: " << HPXUnique << endl;
		fp << "Query " << i+1 << " Matches: " << Matches << " MRH Only: " << MRHUnique << " HPX Only: " << HPXUnique << "\n";
		// Output Query and Query results to output file
		// Write Out Query # and the Query in phi,theta,p,z
		fp << "Query #" << i+1 << "\n";
		fp << "Phi,Theta,P,Z,...\n";
		   fp << neighbors[i].pt.phi << "," << neighbors[i].pt.theta << ","
		      << neighbors[i].pt.phi/pi << "," << cos(neighbors[i].pt.theta) << "\n";
		fp << "\n";

		// Now write out the MRH found results
		fp << "Query # " << i+1 << " MRH Found Results\n";
		fp << "Phi,Theta,P,Z\n";
		for( j=0; j < foundMRH.size(); j++) {
			fp << foundMRH[j].pt.phi << "," << foundMRH[j].pt.theta << ","
			   << foundMRH[j].pt.phi/pi << "," << cos(foundMRH[j].pt.theta) << "\n";
		}
		fp << "\n";
		// Now write out the HPX found results
		fp << "Query # " << i+1 << " HPX Found Results\n";
		fp << "Phi,Theta,P,Z\n";
		for( j=0; j < foundHPX.size(); j++) {
			fp << foundHPX[j].pt.phi << "," << foundHPX[j].pt.theta << ","
			   << foundHPX[j].pt.phi/pi << "," << cos(foundHPX[j].pt.theta) << "\n";	
		}
		fp << "\n\n";	
	}
	cout << "	Finished Computing MRH vs HPX Neighbor Query Accuracy Test!\n";
	fp << "	Finished Computing MRH vs HPX Neighbor Query Accuracy Test!\n";

	// Compute/get statistics of the final data structures 
	nMRHnodes = mMRH.NumNodes();
	nHPXnodes = mHPX.Npix();
    // Get Min,Max,Average Depth of QuadTrees of each data stucture
	int maxDepthMRH = mMRH.MaxDepth();
	int minDepthMRH = mMRH.MinDepth();
	int avgDepthMRH = mMRH.AvgDepth();

	cout << "Number Points: " << ms.size() << endl;
	cout << "Number Neighbor Queries: " << neighbors.size() << endl;
	cout << "Total MRH Nodes: " << nMRHnodes << " Total HPX Nodes: " << nHPXnodes << endl;
	cout << "MRH Min Depth: " << minDepthMRH << " Max Depth: " << maxDepthMRH << " Avg Depth: " << avgDepthMRH << endl;
	cout << "HPX Depth: " << order << endl;
	cout << "Total Matches: " << totalMatches << endl;
	cout << "Total MRH Found: " << totalMrhFound << " Total MRH Unique: " << totalMrhUnique << endl;
	cout << "Total HPX Found: " << totalHpxFound << " Total HPX Unique: " << totalHpxUnique << endl;

	fp << "Number Points: " << ms.size() << "\n";
	fp << "Number Neighbor Queries: " << neighbors.size() << "\n";
	fp << "Total MRH Nodes: " << nMRHnodes << " Total HPX Nodes: " << nHPXnodes << "\n";
	fp << "MRH Min Depth: " << minDepthMRH << " Max Depth: " << maxDepthMRH << " Avg Depth: " << avgDepthMRH << "\n";
	fp << "HPX Depth: " << order << "\n";
	fp << "Total Matches: " << totalMatches << "\n";
	fp << "Total MRH Found: " << totalMrhFound << " Total MRH Unique: " << totalMrhUnique << "\n";
	fp << "Total HPX Found: " << totalHpxFound << " Total HPX Unique: " << totalHpxUnique << "\n";

	fp.close();

}


void TwoPointCorrBinning(std::string dataInputFile,double radius,int verbose)
{
	int i,j;
	ifstream ip;
	int64 order;
	vector<pointing> points;
	std::vector<Measurement> ms;
	std::vector<std::pair<Measurement,Measurement>> found;
	
	// Open and Parse Data Input File
	ip.open(dataInputFile.c_str());
	if(ip.fail()){
	  cout << "Unable to open " << dataInputFile.c_str() << "!" << endl;
	  exit(1);     
	}
	ms = ReadRandomMeasurements(ip,order);
	ip.close();

	// Build MultiResHpx Data Structure
	cout << "Begin adding " << ms.size() << " records to the MRH Data Structure...\n";
	MultiResHpx_Map<Measurement> mMRH(order,NEST);
	for( i = 0; i < ms.size(); i++)
	{
		mMRH.AddRecord(ms[i],ms[i].pt);
	}
	cout << "Completed adding " << ms.size() << " records to the MRH Data Structure!\n";

	if(verbose) {
		for( i = 0; i < mMRH.NumRec(); i++ ) {
			cout << mMRH[i].rec << "," << mMRH[i].pt.phi << ","
				<< mMRH[i].pt.theta << "," << mMRH[i].pt.phi/pi << ","
				<< cos(mMRH[i].pt.theta) << endl;
		}
	}

	// Compute Two-Point Correlation Query
	cout << "Query radius," << radius << endl << endl;

	cout << "Now do Two-Point Corrolation Query!" << endl << endl;
	found = mMRH.TwoPointCorrBin(radius);

	// Output results
	if(verbose) {
		cout << "Two-Point Correlated Results:" << endl;
		for( j = 0; j < found.size(); j++ ) {
			cout << found[j].first.rec << ","
				 << found[j].first.pt.phi << ","
				 << found[j].first.pt.theta << ","
				 << found[j].first.pt.phi/pi << ","
				 << cos(found[j].first.pt.theta) << endl 
				 << found[j].second.rec << ","
				 << found[j].second.pt.phi << ","
				 << found[j].second.pt.theta << ","
				 << found[j].second.pt.phi/pi << ","
				 << cos(found[j].second.pt.theta) << endl << endl;
		}
	}
	cout << endl << "Radius," << radius << ",Count," << found.size(); 
}


int main(int64 argc, char* argv[])

{
 	srand(time(NULL));

	int testNum, numPoints, numQueries, maxDepth, numTrials,duplicateYN;
	int minIdx,maxIdx,order;
	std::string outFile,inFile,dataFile,queryFile;
	float density;
	double minPhi,maxPhi,minTheta,maxTheta,minRad,maxRad;

	testNum = atoi(argv[1]);

	cout << "Running Test Number: " << testNum << endl << endl;

	// Duplicate Point Insertion Test
	if( testNum == 0 ) {
		maxDepth = atoi(argv[2]);
		numPoints = atoi(argv[3]);
		duplicateYN = atoi(argv[4]);
		minPhi = atof(argv[5]);
		maxPhi = atof(argv[6]);
		minTheta = atof(argv[7]);
		maxTheta = atof(argv[8]);
		MRHDuplicatePointInsertionTest(maxDepth,numPoints,duplicateYN,minPhi,maxPhi,minTheta,maxTheta);
	}

	// Benchmark File I/O Test
	if( testNum == -1 ) {
		numPoints = atoi(argv[2]);
		numQueries = atoi(argv[3]);
		minPhi = atof(argv[4]);
		maxPhi = atof(argv[5]);
		minTheta = atof(argv[6]);
		maxTheta = atof(argv[7]);
		minRad = atof(argv[8]);
		maxRad = atof(argv[9]);
		TestFileWritersReaders(numPoints,numQueries,minPhi,maxPhi,minTheta,maxTheta,minRad,maxRad);
	}

	// Create Random Measurements File
	if( testNum == -2 ) {
		numPoints = atoi(argv[2]);
		outFile = argv[3];
		minPhi = atof(argv[4]);
		maxPhi = atof(argv[5]);
		minTheta = atof(argv[6]);
		maxTheta = atof(argv[7]);
		CreateRandomMeasurementsFile(numPoints,outFile,minPhi,maxPhi,minTheta,maxTheta);
	}

	// Create Random Measurements File of HPX order*density
	if( testNum == -3 ) {
		order = atoi(argv[2]);
		density = atof(argv[3]);
		outFile = argv[4];
		CreateRandomMeasurementsFileOfDensity(order,density,outFile);
	}

	// Create Random Query Discs File
	if( testNum == -4 ) {
		numQueries = atoi(argv[2]);
		outFile = argv[3];
		minPhi = atof(argv[4]);
		maxPhi = atof(argv[5]);
		minTheta = atof(argv[6]);
		maxTheta = atof(argv[7]);
		minRad = atof(argv[8]);
		maxRad = atof(argv[9]);
		cout << "Num Queries: " << numQueries << endl;
		cout << "Outfile: " << outFile << endl;
		cout << "Min Phi: " << minPhi << endl;
		cout << "Max Phi: " << maxPhi << endl;
		cout << "Min Theta: " << minTheta << endl;
		cout << "Max Theta: " << maxTheta << endl;
		cout << "Min Radius: " << minRad << endl;
		cout << "Max Radius: " << maxRad << endl;
		CreateRandomDiscsFile(numQueries,outFile,minPhi,maxPhi,minTheta,maxTheta,minRad,maxRad);
	}

	// Create Random Query Polygons File
	if( testNum == -5 ) {
		numQueries = atoi(argv[2]);
		outFile = argv[3];
		minPhi = atof(argv[4]);
		maxPhi = atof(argv[5]);
		minTheta = atof(argv[6]);
		maxTheta = atof(argv[7]);
		minRad = atof(argv[8]);
		maxRad = atof(argv[9]);
		cout << "Num Queries: " << numQueries << endl;
		cout << "Outfile: " << outFile << endl;
		cout << "Min Phi: " << minPhi << endl;
		cout << "Max Phi: " << maxPhi << endl;
		cout << "Min Theta: " << minTheta << endl;
		cout << "Max Theta: " << maxTheta << endl;
		cout << "Min Radius: " << minRad << endl;
		cout << "Max Radius: " << maxRad << endl;
		CreateRandomPolysFile(numQueries,outFile,minPhi,maxPhi,minTheta,maxTheta,minRad,maxRad);
	}

	// Create Random Query Strips File
	if( testNum == -6 ) {
		numQueries = atoi(argv[2]);
		outFile = argv[3];
		double minTheta = atof(argv[4]);
		double maxTheta = atof(argv[5]);
		cout << "Num Queries: " << numQueries << endl;
		cout << "Outfile: " << outFile << endl;
		cout << "Min Theta: " << minTheta << endl;
		cout << "Max Theta: " << maxTheta << endl;
		CreateRandomStripsFile(numQueries,outFile,minTheta,maxTheta);
	}

	// Create Random Neighbors Query File
	if( testNum == -7 ) {
		numQueries = atoi(argv[2]);
		dataFile = argv[3];
		outFile = argv[4];
		cout << "Num Queries: " << numQueries << endl;
		cout << "Datafile: " << dataFile << endl;
		cout << "Outfile: " << outFile << endl;
		CreateRandomNeighborsFile(numQueries,dataFile,outFile);
	}

	// Process NOAA Weather Station Locations File
	if( testNum == -8 ) {
		dataFile = argv[2];
		outFile = argv[3];
		ProcessNOAA(dataFile,outFile);
	}

	// Process SDSSIII Query File
	if( testNum == -9 ) {
		dataFile = argv[2];
		outFile = argv[3];
		ProcessSDSS(dataFile,outFile);
	}

	// Process FRAG Flyout File
	if( testNum == -10 ) {
		dataFile = argv[2];
		outFile = argv[3];
		ProcessFRAG(dataFile,outFile);
	}

	// Process MOON Crater Positions File
	if( testNum == -11 ) {
		dataFile = argv[2];
		outFile = argv[3];
		ProcessMOON(dataFile,outFile);
	}

	// Load & Save: Do MRH Point Insertion then
	// save out MRH data structure into file for
	// benchmark tests.
	if( testNum == -12 ) {
		dataFile = argv[2];
		outFile = argv[3];
		MRHPointInsertionAndSave(dataFile,outFile);
	}

	// Benchmark Test: MRH Point Insertion
	if( testNum == 1 ) {
		dataFile = argv[2];
		outFile = argv[3];
		numTrials = atoi(argv[4]);
		MRHPointInsertionBenchMark(dataFile,outFile,numTrials);
	}
	
	// Benchmark Test: HPX Point Insertion
	if( testNum == 2 ) {
		dataFile = argv[2];
		outFile = argv[3];
		numTrials = atoi(argv[4]);
		HPXPointInsertionBenchMark(dataFile,outFile,numTrials);
	}

	// Benchmark Test: MRH Disc Query
	if( testNum == 3 ) {
		dataFile = argv[2];
		queryFile = argv[3];
		outFile = argv[4];
		numTrials = atoi(argv[5]);
		MRHDiscQueryBenchMark(dataFile,queryFile,outFile,numTrials);
	}

	// Benchmark Test: HPX Disc Query
	if( testNum == 4 ) {
		dataFile = argv[2];
		queryFile = argv[3];
		outFile = argv[4];
		numTrials = atoi(argv[5]);
		HPXDiscQueryBenchMark(dataFile,queryFile,outFile,numTrials);
	}

	// HPX vs MRH Disc Query Accuracy Test
	if( testNum == 5 ) {
		dataFile = argv[2];
		queryFile = argv[3];
		outFile = argv[4];
		CompareDiscQueryOutput(dataFile,queryFile,outFile);
	}

	// Benchmark Test: MRH Poly Query
	if( testNum == 6 ) {
		dataFile = argv[2];
		queryFile = argv[3];
		outFile = argv[4];
		numTrials = atoi(argv[5]);
		MRHPolyQueryBenchMark(dataFile,queryFile,outFile,numTrials);
	}

	// Benchmark Test: HPX Poly Query
	if( testNum == 7 ) {
		dataFile = argv[2];
		queryFile = argv[3];
		outFile = argv[4];
		numTrials = atoi(argv[5]);
		HPXPolyQueryBenchMark(dataFile,queryFile,outFile,numTrials);
	}

	// HPX vs MRH Poly Query Accuracy Test
	if( testNum == 8 ) {
		dataFile = argv[2];
		queryFile = argv[3];
		outFile = argv[4];
		ComparePolyQueryOutput(dataFile,queryFile,outFile);
	}

	// Benchmark Test: MRH Strip Query
	if( testNum == 9 ) {
		dataFile = argv[2];
		queryFile = argv[3];
		outFile = argv[4];
		numTrials = atoi(argv[5]);
		MRHStripQueryBenchMark(dataFile,queryFile,outFile,numTrials);
	}

	// Benchmark Test: HPX Strip Query
	if( testNum == 10 ) {
		dataFile = argv[2];
		queryFile = argv[3];
		outFile = argv[4];
		numTrials = atoi(argv[5]);
		HPXStripQueryBenchMark(dataFile,queryFile,outFile,numTrials);
	}

	// HPX vs MRH Strip Query Accuracy Test
	if( testNum == 11 ) {
		dataFile = argv[2];
		queryFile = argv[3];
		outFile = argv[4];
		CompareStripQueryOutput(dataFile,queryFile,outFile);
	}

	// Benchmark Test: MRH Neighbor Query
	if( testNum == 12 ) {
		dataFile = argv[2];
		queryFile = argv[3];
		outFile = argv[4];
		numTrials = atoi(argv[5]);
		MRHNeighborQueryBenchMark(dataFile,queryFile,outFile,numTrials);
	}

	// Benchmark Test: HPX Neighbor Query
	if( testNum == 13 ) {
		dataFile = argv[2];
		queryFile = argv[3];
		outFile = argv[4];
		numTrials = atoi(argv[5]);
		HPXNeighborQueryBenchMark(dataFile,queryFile,outFile,numTrials);
	}

	// HPX vs MRH Neighbor Query Accuracy Test
	if( testNum == 14 ) {
		dataFile = argv[2];
		queryFile = argv[3];
		outFile = argv[4];
		CompareNeighborQueryOutput(dataFile,queryFile,outFile);
	}

	// Two-Point Correlation Binning
	if( testNum == 15 ) {
		dataFile = argv[2];
		float radius = atof(argv[3]);
		int verbose = atoi(argv[4]);
		TwoPointCorrBinning(dataFile,radius,verbose);
	}
}

