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

//GeographicLib References
#include <GeographicLib/Gnomonic.hpp>

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

void GenRandomData(int numPoints,std::vector<Measurement>& ms)
{
	int MIN_INT = 0; int MAX_INT = 9999;
	double MIN_DOUBLE = -9999.0; double MAX_DOUBLE = 9999.0;
	Measurement m;
	ms.clear();
	int64 recNum = 0;
	for(int i = 0; i < numPoints; i++) {
		m.rec = recNum;
		m.pt.phi = 0.0; //assign later
		m.pt.theta = 0.0; //assign later
		m.data1 = MIN_INT + static_cast <int> (rand()) /( static_cast <double> (RAND_MAX/(MAX_INT-MIN_INT)));
		m.data2 = MIN_INT + static_cast <int> (rand()) /( static_cast <double> (RAND_MAX/(MAX_INT-MIN_INT)));
		m.data3 = MIN_INT + static_cast <int> (rand()) /( static_cast <double> (RAND_MAX/(MAX_INT-MIN_INT)));
		m.data4 = MIN_INT + static_cast <int> (rand()) /( static_cast <double> (RAND_MAX/(MAX_INT-MIN_INT)));

		m.data1 = MIN_DOUBLE + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(MAX_DOUBLE-MIN_DOUBLE)));
		m.data2 = MIN_DOUBLE + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(MAX_DOUBLE-MIN_DOUBLE)));
		m.data3 = MIN_DOUBLE + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(MAX_DOUBLE-MIN_DOUBLE)));
		m.data4 = MIN_DOUBLE + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(MAX_DOUBLE-MIN_DOUBLE)));
		
		ms.push_back(m);		
		recNum += 1;
	}
}

void GenRandomPhiTheta(int numPoints,std::vector<pointing>& points,int64& order,bool verbose)
{
	// GIS Longitude (Phi) Range: [-180.0,+180.0] degrees, relative to Prime Meridian.
	// GIS Latitude (Theta) Range: [-90.0,+90.0] degrees, relative to Equator.
	// HPX Longitude (Phi) Range: [0.0,360.0] degrees, relative to Prime Meridian and progressing East.
	// HPX Colatitude (Theta) Range: [0.0,180.0] degrees, relative to North Pole and progressing South to
	// South Pole.
	double LO_THETA = 0.001*D2R; double LO_PHI = 0.001*D2R;
	double HI_THETA = 179.999*D2R; double HI_PHI = 359.999*D2R;
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

	while( count < numPoints )
	{
		p1.phi = LO_PHI + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(HI_PHI-LO_PHI)));
		p1.theta = LO_THETA + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(HI_THETA-LO_THETA)));

		// Compute the minimum point seperaton
		for( j = 0; j < points.size(); j++ ) {
			p2 = points[j];
			curRad = acos(fabs(cosdist_zphi(cos(p1.theta),p1.phi,cos(p2.theta),p2.phi)));
			if(curRad < minRad) {
				minRad = curRad;
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
	}

	// Select appropriate HEALPix resolution based on Minimal Radial Distance
	order = -1;
	for( i = 0; i <= 29; i++ ) 
	{
		if( HpxOrderResTable[i] < minRad )
		{
			order = i;
			break;
		}
	}
	if(verbose) {
		printf("For random dataset computed: HPX order = %d HPX cell res = %2.15f\n\n",order,minRad);
	}

}

void GenRandomPhiThetaFiltered(int numPoints,float minDist,std::vector<pointing>& points,int64& order,bool verbose)
{
	// GIS Longitude (Phi) Range: [-180.0,+180.0] degrees, relative to Prime Meridian.
	// GIS Latitude (Theta) Range: [-90.0,+90.0] degrees, relative to Equator.
	// HPX Longitude (Phi) Range: [0.0,360.0] degrees, relative to Prime Meridian and progressing East.
	// HPX Colatitude (Theta) Range: [0.0,180.0] degrees, relative to North Pole and progressing South to
	// South Pole.
	double LO_THETA = 0.001*D2R; double LO_PHI = 0.001*D2R;
	double HI_THETA = 179.999*D2R; double HI_PHI = 359.999*D2R;
	pointing p1,p2;
	double minRad,curRad;
	unsigned int i,j;
	int64 nside;
	double pixres;
	bool passMinDist;
	int count = 0;
	int reattempts = 1000;
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
	while( count < numPoints )
	{
		p1.phi = LO_PHI + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(HI_PHI-LO_PHI)));
		p1.theta = LO_THETA + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(HI_THETA-LO_THETA)));

		// Check new point for minimum distance between ALL OTHER points requirement
		passMinDist = true;
		for( j = 0; j < points.size(); j++ ) {
			p2 = points[j];
			curRad = acos(fabs(cosdist_zphi(cos(p1.theta),p1.phi,cos(p2.theta),p2.phi)));
			if(curRad < minDist) {
				passMinDist = false;
				j = points.size();
			} else {
				if(curRad < minRad) {
					minRad = curRad;
				}
			}
		}
		if( passMinDist == true ) {
			points.push_back(p1);
			if(verbose) {
				cout << count+1 << "," 
					 << p1.phi << "," 
					 << p1.theta << "," 
					 << p1.phi/pi << "," 
					 << cos(p1.theta) << endl;
			}
			count++;
			reattempts = 1000; //reset the reattempts counter
		} else {
			reattempts--;
			if(reattempts == 0) {
				cout << "Unable to create new random point outside minimum point distance specified!\n";
				cout << "Specify larger minimum point distance or reduce total number of points! EXITING!\n";
				exit(1);
			}
		}
	}

	// Select appropriate HEALPix resolution based on Minimal Radial Distance
	order = -1;
	for( i = 0; i <= 29; i++ ) 
	{
		if( HpxOrderResTable[i] < minRad )
		{
			order = i;
			break;
		}
	}
	if(verbose) {
		printf("For random dataset computed: HPX order = %d HPX cell res = %2.15f\n\n",order,minRad);
	}

}

void RandomDisc(double& phi, double& theta, double& radius){
	// GIS Longitude (Phi) Range: [-180.0,+180.0] degrees, relative to Prime Meridian.
	// GIS Latitude (Theta) Range: [-90.0,+90.0] degrees, relative to Equator.
	// HPX Longitude (Phi) Range: [0.0,360.0] degrees, relative to Prime Meridian and progressing East.
	// HPX Colatitude (Theta) Range: [0.0,180.0] degrees, relative to North Pole and progressing South to
	// South Pole.
	pointing pt;
	bool validDisc = false;
	double LO_THETA = 0.01*D2R; double LO_PHI = 0.01*D2R;
	double HI_THETA = 179.99*D2R; double HI_PHI = 359.99*D2R;
	double MIN_RAD = 1.0*D2R; double MAX_RAD = 45.0*D2R;
	double THETA_BUFF = 20.0*D2R;
	double PHI_BUFF = 20.0*D2R;

	while( validDisc == false ) {
		// Generate random query center location in HPX Coordinate System
		phi = LO_PHI + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(HI_PHI-LO_PHI)));
		theta = LO_THETA + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(HI_THETA-LO_THETA)));

		// Generate random query radius, making sure radius doesn't sweep
		// off of allowed latitude,longitude limits!
		radius = MIN_RAD + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(MAX_RAD-MIN_RAD)));

		if( phi+radius < HI_PHI-PHI_BUFF &&
			phi-radius > LO_PHI+PHI_BUFF &&
			theta+radius < HI_THETA-THETA_BUFF &&
			theta-radius > LO_THETA+THETA_BUFF ) {
			validDisc = true;
		}
	}
}

vector<pointing>  RandomConvexPoly(bool RANDPOLYVERBOSE) {
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
	RandomDisc(qPhi,qTheta,qRadius);
	
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
	int i,j,k;
	Matches = 0;
	HPXUnique = 0;
	MRHUnique = 0;

	if(verbose){
		cout << "\nHPX Found Measurement Indices:\n";
		cout << "#,REC,PHI,THETA,P,Z\n";
	}

	// Now compare each HPX reported data index to list of
	// those found by MRH. Will get number of HPX-MRH matches
	// as well as unique HPX found.
	for( j = 0; j < foundHPX.size(); j++ ) {
		matchedYN = false;
		for( k = 0; k < foundMRH.size(); k++ ) {
			// Check for match
			if( foundMRH[k].equals(foundHPX[j]) ) {
				Matches++;
				matchedYN = true;
			}
		}
		// Measurement index found in HPX but NOT MRH!
		if( matchedYN == false ) {
			HPXUnique++;
		}
		if(verbose){
			// Output the HPX poly query found data points
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
				matchedYN = true;
			}
		}
		// Measurement index found in HPX but NOT MRH!
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

void TestMortonSort(int maxDepth,int numPoints)
{
	uint64 start_s,stop_s,timeInsertMRH;
	int i,nMRHnodes = 0;
	int64 order;
	pointing p1,p2;
	vec3 v1,v2;
	vector<pointing> points;
	vector<Measurement> measurements;

	MultiResHpx_Map<Measurement> mMRH(maxDepth,NEST);
		
	// Generate random GIS longitude, latitude pairs
	GenRandomPhiTheta(numPoints,points,order,false);

	// Generate random measurement data
	GenRandomData(numPoints,measurements);

	// Insert the GIS longitude, latitude pairs into MRH data structure
	cout << "Begin insertion of points into MRH data structure..." << endl;
	
	start_s=GetTimeMs64();
	for( i = 0; i < numPoints; i++)
	{
		// Insert next GIS longitude,latitude, data index tuple into MRH
		measurements[i].pt = points[i]; // Set spatial location of measurment
		mMRH.AddRecord(measurements[i],points[i]);
		if( i%1000 == 0 ) {
			cout << "Inserting Point " << i+1 << " of " << numPoints << endl;
		}
	}
	stop_s=GetTimeMs64();
	timeInsertMRH = stop_s-start_s;
	cout << "End insertion of points into MRH data structure in " << timeInsertMRH
			 << " ms." << endl;

	//Print all trees as simple lists
	for( i = 0; i < 12; i++ ) 
	{
		cout << "\n\n######################\n";
		cout << "#### BASE CELL " << i << " ####\n";
		cout << "######################" << endl;
		mMRH.PrintListAtIndex(i); 
	}

}


//####################################################################################
//####################################################################################
//######																		######
//######   MRH AND HPX DATA POINT INSERTION AND MEMORY FOOT PRINT BENCHMARKING	######
//######																		######
//####################################################################################
//####################################################################################

void MRHRandomPointInsertion(int maxDepth,int numPoints,float minDist)
{
	uint64 start_s,stop_s,timeInsertMRH;
	int i,nMRHnodes = 0;
	int64 order;
	pointing p1,p2;
	vec3 v1,v2;
	vector<pointing> points;
	vector<Measurement> measurements;

	MultiResHpx_Map<Measurement> mMRH(maxDepth,NEST);
		
	// Generate random GIS longitude, latitude pairs
	GenRandomPhiTheta(numPoints,points,order,false);

	// Generate random measurement data
	GenRandomData(numPoints,measurements);

	// Insert the GIS longitude, latitude pairs into MRH data structure
	cout << "Begin insertion of points into MRH data structure..." << endl;
	
	start_s=GetTimeMs64();
	for( i = 0; i < numPoints; i++)
	{
		// Insert next GIS longitude,latitude, data index tuple into MRH
		measurements[i].pt = points[i]; // Set spatial location of measurment
		mMRH.AddRecord(measurements[i],points[i]);
		if( i%1000 == 0 ) {
			cout << "Inserting Point " << i+1 << " of " << numPoints << endl;
		}
	}
	stop_s=GetTimeMs64();
	timeInsertMRH = stop_s-start_s;
	cout << "End insertion of points into MRH data structure in " << timeInsertMRH
			 << " ms." << endl;

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
	for( i = 0; i < mMap.size(); i++ ) {
		cout << i << ","
			<< mMap[i].rec << ","
			<< mMap[i].pt.phi << ","
			<< mMap[i].pt.theta << ","
			<< mMap[i].pt.phi/pi << ","
			<< cos(mMap[i].pt.theta) << endl;
	}


}


void RandomPointInsertionBenchMark(int maxDepth,int numPoints,float minDist,int numTrials,std::string trialsOut)
{
	int maxDepthMRH,minDepthMRH,maxDepthHPX,minDepthHPX;
	uint64 start_s,stop_s,timeInsertMRH,timeInsertHPX;
	double avgDepthMRH,avgDepthHPX;
	int64 memFootPrintMRH,memFootPrintHPX;
	int HPXcollision,i,k,nMRHnodes,nHPXnodes = 0;
	int64 order,hpxid,sizeOfMRHNode,sizeOfHPXNode;
	bool HPXOrderNotSet = true;
	int numAttempts = 0;
	pointing p1,p2;
	vec3 v1,v2;
	vector<pointing> points;
	vector<Measurement> measurements;

	// Open trials output file
	FILE* fp = fopen(trialsOut.c_str(), "w");
	if(fp == NULL)
	{
	  cout << "Unable to open trials output file: " << trialsOut.c_str() << endl;
	  return;        
	}
	// Write out header
    fprintf(fp,"COUNT,MINDIST,,MRHNODES,MRHNODESIZE,MRHMINDEPTH,MRHMAXDEPTH,MRHAVGDEPTH,MRHFOOTPRINT,MRHTIMEINSERT,,");
    fprintf(fp,"HPXNODES,HPXNODESIZE,HPXMINDEPTH,HPXMAXDEPTH,HPXAVGDEPTH,HPXFOOTPRINT,HPXTIMEINSERT,HPXCOLLISION,");
    fprintf(fp,"\n");

	for( k = 0; k < numTrials; k++ )
	{
		cout << "\nRun " << k+1 << " of " << numTrials << endl;

		//// Instantiate MRH_Map data structure
		MultiResHpx_Map<Measurement> mMRH(maxDepth,NEST);
		
		//// Instantiate HPX_Map data structure
		Healpix_Map<Measurement> mHPX(1,NEST);
			
		// Generate random GIS longitude, latitude pairs
		GenRandomPhiTheta(numPoints,points,order,false);

		// Generate random data
		GenRandomData(numPoints,measurements);

		// Insert the GIS longitude, latitude pairs into MRH data structure
		if(TA1_DEBUG) {
			cout << "Begin insertion of points into MRH data structure..." << endl;
		}
		start_s=GetTimeMs64();
		for( i = 0; i < numPoints; i++)
		{
			// Insert next GIS longitude,latitude, data index tuple into MRH
			measurements[i].pt = points[i]; // Set spatial location of measurment
			mMRH.AddRecord(measurements[i],points[i]);
			if(TA1_DEBUG) {
				if( i%1000 == 0 ) {
					cout << "Inserting Point " << i+1 << " of " << numPoints << endl;
				}
			}
		}
		stop_s=GetTimeMs64();
		timeInsertMRH = stop_s-start_s;
		if(TA1_DEBUG) {
			cout << "End insertion of points into MRH data structure in " << timeInsertMRH
				 << " ms." << endl;
		}
		// Insert the GIS longitude, latitude pairs into HEALPix data structure
		HPXcollision = 0;
		mHPX.Set(order,NEST);
		if(TA1_DEBUG) {
			cout << "Begin insertion of points into HPX data structure..." << endl;
		}
		start_s=GetTimeMs64();
		for( i = 0; i < numPoints; i++)
		{
			// Compute HEALPix index of next GIS longitude,latitude,data index tuple
			hpxid = mHPX.ang2pix(points[i]);

			// Insert data index associated with
			measurements[i].pt = points[i];// Set spatial location of measurment
			mHPX[hpxid] = measurements[i]; 
		}
		stop_s=GetTimeMs64();
		timeInsertHPX = stop_s-start_s;
		if(TA1_DEBUG) {
			cout << "End insertion of points into HPX data structure in " << timeInsertHPX
				 << " ms." << endl;
		}
		if(TA1_DEBUG) {
			cout << "Time Insert Nodes (ms): MRH = " << timeInsertMRH
				 << " HPX = " << timeInsertHPX << endl;
		}
		// Compute statistics of the final data structures 

		// Compute memory footprint of each data structure
		nMRHnodes = mMRH.NumNodes();
		nHPXnodes = mHPX.Npix();
		sizeOfMRHNode = sizeof(MortonNode);
		sizeOfHPXNode = sizeof(Measurement);

		memFootPrintMRH = sizeOfMRHNode*nMRHnodes + sizeof(Measurement)*mMRH.NumRec();
		memFootPrintHPX = sizeOfHPXNode*nHPXnodes;

		// Get Min,Max,Average Depth of QuadTrees of each data stucture
		maxDepthMRH = mMRH.MaxDepth();
		minDepthMRH = mMRH.MinDepth();
		avgDepthMRH = mMRH.AvgDepth();

		  // Depth of HEALPix data structure is always the same = order, because
		  // HEALPix is always a full QuadTree forest deteremined at instantiation and
		  // is the size of the HEALPix Map or vector<T> created.
		minDepthHPX = order;
		maxDepthHPX = order;
		avgDepthHPX = order;

		if(TA1_DEBUG) {
			cout << "\nMRH stats: numPoints = " << numPoints << endl
				           << " numNodes = " << nMRHnodes << endl
				           << " sizeNode = " << sizeOfMRHNode << endl
				           << " minDepth = " << minDepthMRH << endl 
						   << " maxDepth = " << maxDepthMRH << endl 
						   << " avgDepth = " << avgDepthMRH << endl 
						   << " footprint = " << memFootPrintMRH << endl;
			cout << "\nHPX stats: numPoints = " << numPoints << endl
				           << " numNodes = " << nHPXnodes << endl
				           << " sizeNode = " << sizeOfHPXNode << endl
				           << " minDepth = " << minDepthHPX << endl 
						   << " maxDepth = " << maxDepthHPX << endl 
						   << " avgDepth = " << avgDepthHPX << endl 
						   << " footprint = " << memFootPrintHPX << endl;	
		}

		// Write trial result to file
		fprintf(fp,"%d,%d,%f,%d,%d,%f,%f,%d,%f,%d,%d,%f,%f,%d\n",numPoints,
			nMRHnodes,sizeOfMRHNode,minDepthMRH,maxDepthMRH,avgDepthMRH,memFootPrintMRH,timeInsertMRH,
			nHPXnodes,sizeOfHPXNode,minDepthHPX,maxDepthHPX,avgDepthHPX,memFootPrintHPX,timeInsertHPX,HPXcollision);
	
		//Print all trees
		for(int i = 0; i < 12; i++ ) 
		{
			cout << "\n\n######################\n";
			cout << "#### BASE CELL " << i << " ####\n";
			cout << "######################" << endl;
			mMRH.PrintTreeAtIndex(i); 
		}
	}
	// End of trials, close file
	fclose(fp); 
}

//####################################################################################
//####################################################################################
//######																		######
//######			MRH AND HPX DATA RANGE QUERY BENCHMARKING					######
//######																		######
//####################################################################################
//####################################################################################


void RandomDiscQueryBenchMark(int maxDepth,int numPoints,float minDist,int numQueries,int numTrials,std::string trialsOut)
{
	pointing pt;
	std::vector<pointing> points,points2;
	uint64 start_s,stop_s;
	int HPXcollision,i,j,nTrial,nMRHnodes,nHPXnodes = 0;
	int totalMrhFound,totalHpxFound,totalMatches;
	int64 hpxid,order;
	double avgQueryTimeMRH,avgQueryTimeHPX;
	double qPhi,qTheta,qRadius;
	pointing p1,p2;
	vec3 v1,v2;
	vector<double> dataRad;
	vector<double> HpxOrderResTable;
	rangeset<int64> pixset;
	std::vector<int64> foundV; 
	vector<Measurement> measurements;
	std::vector<Measurement> foundMRH;
	std::vector<Measurement> foundHPX; 

	// Open trials output file
	FILE* fp = fopen(trialsOut.c_str(), "w");
	if(fp == NULL)
	{
	  cout << "Unable to open trials output file: " << trialsOut.c_str() << endl;
	  return;        
	}
	// Write out header
	fprintf(fp,"COUNT,NUMQUERIES,MRHNODES,MRHAVGTIMEQ,HPXNODES,HPXAVGTIMEQ\n");

	totalMrhFound = 0;
	totalHpxFound = 0;
	totalMatches = 0;
	for( nTrial = 0; nTrial < numTrials; nTrial++ )
	{

		cout << "\nRun " << nTrial+1 << " of " << numTrials << endl;
		HPXcollision = 0;
		//// Instantiate MRH_Map data structure
		MultiResHpx_Map<Measurement> mMRH(maxDepth,NEST);
		
		//// Instantiate HPX_Map data structure
		Healpix_Map<Measurement> mHPX(1,NEST);
		
		// Generate random GIS longitude, latitude pairs
		GenRandomPhiTheta(numPoints,points,order,false);

		// Generate random data
		GenRandomData(numPoints,measurements);

		// Insert the GIS longitude, latitude pairs into MRH & HPX data structures
		mHPX.Set(order,NEST);
		for( i = 0; i < numPoints; i++)
		{
			// Insert next GIS longitude,latitude, data index tuple into MRH
			measurements[i].pt = points[i];  // Set spatial location of measurment
			mMRH.AddRecord(measurements[i],points[i]);

			// Compute HEALPix index of next GIS longitude,latitude,data index tuple
			hpxid = mHPX.ang2pix(points[i]);
			measurements[i].pt = points[i];// Set spatial location of measurment
			mHPX[hpxid] = measurements[i]; 
		}
		cout << "	Loaded GIS Points into MRH and HPX data structures!\n";

		// Now generate inputs for numQueries of random disc queries on data
		points.clear();
		dataRad.clear();
		for( i = 0; i < numQueries; i++ )
		{
			RandomDisc(qPhi,qTheta,qRadius);
			pt.phi = qPhi;
			pt.theta = qTheta;
			points.push_back(pt);
            dataRad.push_back(qRadius);
		}
		cout << "	Generated " << numQueries << " Disc Query Inputs!\n";

		// Compute MRH QueryDisc Benchmark
		start_s=GetTimeMs64();
		for( i = 0; i < numQueries; i++)
		{          
			foundMRH = mMRH.QueryDisc(points[i],dataRad[i]);
		}
		stop_s=GetTimeMs64();
        avgQueryTimeMRH = double(stop_s-start_s)/double(numQueries);
		cout << "	Completed " << numQueries << " MRH Disc Queries in " << double(stop_s-start_s) 
			<< " ms! Average Query: " << avgQueryTimeMRH << "\n";

		// Compute HPX query_disc Benchmark
		start_s=GetTimeMs64();
		for( i = 0; i < numQueries; i++)
		{					
			mHPX.query_disc(points[0],dataRad[0],pixset);
			
			// Convert cell rangeset to list of individual
			// cells.
			pixset.toVector(foundV);

			// Spin through HPX's found indices keep only the
			// non empty indices that pass the point-in-disc test
			for( j=0; j < foundV.size(); j++ ) {
				if( mHPX[foundV[j]].rec != EMPTY ) {
					if( IsPointInDisc(points[i],dataRad[i],mHPX[foundV[j]].pt) ){
						foundHPX.push_back(mHPX[foundV[j]]);
					}
				} 
			}
		}
        stop_s=GetTimeMs64();
        avgQueryTimeHPX = double(stop_s-start_s)/double(numQueries);
		cout << "	Completed " << numQueries << " HPX Disc Queries in " << double(stop_s-start_s) 
			<< " ms! Average Query: " << avgQueryTimeHPX << "\n";

		nMRHnodes = mMRH.NumNodes();
		nHPXnodes = mHPX.Npix();

		// Write trial result to file
  		fprintf(fp,"%d,%d,%d,%f,%d,%f\n",
			numPoints,numQueries,nMRHnodes,avgQueryTimeMRH,nHPXnodes,avgQueryTimeHPX);
		cout << "	Wrote Trial Results to " << trialsOut.c_str() << "!\n";

	}
	// End of trials, close file
	fclose(fp); 

	// Finished!
	cout << "Trials Complete!\n\n";
}


void RandomDiscQueryAccuracyBenchMark(int maxDepth,int numPoints,float minDist,int numQueries,int numTrials,std::string trialsOut)
{
	int HPXcollision,i,j,nTrial,nMRHnodes,nHPXnodes = 0;
	int	totalMrhFound,totalHpxFound,totalMatches,totalMrhUnique,totalHpxUnique;

	int64 hpxid,order;
	double qPhi,qTheta,qRadius;
	pointing p1,p2;
	pointing qCent;
	pointing pt;

	vec3 v1,v2;
	vector<double> dataRad;
	vector<pointing> points,points2;
	vector<double> HpxOrderResTable;
	rangeset<int64> pixset;
	std::vector<int64> foundV;
	int ptsFoundHPX,ptsFoundMRH,ptsFound,numMatches;
	std::vector<Measurement> foundMRH;
	std::vector<Measurement> foundHPX; 
	bool matchedYN = false;

	int Matches = 0;
	int MRHUnique = 0;
	int HPXUnique = 0;

	double topZ,bottomZ,leftP,rightP;
	std::vector<Measurement> measurements;

	// Open trials output file
	FILE* fp = fopen(trialsOut.c_str(), "w");
	if(fp == NULL)
	{
	  cout << "Unable to open trials output file: " << trialsOut.c_str() << endl;
	  return;        
	}
	// Write out header
	fprintf(fp,"COUNT,NUMQUERIES,NUMMATCHES,MRHNODES,MRHMAXDEPTH,MRHMINDEPTH,MRHAVGDEPTH,MRHFOUNDTOTAL,MRHUNIQUE,HPXNODES,HPXDEPTH,HPXFOUNDTOTAL,HPXUNIQUE\n");

	totalMrhFound = 0;
	totalHpxFound = 0;
	for( nTrial = 0; nTrial < numTrials; nTrial++ )
	{

		cout << "\nRun " << nTrial+1 << " of " << numTrials << endl;
		HPXcollision = 0;
		//// Instantiate MRH_Map data structure
		MultiResHpx_Map<Measurement> mMRH(maxDepth,NEST);
		
		//// Instantiate HPX_Map data structure
		Healpix_Map<Measurement> mHPX(1,NEST);
		
		// Generate random GIS longitude, latitude pairs
		GenRandomPhiTheta(numPoints,points,order,true);

		// Generate random data
		GenRandomData(numPoints,measurements);

		// Insert the GIS longitude, latitude pairs into MRH & HPX data structures
		mHPX.Set(order,NEST);
		for( i = 0; i < numPoints; i++)
		{
			// Insert next GIS longitude,latitude, data index tuple into MRH
			measurements[i].pt = points[i];  // Set spatial location of measurment
			mMRH.AddRecord(measurements[i],points[i]);

			// Compute HEALPix index of next GIS longitude,latitude,data index tuple
			hpxid = mHPX.ang2pix(points[i]);
			measurements[i].pt = points[i];// Set spatial location of measurment
			mHPX[hpxid] = measurements[i]; 
		}
		cout << "	Loaded GIS Points into MRH and HPX data structures!\n";

		// Now generate inputs for numQueries of random disc queries on data
		points.clear();
		dataRad.clear();
		for( i = 0; i < numQueries; i++ )
		{
			RandomDisc(qPhi,qTheta,qRadius);
			pt.phi = qPhi;
			pt.theta = qTheta;
			points.push_back(pt);
            dataRad.push_back(qRadius);
			if(DISC_ACC_DEBUG){
				topZ = cos(pt.theta+qRadius);
				bottomZ = cos(pt.theta-qRadius);
				leftP = (pt.phi-qRadius)/pi;
				rightP = (pt.phi+qRadius)/pi;
				cout << "Disc Boundaries:" << endl;
				cout << leftP << "," << topZ << endl;
				cout << rightP << "," << topZ << endl;
				cout << rightP << "," << bottomZ << endl;
				cout << leftP << "," << bottomZ << endl;
				cout << leftP << "," << topZ << endl << endl;
			}
		}
		cout << "	Generated " << numQueries << " Disc Query Inputs!\n";

		// Compute MRH to HPX Query Disc Accuracy
		totalMrhFound = 0;totalHpxFound=0;totalMatches=0;totalMrhUnique=0;totalHpxUnique=0;
		for( i = 0; i < numQueries; i++)
		{ 
			cout << "	Analyzing Query " << i+1 << " of " << numQueries << endl;
			ptsFoundMRH = 0;
			ptsFoundHPX = 0;
			numMatches = 0;
			ptsFound = 0;
			foundMRH.clear();
			foundHPX.clear();
			pixset.clear();
			foundV.clear();

			qCent = points[i];
			
			// Do MRH Disc Query
	 		foundMRH = mMRH.QueryDisc(points[i],dataRad[i]);

			// Do HPX Disc Query 
			mHPX.query_disc(points[i],dataRad[i],pixset);

			// Get count of data points found by MRH query
			ptsFoundMRH = foundMRH.size();

			// Convert cell rangeset to list of individual
			// cells.
			pixset.toVector(foundV);

			// First let's spin through HPX's found indices keep only the
			// non empty indices that pass the point-in-poly test
			for( j=0; j < foundV.size(); j++ ) {
				if( mHPX[foundV[j]].rec != EMPTY ) {
					// Apply point-in-disc test to filter out possible
					// false positives
					if( IsPointInDisc(points[i],dataRad[i],mHPX[foundV[j]].pt) ){
						foundHPX.push_back(mHPX[foundV[j]]);
					}
				} 
			}
			AnalyzeResults(foundHPX,foundMRH,Matches,HPXUnique,MRHUnique,true);
			totalMatches += Matches;
			totalMrhUnique += MRHUnique;
			totalHpxUnique += HPXUnique;
			totalMrhFound += foundMRH.size();
			totalHpxFound += foundHPX.size();

		}
		cout << "	Finished Computing MRH vs HPX Disc Query Accuracy Test!\n";

		// Compute/get statistics of the final data structures 
		nMRHnodes = mMRH.NumNodes();
		nHPXnodes = mHPX.Npix();
  	    // Get Min,Max,Average Depth of QuadTrees of each data stucture
		int maxDepthMRH = mMRH.MaxDepth();
		int minDepthMRH = mMRH.MinDepth();
		int avgDepthMRH = mMRH.AvgDepth();

		// Write trial result to file
	    //"COUNT,NUMQUERIES,NUMMATCHES,MRHNODES,MRHMAXDEPTH,MRHMINDEPTH,MRHAVGDEPTH,MRHFOUNDTOTAL,MRHUNIQUE,HPXNODES,HPXDEPTH,HPXFOUNDTOTAL,HPXUNIQUE
		fprintf(fp,"%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",
			numPoints,numQueries,totalMatches,
			nMRHnodes,maxDepthMRH,minDepthMRH,avgDepthMRH,totalMrhFound,totalMrhUnique,
			nHPXnodes,order,totalHpxFound,totalHpxUnique);
		//DEBUG
		cout << "	Wrote Trial Results to " << trialsOut.c_str() << "!\n";

	}
	// End of trials, close file
	fclose(fp); 

	// Finished!
	cout << "Trials Complete!\n\n";

}


//
//
// USER DEFINED POLYGON VS RANDOM DATA BENCH MARK
//
//

void RandomPolyQueryBenchMark(int maxDepth,int numPoints,float minDist,int numQueries,int numTrials,std::string trialsOut)
{

	uint64 start_s,stop_s;
	int i,j,nTrial,nMRHnodes,nHPXnodes = 0;
	int totalMrhFound,totalHpxFound,totalMatches;
	int64 hpxid,order;
	double avgQueryTimeMRH,avgQueryTimeHPX;
	pointing p1,p2,pt;
	vec3 v1,v2;
	vector<pointing> poly,dataRad,points,tempPoly,nextPoly;
	vector<vector<pointing>> polys,polys2;
	vector<double> HpxOrderResTable;
	rangeset<int64> pixset;
	std::vector<int64> foundV; 
	vector<MortonNode> mMrh;
    bool HPXOrderNotSet = true;
	std::vector<Measurement> foundMRH,foundHPX;
	std::vector<Measurement> measurements;

	// Open trials output file
	FILE* fp = fopen(trialsOut.c_str(), "w");
	if(fp == NULL)
	{
	  cout << "Unable to open trials output file: " << trialsOut.c_str() << endl;
	  return;        
	}
	// Write out header
	fprintf(fp,"COUNT,NUMQUERIES,MRHNODES,MRHAVGTIMEQ,HPXNODES,HPXAVGTIMEQ\n");

	totalMrhFound = 0;
	totalHpxFound = 0;
	totalMatches = 0;
	for( nTrial = 0; nTrial < numTrials; nTrial++ )
	{

		cout << "\nTrial " << nTrial+1 << " of " << numTrials << endl;

		//// Instantiate MRH_Map data structure
		MultiResHpx_Map<Measurement> mMRH(maxDepth,NEST);
		
		//// Instantiate HPX_Map data structure
		Healpix_Map<Measurement> mHPX(1,NEST);
		
		// Generate random GIS longitude, latitude pairs
		GenRandomPhiTheta(numPoints,points,order,false);

		// Generate random data
		GenRandomData(numPoints,measurements);

		// Insert the GIS longitude, latitude pairs into MRH & HPX data structures
		mHPX.Set(order,NEST);
		for( i = 0; i < numPoints; i++)
		{
			// Insert next GIS longitude,latitude, data index tuple into MRH
			measurements[i].pt = points[i];  // Set spatial location of measurment
			mMRH.AddRecord(measurements[i],points[i]);

			// Compute HEALPix index of next GIS longitude,latitude,data index tuple
			hpxid = mHPX.ang2pix(points[i]);
			measurements[i].pt = points[i];// Set spatial location of measurment
			mHPX[hpxid] = measurements[i]; 
		}
		if(POLY_DEBUG){
			cout << "	Loaded GIS Points into MRH and HPX data structures!\n";
		}

		// Generate random convex, ccw defined, polygons 
		int index = 0;
		
		for( i = 0; i < numQueries; i++ )
		{
			tempPoly.clear();
			nextPoly = RandomConvexPoly(false);

			// Store poly points for MRH queries
			polys.push_back(nextPoly);
		}

		// Compute MRH QueryPolygon Benchmark
		start_s=GetTimeMs64();
		for( i = 0; i < numQueries; i++)
		{          
 			mMRH.QueryPolygon(polys[i]);
		}
		stop_s=GetTimeMs64();
		avgQueryTimeMRH = double(stop_s-start_s)/double(numQueries);

		if(POLY_DEBUG){
			cout << "	Completed " << numQueries << " MRH Polygon Queries in " << double(stop_s-start_s) 
				<< " ms! Average Query: " << avgQueryTimeMRH << "\n";
		}
		//
		//	Compute HEALPIX QueryPolygon Benchmark
		//
		start_s=GetTimeMs64();
		for( i = 0; i < numQueries; i++)
		{					
			mHPX.query_polygon(polys[i],pixset);
			// Convert cell rangeset to list of individual
			// cells.
			pixset.toVector(foundV);

			// Spin through HPX's found indices keep only the
			// non empty indices that pass the point-in-poly test
			for( j=0; j < foundV.size(); j++ ) {
				if( mHPX[foundV[j]].rec != EMPTY ) {
					if( IsPointInPoly(mHPX[foundV[j]].pt,polys[i]) ){
						foundHPX.push_back(mHPX[foundV[j]]);
					}
				} 
			}
		}
		stop_s=GetTimeMs64();
		avgQueryTimeHPX = double(stop_s-start_s)/double(numQueries);
		if(POLY_DEBUG){
			cout << "	Completed " << numQueries << " HPX Polygon Queries in " << double(stop_s-start_s) 
				<< " ms! Average Query: " << avgQueryTimeHPX << "\n";
		}
     	nMRHnodes = mMRH.NumNodes();
		nHPXnodes = mHPX.Npix();

		// Write trial result to file
  		fprintf(fp,"%d,%d,%d,%f,%d,%f\n",
			numPoints,numTrials,nMRHnodes,avgQueryTimeMRH,nHPXnodes,avgQueryTimeHPX);
		if(POLY_DEBUG){
			cout << "	Wrote Trial Results to " << trialsOut.c_str() << "!\n";
		}
	}
	// End of trials, close file
	fclose(fp); 

	// Finished!
	cout << "Trials Complete!\n\n";

}


void RandomPolyQueryAccuracyBenchMark(int maxDepth,int numPoints,float minDist,int numQueries,int numTrials,std::string trialsOut)
{
	FILE* fp;
	int i,j,k,nQuery,nTrial,nMRHnodes,nHPXnodes = 0;
	int totalMrhFound,totalHpxFound,totalMatches,totalMrhUnique,totalHpxUnique;

	int64 hpxid,order;
	pointing pt;
	vec3 v1,v2;
	vector<double> dataRad;
	vector<pointing> points,tempPoly;
	vector<vector<pointing>> polys,polys2;
	vector<pointing> nextPoly;
	vector<double> HpxOrderResTable;
	rangeset<int64> pixset;
	std::vector<int64> foundV,foundVculled; 
	vector<MortonNode> found;
	std::vector<Measurement> foundMRH,foundHPX;
	std::vector<Measurement> measurements;
	bool matchedYN = false;
	int Matches = 0;
	int MRHUnique = 0;
	int HPXUnique = 0;

	// Open trials output file
	fp = fopen(trialsOut.c_str(), "w");
	if(fp == NULL)
	{
	  cout << "Unable to open trials output file: " << trialsOut.c_str() << endl;
	  return;        
	}
	// Write out header
	fprintf(fp,"COUNT,NUMQUERIES,MRHNODES,HPXNODES,MATCHES,MRHONLY,HPXONLY\n");

	totalMrhFound = 0;
	totalHpxFound = 0;
	totalMatches = 0;
	for( nTrial = 0; nTrial < numTrials; nTrial++ )
	{
		cout << "\nTrial " << nTrial+1 << " of " << numTrials << endl;

		//// Instantiate MRH_Map data structure
		MultiResHpx_Map<Measurement> mMRH(maxDepth,NEST);
		
		//// Instantiate HPX_Map data structure
		Healpix_Map<Measurement> mHPX(1,NEST);
		
		// Generate random GIS longitude, latitude pairs
		GenRandomPhiTheta(numPoints,points,order,true);

		// Generate random data
		GenRandomData(numPoints,measurements);

		// Insert the GIS longitude, latitude pairs into MRH & HPX data structures
		mHPX.Set(order,NEST);
		for( i = 0; i < numPoints; i++)
		{
			// Insert next GIS longitude,latitude, data index tuple into MRH
			measurements[i].pt = points[i];  // Set spatial location of measurment
			mMRH.AddRecord(measurements[i],points[i]);

			// Compute HEALPix index of next GIS longitude,latitude,data index tuple
			hpxid = mHPX.ang2pix(points[i]);
			measurements[i].pt = points[i];// Set spatial location of measurment
			mHPX[hpxid] = measurements[i]; 
		}
		if(POLY_DEBUG){
			cout << "	Loaded GIS Points into MRH and HPX data structures!\n";
		}

		// Generate random convex, cw defined, polygons 
		int index = 0;
		
		polys.clear();
		for( i = 0; i < numQueries; i++ )
		{
			nextPoly.clear();
			nextPoly = RandomConvexPoly(true);

			// Store poly points for MRH queries
			polys.push_back(nextPoly);
		}

		// Compute QueryPolygon Accuracy Benchmark
		totalMrhFound = 0;totalHpxFound=0;totalMatches=0;totalMrhUnique=0;totalHpxUnique=0;
		for( nQuery = 0; nQuery < numQueries; nQuery++)
		{          
			Matches = 0;
			MRHUnique = 0;
			HPXUnique = 0;
			foundV.clear();
			foundMRH.clear();
			foundHPX.clear();

			// MRH Polygon Query
			foundMRH = mMRH.QueryPolygon(polys[nQuery]);

			// HPX Polygon Query
			mHPX.query_polygon(polys[nQuery],pixset);
		
			// Convert cell rangeset to list of individual
			// cells.
			pixset.toVector(foundV);

			// Spin through HPX's found indices keep only the
			// non empty indices that pass the point-in-poly test
			for( j=0; j < foundV.size(); j++ ) {
				if( mHPX[foundV[j]].rec != EMPTY ) {
					if( IsPointInPoly(mHPX[foundV[j]].pt,polys[nQuery]) ){
						foundHPX.push_back(mHPX[foundV[j]]);
					}
				} 
			}

			AnalyzeResults(foundHPX,foundMRH,Matches,HPXUnique,MRHUnique,true);
			totalMatches += Matches;
			totalMrhUnique += MRHUnique;
			totalHpxUnique += HPXUnique;
			totalMrhFound += foundMRH.size();
			totalHpxFound += foundHPX.size();
		}

     	nMRHnodes = mMRH.NumNodes();
		nHPXnodes = mHPX.Npix();

		// Write trial result to file
  		fprintf(fp,"%d,%d,%d,%d,%d,%d,%d\n",
			numPoints,numTrials,nMRHnodes,nHPXnodes,Matches,MRHUnique,HPXUnique);
		if(POLY_ACC_DEBUG){
			cout << "	Wrote Trial Results to " << trialsOut.c_str() << "!\n";
		}
	}
	// End of trials, close file
	fclose(fp); 

	// Finished!
	cout << "Trials Complete!\n\n";
}


void RandomNeighborsQueryAccuracyBenchMark(int maxDepth,int numPoints,int numQueries,int numTrials,std::string trialsOut)
{
	int LO_IDX = 0; 
	int HI_IDX = numPoints-1; 
	int i,j,k,nQuery,nTrial,nMRHnodes,nHPXnodes = 0;
	int totalMrhFound,totalHpxFound,totalMatches,totalMrhUnique,totalHpxUnique;
	uint64 start_s,stop_s;
	int64 hpxid,order;
	pointing pt;
	vec3 v1,v2;
	vector<double> dataRad;
	vector<pointing> points,tempPoly;
	vector<int> neighbors;
	vector<pointing> nextPoly;
	vector<double> HpxOrderResTable;
	fix_arr<int64,8> pixset;
	std::vector<int64> foundV,foundVculled; 
	std::vector<Measurement> foundMRH,foundHPX;
	std::vector<Measurement> measurements;
	bool matchedYN = false;
	int Matches = 0;
	int MRHUnique = 0;
	int HPXUnique = 0;

	// Open trials output file
	FILE* fp = fopen(trialsOut.c_str(), "w");
	if(fp == NULL)
	{
	  cout << "Unable to open trials output file: " << trialsOut.c_str() << endl;
	  return;        
	}
	// Write out header
	fprintf(fp,"COUNT,NUMQUERIES,MRHNODES,MRHAVGTIMEQ,HPXNODES,HPXAVGTIMEQ\n");

	totalMrhFound = 0;
	totalHpxFound = 0;
	totalMatches = 0;
	for( nTrial = 0; nTrial < numTrials; nTrial++ )
	{

		cout << "\nTrial " << nTrial+1 << " of " << numTrials << endl;

		//// Instantiate MRH_Map data structure
		MultiResHpx_Map<Measurement> mMRH(maxDepth,NEST);
		
		//// Instantiate HPX_Map data structure
		Healpix_Map<Measurement> mHPX(1,NEST);
		
		// Generate random GIS longitude, latitude pairs
		GenRandomPhiTheta(numPoints,points,order,true);

		// Generate random data
		GenRandomData(numPoints,measurements);

		// Insert the GIS longitude, latitude pairs into MRH & HPX data structures
		mHPX.Set(order,NEST);
		for( i = 0; i < numPoints; i++)
		{
			// Insert next GIS longitude,latitude, data index tuple into MRH
			measurements[i].pt = points[i];  // Set spatial location of measurment
			mMRH.AddRecord(measurements[i],points[i]);

			// Compute HEALPix index of next GIS longitude,latitude,data index tuple
			hpxid = mHPX.ang2pix(points[i]);
			measurements[i].pt = points[i];// Set spatial location of measurment
			mHPX[hpxid] = measurements[i]; 
		}
		if(NEIGHBORS_DEBUG){
			cout << "	Loaded GIS Points into MRH and HPX data structures!\n";
		}

		// Compute random neighbor indices
		for( i = 0; i < numQueries; i++ ) {
			neighbors.push_back(LO_IDX + static_cast <int> (rand()) /( static_cast <double> (RAND_MAX/((double)HI_IDX-(double)LO_IDX))));
		}

		// Compute Neighbors Accuracy Benchmark
		totalMrhFound = 0;totalHpxFound=0;totalMatches=0;totalMrhUnique=0;totalHpxUnique=0;
		for( i = 0; i < numQueries; i++)
		{          
			Matches = 0;
			MRHUnique = 0;
			HPXUnique = 0;
			foundV.clear();
			foundMRH.clear();
			foundHPX.clear();

			if(NEIGHBORS_DEBUG){
				cout << "\n\nQuery Point:," << mHPX.ang2pix(measurements[neighbors[i]].pt) << ","
					<< order << ","
					<< measurements[neighbors[i]].pt.phi << ","
					<< measurements[neighbors[i]].pt.theta << ","
					<< measurements[neighbors[i]].pt.phi/pi << ","
					<< cos(measurements[neighbors[i]].pt.theta) << endl << endl;
			}
			// MRH Neighbor Query
			foundMRH = mMRH.Neighbors(measurements[neighbors[i]].pt,order);

			// HPX Neighbor Query
			mHPX.Set(order,NEST);
			hpxid = mHPX.ang2pix(measurements[neighbors[i]].pt);
			mHPX.neighbors(hpxid,pixset);
		
			// Spin through HPX's found indices keep only the
			// non empty indices that pass the point-in-strip test
			for( j=0; j < 8; j++ ) {
				if(NEIGHBORS_DEBUG){
					cout << j << "," << pixset[j] << endl;
				}
				if( mHPX[pixset[j]].rec != EMPTY ) {
					foundHPX.push_back(mHPX[pixset[j]]);
				} 
			}

			AnalyzeResults(foundHPX,foundMRH,Matches,HPXUnique,MRHUnique,true);
			totalMatches += Matches;
			totalMrhUnique += MRHUnique;
			totalHpxUnique += HPXUnique;
			totalMrhFound += foundMRH.size();
			totalHpxFound += foundHPX.size();
		}

     	nMRHnodes = mMRH.NumNodes();
		nHPXnodes = mHPX.Npix();

		// Write trial result to file
  		fprintf(fp,"%d,%d,%d,%d,%d,%d,%d\n",
			numPoints,numTrials,nMRHnodes,nHPXnodes,Matches,MRHUnique,HPXUnique);
		if(NEIGHBORS_DEBUG){
			cout << "	Wrote Trial Results to " << trialsOut.c_str() << "!\n";
		}
	}
	// End of trials, close file
	fclose(fp); 

	// Finished!
	cout << "Trials Complete!\n\n";
}


void RandomStripQueryAccuracyBenchMark(int maxDepth,int numPoints,int numQueries,int numTrials,std::string trialsOut)
{
	double LO_THETA = 0.01*D2R; 
	double HI_THETA = 179.99*D2R; 
	int i,j,k,nQuery,nTrial,nMRHnodes,nHPXnodes = 0;
	int totalMrhFound,totalHpxFound,totalMatches,totalMrhUnique,totalHpxUnique;
	uint64 start_s,stop_s;
	int64 hpxid,order;
	pointing pt;
	vec3 v1,v2;
	vector<double> dataRad;
	vector<pointing> points,tempPoly;
	vector<double> qTheta1,qTheta2;
	vector<pointing> nextPoly;
	vector<double> HpxOrderResTable;
	rangeset<int64> pixset;
	std::vector<int64> foundV,foundVculled; 
	std::vector<Measurement> foundMRH,foundHPX;
	std::vector<Measurement> measurements;
	bool matchedYN = false;
	int Matches = 0;
	int MRHUnique = 0;
	int HPXUnique = 0;

	// Open trials output file
	FILE* fp = fopen(trialsOut.c_str(), "w");
	if(fp == NULL)
	{
	  cout << "Unable to open trials output file: " << trialsOut.c_str() << endl;
	  return;        
	}
	// Write out header
	fprintf(fp,"COUNT,NUMQUERIES,MRHNODES,MRHAVGTIMEQ,HPXNODES,HPXAVGTIMEQ\n");

	totalMrhFound = 0;
	totalHpxFound = 0;
	totalMatches = 0;
	for( nTrial = 0; nTrial < numTrials; nTrial++ )
	{

		cout << "\nTrial " << nTrial+1 << " of " << numTrials << endl;

		//// Instantiate MRH_Map data structure
		MultiResHpx_Map<Measurement> mMRH(maxDepth,NEST);
		
		//// Instantiate HPX_Map data structure
		Healpix_Map<Measurement> mHPX(1,NEST);
		
		// Generate random GIS longitude, latitude pairs
		GenRandomPhiTheta(numPoints,points,order,true);

		// Generate random data
		GenRandomData(numPoints,measurements);

		// Insert the GIS longitude, latitude pairs into MRH & HPX data structures
		mHPX.Set(order,NEST);
		for( i = 0; i < numPoints; i++)
		{
			// Insert next GIS longitude,latitude, data index tuple into MRH
			measurements[i].pt = points[i];  // Set spatial location of measurment
			mMRH.AddRecord(measurements[i],points[i]);

			// Compute HEALPix index of next GIS longitude,latitude,data index tuple
			hpxid = mHPX.ang2pix(points[i]);
			measurements[i].pt = points[i];// Set spatial location of measurment
			mHPX[hpxid] = measurements[i]; 
		}
		if(STRIP_DEBUG){
			cout << "	Loaded GIS Points into MRH and HPX data structures!\n";
		}

		// Compute random strip ranges
		for( i = 0; i < numQueries; i++ ) {
			qTheta1.push_back(LO_THETA + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(HI_THETA-LO_THETA))));
			qTheta2.push_back(LO_THETA + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(HI_THETA-LO_THETA))));
		}

		// Compute QueryStrip Accuracy Benchmark
		totalMrhFound = 0;totalHpxFound=0;totalMatches=0;totalMrhUnique=0;totalHpxUnique=0;
		for( nQuery = 0; nQuery < numQueries; nQuery++)
		{          
			Matches = 0;
			MRHUnique = 0;
			HPXUnique = 0;
			foundV.clear();
			foundMRH.clear();
			foundHPX.clear();

			if(STRIP_DEBUG){
				cout << "\nTHETA1," << qTheta1[i] << endl;
				cout << "THETA2," << qTheta2[i] << endl << endl;
			}
			// MRH Polygon Strip
			foundMRH = mMRH.QueryStrip(qTheta1[nQuery],qTheta2[nQuery]);

			// HPX Polygon Query
			mHPX.Set(order,RING);
			mHPX.query_strip(qTheta1[nQuery],qTheta2[nQuery],true,pixset);
		
			// Convert cell rangeset to list of individual
			// cells.
			pixset.toVector(foundV);

			// Spin through HPX's found indices keep only the
			// non empty indices that pass the point-in-strip test
			for( j=0; j < foundV.size(); j++ ) {
				// Convert RING index to NESTED index
				hpxid = mHPX.ring2nest(foundV[j]);
				if( mHPX[hpxid].rec != EMPTY ) {
					if( IsPointInStrip(qTheta1[nQuery],qTheta2[nQuery],mHPX[hpxid].pt) ){
						foundHPX.push_back(mHPX[hpxid]);
					}
				} 
			}

			AnalyzeResults(foundHPX,foundMRH,Matches,HPXUnique,MRHUnique,true);
			totalMatches += Matches;
			totalMrhUnique += MRHUnique;
			totalHpxUnique += HPXUnique;
			totalMrhFound += foundMRH.size();
			totalHpxFound += foundHPX.size();
		}

     	nMRHnodes = mMRH.NumNodes();
		nHPXnodes = mHPX.Npix();

		// Write trial result to file
  		fprintf(fp,"%d,%d,%d,%d,%d,%d,%d\n",
			numPoints,numTrials,nMRHnodes,nHPXnodes,Matches,MRHUnique,HPXUnique);
		if(POLY_ACC_DEBUG){
			cout << "	Wrote Trial Results to " << trialsOut.c_str() << "!\n";
		}
	}
	// End of trials, close file
	fclose(fp); 

	// Finished!
	cout << "Trials Complete!\n\n";

}



void TestTwoPointCorr(int maxDepth,int numPoints,int numQueries,double radius)
{
	int LO_IDX = 0; 
	int HI_IDX = numPoints-1; 
	int i,j;
	int64 order;
	vector<pointing> points;
	std::vector<Measurement> measurements;
	std::vector<std::pair<Measurement,Measurement>> found;
	
	// Instantiate MRH_Map data structure
	MultiResHpx_Map<Measurement> mMRH(maxDepth,NEST);
	
	// Generate random GIS longitude, latitude pairs
	GenRandomPhiTheta(numPoints,points,order,true);

	// Generate random data
	GenRandomData(numPoints,measurements);

	// Insert the GIS longitude, latitude pairs into MRH & HPX data structures
	for( i = 0; i < numPoints; i++)
	{
		// Insert next GIS longitude,latitude, data index tuple into MRH
		measurements[i].pt = points[i];  // Set spatial location of measurment
		mMRH.AddRecord(measurements[i],points[i]); 
	}

	// Compute Two-Point Correlation Query
	for( i = 0; i < numQueries; i++ ) {
		cout << "Query radius," << radius << endl << endl;

		cout << "Now do Two-Point Corrolation Query!" << endl << endl;
		found = mMRH.TwoPointCorrBin(radius);
	
		// Output results
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
}








//
//
//
//
//void RandPointInsertSaveMRHToFile(int maxDepth,int numPoints,std::string MRHOut)
//{
//	int maxDepthMRH,minDepthMRH;
//	uint64 start_s,stop_s,timeInsertMRH;
//	double avgDepthMRH;
//	int i,nMRHnodes = 0;
//	int64 sizeOfMRHNode;
//	pointing p1,p2;
//	vec3 v1,v2;
//	vector<pointing> points;
//	int64 order;
//
//	MultiResHpx mMRH = MultiResHpx(maxDepth,NEST);
//		
//	// Generate random GIS longitude, latitude pairs
//    GenRandomPhiTheta(numPoints,points,order,false);
//
//	// Insert the GIS longitude, latitude pairs into MRH data structure
//	cout << "Begin insertion of points into MRH data structure..." << endl;
//	
//	start_s=GetTimeMs64();
//	for( i = 0; i < numPoints; i++)
//	{
//		// Insert next GIS longitude,latitude, data index tuple into MRH
//		mMRH.Insert(points[i],i);
//		if(TA1_DEBUG) {
//			if( i%1000 == 0 ) {
//				cout << "Inserting Point " << i+1 << " of " << numPoints << endl;
//			}
//		}
//	}
//	stop_s=GetTimeMs64();
//	timeInsertMRH = stop_s-start_s;
//	cout << "End insertion of points into MRH data structure in " << timeInsertMRH
//			 << " ms." << endl;
//	// Compute statistics of the final data structures 
//
//	// Compute memory footprint of each data structure
//	nMRHnodes = mMRH.NumNodes();
//	sizeOfMRHNode = sizeof(MortonNode);
//
//	// Get Min,Max,Average Depth of QuadTrees of each data stucture
//	maxDepthMRH = mMRH.MaxDepth();
//	minDepthMRH = mMRH.MinDepth();
//	avgDepthMRH = mMRH.AvgDepth();
//
//	cout << "MRH stats: numPoints = " << numPoints
//		           << " numNodes = " << nMRHnodes
//		           << " sizeNode = " << sizeOfMRHNode
//		           << " minDepth = " << minDepthMRH 
//				   << " maxDepth = " << maxDepthMRH 
//				   << " avgDepth = " << avgDepthMRH << endl;
//
//	// Now output the resultant MRH data structure to file
//	mMRH.SaveToFile(MRHOut);
//
//}
//
//void LoadMRHFromFile(std::string MRHIn)
//{
//	MultiResHpx mMRH = MultiResHpx(29,NEST);	
//
//	cout << "Load MRH data structure from file: " << MRHIn.c_str() << endl;
//	if( !mMRH.LoadFromFile(MRHIn) ) {
//		cout << "MRH data structure load failed!\n";
//		exit(1);
//	}
//	//Print all trees
//	for( unsigned int i = 0; i < 12; i++ ) 
//	{
//		cout << "\n\n######################\n";
//		cout << "#### BASE CELL " << i << " ####\n";
//		cout << "######################" << endl;
//		mMRH.PrintTreeAtIndex(i); 
//	}
//
//}

struct Point {
  double x, y;
  Point(double _x, double _y) : x(_x), y(_y) {}
};


// Is p0 inside p?  Polygon 
bool my_inside(const pointing p0, const std::vector<pointing>& p) {
  size_t n = p.size();
  bool result = false;
  for (size_t i = 0; i < n; ++i) {
    size_t j = (i + 1) % n;
    if (
        // Does p0.y lies in half open y range of edge.
        // N.B., horizontal edges never contribute
        ( (p[j].phi <= p0.phi && p0.phi < p[i].phi) || 
          (p[i].phi <= p0.phi && p0.phi < p[j].phi) ) &&
        // is p to the left of edge?
        ( p0.phi < p[j].phi + (p[i].phi - p[j].phi) * (p0.theta - p[j].theta) /
          (p[i].theta - p[j].theta) )
        )
      result = !result;
  }
  return result;
}


// Improved point in polygon test which includes edge
// and vertex points

bool _equals(double a,double b) {
	double EPS = 0.00000000001;
	if( fabs(a-b) < EPS ) {
		return true;
	}
	return false;
}

bool _pointing_equals(pointing a, pointing b) {
	double EPS = 0.0000000000001;
	if( fabs(a.phi-b.phi) < EPS &&
		fabs(a.theta-b.theta) < EPS ) {
			return true;
	}
	return false;
}

bool point_in_poly(pointing pt,std::vector<pointing> poly)
{
	int i,n;
	bool inside;
	pointing p1,p2;
	double p1x,p1y,p2x,p2y,xints;
	// check if point is a vertex
	for( i = 0; i < poly.size(); i++)
	{
		if( _pointing_equals(pt,poly[i]) ) {
			return true;
		}
	}

	// check if point is on a boundary
	for( i=0; i < poly.size(); i++) {
		if( i == 0 ) {
			p1 = poly[0];
			p2 = poly[1];
		}
		else {
			p1 = poly[i-1];
			p2 = poly[i];
		}

		if( _equals(p1.phi,p2.phi) &&
			_equals(p1.phi,pt.phi) &&
			pt.theta > min(p1.theta,p2.theta) &&
			pt.theta < max(p1.theta,p2.theta) ) {
				return true;
		}
	}

	n = poly.size();
	inside = false;

	p1x = poly[0].theta;
	p1y = poly[0].phi;

	for(i = 0; i < n+1; i++) {
		p2x = poly[i%n].theta;
		p2y = poly[i%n].phi;
		if( pt.phi > min(p1y,p2y) ) {
			if( pt.phi <= max(p1y,p2y) ) {
				if( pt.theta <= max(p1x,p2x) ) {
					if( !_equals(p1y,p2y) ) {
						xints = (pt.phi-p1y)*(p2x-p1x)/(p2y-p1y)+p1x;
					}
					if( _equals(p1x,p2x) || pt.theta <= xints ) {
						inside = !inside;
					}
				}
			}
		}
		p1x = p2x;
		p1y = p2y;
	}
	if( inside == true ) {
		return true;
	}
	return false;
}


void TestGeographicLibPointInDisc
(
 double test_long, 
 double test_lat,
 double cent_long,
 double cent_lat,
 double radius
)
{
	double s12,azi1,azi2,m12,M12,M21,S12,a12;
	GeographicLib::Geodesic g(1.0,0.0);
	pointing pt,center;

	pt.theta = test_lat;
	pt.phi = test_long;
	center.theta = cent_lat;
	center.phi = cent_long;

	// Convert center & query point coordinates to GIS (degrees)
	pt = HPXtoGIS(pt);
	center = HPXtoGIS(center);
	pt.phi *= rad2degr;
	pt.theta *= rad2degr;
	center.phi *= rad2degr;
	center.theta *= rad2degr;

	// Compute geodesic given longitude,latitude (degrees) end points.
	a12 = g.Inverse(center.theta,center.phi,pt.theta,pt.phi,s12,azi1,azi2,m12,M12,M21,S12);
	
	// Convert arc length to radians
	a12 /= rad2degr;

	// Compare arc length distance between center point and query point to query radius
	if( a12 <= radius ) {
		cout << "Dist: " << a12 << " Radius: " << radius << " Point-In_Disc = TRUE!" << endl;;
	} else {
		cout << "Dist: " << a12 << " Radius: " << radius << " Point-In_Disc = FALSE!" << endl;;
	}
}

void TestGeographicLibPointInPoly(double latitude, double longitude,std::string points)
{
	// Instance GeographicLib
	GeographicLib::Gnomonic g(GeographicLib::Geodesic(1.0,0.0));
	ifstream ifp;
	double x,y,azi,rk,_lat,_long,min_lat,max_lat,min_long,max_long;
	std::vector<pointing> poly,p;
	int num_points,i;

	min_long = 99999.0;
	max_long = -99999.0;
	min_lat = 99999.0;
	max_lat = -99999.0;

	//const double lat0 = 48 + 50/60.0, lon0 = 2 + 20/60.0; // Paris
	//const double lat0 = -90.0, lon0 = 0.0;
	const double lat0 = latitude;
	const double lon0 = longitude;

	//pointing p0(latitude,longitude);
	//pointing p0(0.0,120.0);
	pointing p0proj;
	g.Forward(lat0,lon0,latitude,longitude,x,y,azi,rk);
    p0proj.theta = x;
	p0proj.phi = y;

	// Push back test polygon (degrees) (latitude,longitude)
	ifp.open(points.c_str());
	if(ifp.fail()){
	  cout << "Unable to open " << points.c_str() << "!" << endl;
	  exit(1);     
	}
	// First line is number of records
	ifp >> num_points;  
	cout << "Opened " << points.c_str() << " found " << num_points << " poly points!\n";
	skip_line(ifp); //skip rest of line
	for( i = 0; i < num_points; i++)
	{
		ifp >> _lat >> _long;
		// Skip the rest of the line
		skip_line(ifp);
		cout << _lat << " " << _long << endl;

		if(_lat < min_lat ){ min_lat = _lat; }
		if(_lat > max_lat ){ max_lat = _lat; }
		if(_long < min_long ){ min_long = _long; }
		if(_long > max_long ){ max_long = _long; }

		poly.push_back(pointing(_lat,_long));
	}
	ifp.close();

	//Set "outside" reference to be just beyond polygon bounding box
	pointing p0(max_lat+0.01,max_long+0.01);

	// Project polygon
	cout << "Theta,Phi,Easting,Northing,Azimuth,Reciprocal" << endl;
	for( int i = 0; i < poly.size(); i++ ) {
		g.Forward(lat0,lon0,poly[i].theta,poly[i].phi,x,y,azi,rk);
		p.push_back(pointing(y,x));
		cout << poly[i].theta << "," << poly[i].phi << "," << y << "," << x << "," << azi << "," << rk << endl;
	}

	// Determine if longitude,latitude point is inside/outside the polygon
	std::cout << "\nMy_Inside? " << my_inside(p0,p) << endl;
/*winner==>*/std::cout << "\nIsPointInPoly2? " << IsPointInPoly2(p0proj,p) << endl;

	std::cout << "\nPoint_In_Poly? " << point_in_poly(pointing(latitude,longitude),poly);
}


int main(int64 argc, char* argv[])

{
 	srand(time(NULL));

	int testNum, numPoints, numQueries, maxDepth, numTrials;
	std::string outFile,inFile;
	float minDist;

	testNum = atoi(argv[1]);

	cout << "Running Test Number: " << testNum << endl << endl;

	if(testNum == -1) {
		maxDepth = atoi(argv[2]);
		numPoints = atoi(argv[3]);
		TestMortonSort(maxDepth,numPoints);
	}

	// Random Point Insertion (MRH Data Structure ONLY)
	if( testNum == 0 ) {
		maxDepth = atoi(argv[2]);
		numPoints = atoi(argv[3]);
		minDist = atof(argv[4]);
		MRHRandomPointInsertion(maxDepth,numPoints,minDist);
	}

	// Random Point Insertion Bench Mark
	if( testNum == 1 ) {
		maxDepth = atoi(argv[2]);
		numPoints = atoi(argv[3]);
		minDist = atof(argv[4]);
		numTrials = atoi(argv[5]);
		outFile = argv[6];
		RandomPointInsertionBenchMark(maxDepth,numPoints,minDist,numTrials,outFile);
	}

	// Random Disc Query Bench Mark
	if( testNum == 2 ) {
		maxDepth = atoi(argv[2]);
		numPoints = atoi(argv[3]);
		minDist = atof(argv[4]);
		numQueries = atoi(argv[5]);
		numTrials = atoi(argv[6]);
		outFile = argv[7];
		RandomDiscQueryBenchMark(maxDepth,numPoints,minDist,numQueries,numTrials,outFile);
	}

	// Random Disc Query Accuracy Bench Mark
	if( testNum == 3 ) {
		maxDepth = atoi(argv[2]);
		numPoints = atoi(argv[3]);
		minDist = atof(argv[4]);
		numQueries = atoi(argv[5]);
		numTrials = atoi(argv[6]);
		outFile = argv[7];
		RandomDiscQueryAccuracyBenchMark(maxDepth,numPoints,minDist,numQueries,numTrials,outFile);
	}

	// Random Polygon Query Bench Mark
	if( testNum == 4 ) {
		maxDepth = atoi(argv[2]);
		numPoints = atoi(argv[3]);
		minDist = atof(argv[4]);
		numQueries = atoi(argv[5]);
		numTrials = atoi(argv[6]);
		outFile = argv[7];
		RandomPolyQueryBenchMark(maxDepth,numPoints,minDist,numQueries,numTrials,outFile);
	}

	// Random Polygon Query Accuracy Bench Mark
	if( testNum == 5 ) {
		maxDepth = atoi(argv[2]);
		numPoints = atoi(argv[3]);
		minDist = atof(argv[4]);
		numQueries = atoi(argv[5]);
		numTrials = atoi(argv[6]);
		outFile = argv[7];
		RandomPolyQueryAccuracyBenchMark(maxDepth,numPoints,minDist,numQueries,numTrials,outFile);
	}

	// Random Neighbor Query Bench Mark
	if( testNum == 6 ) {
		maxDepth = atoi(argv[2]);
		numPoints = atoi(argv[3]);
		numQueries = atoi(argv[4]);
		numTrials = atoi(argv[5]);
		outFile = argv[6];
		RandomNeighborsQueryAccuracyBenchMark(maxDepth,numPoints,numQueries,numTrials,outFile);
	}


	// Random Latitude Strip Query Bench Mark
	if( testNum == 7 ) {
		maxDepth = atoi(argv[2]);
		numPoints = atoi(argv[3]);
		numQueries = atoi(argv[4]);
		numTrials = atoi(argv[5]);
		outFile = argv[6];
		RandomStripQueryAccuracyBenchMark(maxDepth,numPoints,numQueries,numTrials,outFile);
	}

	// Random Two-Point Correlation Query
	if( testNum == 8 ) {
		maxDepth = atoi(argv[2]);
		numPoints = atoi(argv[3]);
		numQueries = atoi(argv[4]);
		double radius = atof(argv[5]);
		TestTwoPointCorr(maxDepth,numPoints,numQueries,radius);
	}

	if( testNum == 9 ) {
		double latitude = atof(argv[2]);
		double longitude = atof(argv[3]);
		inFile = argv[4];
		TestGeographicLibPointInPoly(latitude,longitude,inFile);

	}

	if( testNum == 10 ) {
		double test_long = atof(argv[2]);
		double test_lat = atof(argv[3]);
		double center_long = atof(argv[4]);
		double center_lat = atof(argv[5]);
		double radius = atof(argv[6]);
		TestGeographicLibPointInDisc(test_long,test_lat,center_long,center_lat,radius);

	}

	//// Random Point Insertion, Save Data Structure to File (MRH Data Structure ONLY)
	//if( testNum == 7 ) {
	//	maxDepth = atoi(argv[2]);
	//	numPoints = atoi(argv[3]);
	//	outFile = argv[4];		
	//	RandPointInsertSaveMRHToFile(maxDepth,numPoints,outFile);
	//}

	//if( testNum == 9 ) {
	//	inFile = argv[2];		
	//	LoadMRHFromFile(inFile);
	//}

	//// Test Random Polygon Generator
	//if( testNum == 10 ) {
	//	RandomConvexPoly(true);
	//}

}

