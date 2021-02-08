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

	// Create some random data with random point locations to input into MultiResHpx_Map Objects

	// Create some random queries 



}

