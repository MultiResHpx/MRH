#include "MultiResHpx.h"
// STD includes
#include <time.h>
#include <vector>
#include <math.h>
#include <cstdlib>
#include "lsconstants.h"
#include "healpix_map.h"


#define MRHDEBUG true
#define EMPTY -99999


class HPXNode 
{
public:
	HPXNode(){longitude = -99999.0; latitude = -99999.0; data = EMPTY;};
	HPXNode(double _longitude, double _latitude, int _data)
	{
		longitude = _longitude;
		latitude = _latitude;
		data = _data;
	}
	~HPXNode(){};

	double longitude;
	double latitude;
	int data;
};


void TestBuildHpx();
//void TestBitPackingUnpacking(int order,double phi, double theta);
//void TestNewHpx2MortonMapping(int order,double phi1, double z1, double phi2, double z2);
void TestMRHInsertNodeBase0(int maxdepth);
void TestMRHInsertOneNodeSearchOneNode(int maxdepth,double longitude, double latitude, int data, int64 hpxid, int order);
void TestMRHInsertNodeAllBase(int maxdepth,int seed,int num_points);
void TestMRHDeleteNode(int maxdepth,int num_insert, char* argv[]);
void TestMRHInsertNodeAllBaseThenDeleteAll(int maxdepth,int seed,int num_points);

void TestMRHTransformersAndInfo(int maxdepth,double longitude,double latitude,int order,int step);
void TestHPXIndexToBase0Normalization(int maxDepth,int numPoints);

void TestMRHQueryDiscRandom(int maxdepth,int seed,int num_points,double longitude, double latitude, double radius);
void TestMRHQueryDiscInclusiveRandom(int maxdepth,int seed,int num_points,double longitude, double latitude, double radius, int fact);
void TestMRHInsertOneNodeQueryDisc(int max_depth,double longitude,double latitude,int data,double qLong,double qLat,double qRad);
void TestMRHInsertOneNodeQueryDiscInclusive(int max_depth,double longitude,double latitude,int data,double qLong,double qLat,double qRad,int fact);
void TestMRHInsertOneNodeQueryPolygon(int max_depth,double longitude,double latitude,int data,int num_points,char* argv[]);
void TestMRHInsertOneNodeQueryPolygonInclusive(int max_depth,double longitude,double latitude,int data,int num_points,int fact,char* argv[]);
void TestMRHInsertOneNodeQueryStrip(int max_depth,double longitude,double latitude,int data,double lat1,double lat2,bool inclusive);
void TestMRHNeighbors(int max_depth,double longitude,double latitude,int data,double qLong,double qLat);

//void OutputBaseCellBoundary(int order,int fn,int step);
//void PointInWhichQuad(int hpx1,int hpx2, int hpx3, int hpx4, int order, double pPhi, double pTheta);
//void TestHpxToMorton(int hpxID, int order);
//void TestMortonToHpx(int morton);
//void TestAllNSIDE8HpxToMorton(int order);
//void TestPhiThetaToMorton(double longitude, double latitude, int level);
//void TestMortonToPhiTheta(int morton);
//void TestParentOfMorton(int morton);
//void TestChildrenOfMorton(int morton);
//void TestSiblingsOfMorton(int morton);
//void TestInsertNode(int maxdepth,int num_insert, char* argv[]);
//void TestSearchNode(int maxdepth,int num_insert, char* argv[]);
//void TestWriteTree(int maxdepth,int num_insert, std::string filename,char* argv[]);
//void TestLoadTree(int maxdepth,std::string filename);

//Utility Methods
double ArcDist(double th1,double phi1,double th2,double phi2);
double GSDist(double th1,double phi1,double th2,double phi2);
double CordDist(double th1,double phi1,double th2,double phi2);
