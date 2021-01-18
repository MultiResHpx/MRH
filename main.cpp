#include <stdlib.h>
#include <string>
//TESTING SUITE INCLUDES
#include "TestMortonLQT.h"
#include "TestMultiResHpx.h"

int main(int64 argc, char* argv[])

{
	int testNum, order, level,max_depth,num_insert,data;
	double longitude,latitude;
	Morton morton;
	int64 hpxIdx;
	std::string filename;

	testNum = atoi(argv[1]);

	if(MORTONLQTDEBUG) {
		cout << "Test #" << testNum << endl;
	}

	if(testNum == -2) {
		cout << "Set HPX Map to different orders and sizes:" << endl;
		TestBuildHpx();
	}

	if(testNum == -1) {
		morton = StringToMorton(argv[2]);
		cout << "Morton: "; PrintMorton(morton); cout << endl;
		cout << "Morton Str: " << MortonToString(morton) << endl;

	}

	if(testNum == 0) {
		//order = atoi(argv[2]);
		//int face_num = atoi(argv[3]);
		//int step = atoi(argv[4]);
  //      OutputBaseCellBoundary(order,face_num,step);
		TestMask();
	}

//####
//#### MORTONLQT TESTS
//####

    // HEALPix to Morton
	if(testNum == 1 ) {
		hpxIdx = StringToHpx(argv[2]);
		order = atoi(argv[3]);
      TestHpxToMorton(hpxIdx,order);
	}
	
	// Morton to HEALPix
	if(testNum == 2 ) {
		morton = StringToMorton(argv[2]);
        TestMortonToHpx(morton);
	}

	// All NSIDE HEALPix to Morton
	if(testNum == 3 ) {
		order = atoi(argv[2]);
        TestAllNSIDE8HpxToMorton(order);
	}

	// Longitude, Latitude to Morton
	if(testNum == 4 ) {
	   longitude = atof(argv[2]);
	   latitude = atof(argv[3]);
	   level = atoi(argv[4]);
       TestPhiThetaToMorton(longitude,latitude,level);
	}

	// Morton to Longitude, Latitude
	if(testNum == 5 ) {
	   morton = StringToMorton(argv[2]);
       TestMortonToPhiTheta(morton);
	}

	// Parent of Morton
	if(testNum == 6) {
	   morton = StringToMorton(argv[2]);
       TestParentOfMorton(morton);
	}

	// Children of Morton
	if(testNum == 7) {
	   morton = StringToMorton(argv[2]);
       TestChildrenOfMorton(morton);
	}

	// Siblings of Morton
	if(testNum == 8) {
	   morton = StringToMorton(argv[2]);
       TestSiblingsOfMorton(morton);
	}

	// Insert Morton Node
	if(testNum == 9) {
	   max_depth = atoi(argv[2]);
	   num_insert = atoi(argv[3]);
      TestInsertNode(max_depth,num_insert,argv,1);
	}

	// Insert Morton Node DEGREES
	if(testNum == 10) {
	   max_depth = atoi(argv[2]);
	   num_insert = atoi(argv[3]);
       TestInsertNode(max_depth,num_insert,argv,0);
	}

	// Search Morton Node
	if(testNum == 11) {
	   max_depth = atoi(argv[2]);
	   num_insert = atoi(argv[3]);
       TestSearchNode(max_depth,num_insert,argv);
	}

	// Delete Morton Node
	if(testNum == 12) {
	   max_depth = atoi(argv[2]);
	   num_insert = atoi(argv[3]);
      TestDeleteNode(max_depth,num_insert,argv);
	}

	// Write Morton LQT
	if(testNum == 13) {
	   max_depth = atoi(argv[2]);
	   num_insert = atoi(argv[3]);
	   filename = argv[4];
      TestWriteTree(max_depth,num_insert,filename,argv);
	}

	// Load Morton LQT
	if(testNum == 14) {
	   max_depth = atoi(argv[2]);
	   filename = argv[3];
      TestLoadTree(max_depth,filename);
	}

//####
//#### MULTIRESHPX TESTS
//####
	// Insert Data Index into Base 0 Cell
	if(testNum == 15) {
  	   max_depth = atoi(argv[2]);
       TestMRHInsertNodeBase0(max_depth);
	}

	// Insert Data Index into any Base Cell and Search for it
	if(testNum == 16) {
  	   max_depth = atoi(argv[2]);
	   longitude = atof(argv[3]);
	   latitude = atof(argv[4]);
       data = atoi(argv[5]);
       hpxIdx = StringToHpx(argv[6]);
	   order = atoi(argv[7]);
       TestMRHInsertOneNodeSearchOneNode(max_depth,longitude,latitude,data,hpxIdx,order);
	}

	// Insert N Data Indices from any Base Cell
	if(testNum == 17) {
  	   max_depth = atoi(argv[2]);
	   int seed = atoi(argv[3]);
       int num_points = atoi(argv[4]);
	   TestMRHInsertNodeAllBase(max_depth,seed,num_points);
	}

	// Insert N Data Indices from any Base Cell then Delete them all, one by one.
	if(testNum == 18) {
  	   max_depth = atoi(argv[2]);
	   int seed = atoi(argv[3]);
       int num_points = atoi(argv[4]);
	   TestMRHInsertNodeAllBaseThenDeleteAll(max_depth,seed,num_points);
	}

	// Insert N Data Indices from any Base Cell then do Disc Query and output results
	if(testNum == 19) {
  	   max_depth = atoi(argv[2]);
	   int seed = atoi(argv[3]);
       int num_points = atoi(argv[4]);
	   double longitude = atof(argv[5]);
	   double latitude = atof(argv[6]);
	   double radius = atof(argv[7]);
  	   TestMRHQueryDiscRandom(max_depth,seed,num_points,longitude,latitude,radius);
	}

	// Insert N Data Indices from any Base Cell then do Disc Query Inclusive and output results
	if(testNum == 20) {
  	   max_depth = atoi(argv[2]);
	   int seed = atoi(argv[3]);
       int num_points = atoi(argv[4]);
	   double longitude = atof(argv[5]);
	   double latitude = atof(argv[6]);
	   double radius = atof(argv[7]);
	   int fact = atoi(argv[8]);
  	   TestMRHQueryDiscInclusiveRandom(max_depth,seed,num_points,longitude,latitude,radius,fact);
	}


	// Insert One Data Index then do Disc Query and output results
	if(testNum == 21) {
  	   max_depth = atoi(argv[2]);
	   double longitude = atof(argv[3]);
	   double latitude = atof(argv[4]);
	   int data = atoi(argv[5]);
	   double qlong = atof(argv[6]);
	   double qlat = atof(argv[7]);
	   double radius = atof(argv[8]);
  	   TestMRHInsertOneNodeQueryDisc(max_depth,longitude,latitude,data,qlong,qlat,radius);
	}

	// Insert One Data Index then do Disc Query Inclusive and output results
	if(testNum == 22) {
  	   max_depth = atoi(argv[2]);
	   double longitude = atof(argv[3]);
	   double latitude = atof(argv[4]);
	   int data = atoi(argv[5]);
	   double qlong = atof(argv[6]);
	   double qlat = atof(argv[7]);
	   double radius = atof(argv[8]);
	   int fact = atoi(argv[9]);
  	   TestMRHInsertOneNodeQueryDiscInclusive(max_depth,longitude,latitude,data,qlong,qlat,radius,fact);
	}

	// Insert One Data Index then do Polygon Query and output results
	if(testNum == 23) {
  	   max_depth = atoi(argv[2]);
	   double longitude = atof(argv[3]);
	   double latitude = atof(argv[4]);
	   int data = atoi(argv[5]);
	   int numpoints = atoi(argv[6]);
  	   TestMRHInsertOneNodeQueryPolygon(max_depth,longitude,latitude,data,numpoints,argv);
	}


	// Insert One Data Index then do Polygon Inclusive Query and output results
	if(testNum == 24) {
  	   max_depth = atoi(argv[2]);
	   double longitude = atof(argv[3]);
	   double latitude = atof(argv[4]);
	   int data = atoi(argv[5]);
	   int numpoints = atoi(argv[6]);
	   int fact = atoi(argv[7]);
  	   TestMRHInsertOneNodeQueryPolygonInclusive(max_depth,longitude,latitude,data,numpoints,fact,argv);
	}

	// Insert One Data Index then do Strip Query and output results
	if(testNum == 25) {
	   bool incYN;
  	   max_depth = atoi(argv[2]);
	   double longitude = atof(argv[3]);
	   double latitude = atof(argv[4]);
	   int data = atoi(argv[5]);
	   double lat1 = atof(argv[6]);
	   double lat2 = atof(argv[7]);
	   int inclusive = atoi(argv[8]);
	   if( inclusive == 0)
	     incYN = false;
	   else
	     incYN = true;
  	   TestMRHInsertOneNodeQueryStrip(max_depth,longitude,latitude,data,lat1,lat2,incYN);
	}


	// Insert One Data Index then do Neighbors Query and output results
	if(testNum == 26) {
  	   max_depth = atoi(argv[2]);
	   double longitude = atof(argv[3]);
	   double latitude = atof(argv[4]);
	   int data = atoi(argv[5]);
	   double qLong = atof(argv[6]);
	   double qLat = atof(argv[7]);
  	   TestMRHNeighbors(max_depth,longitude,latitude,data,qLong,qLat);
	}

	// Test of various HEALPix coordinate transformers and informers
	if(testNum == 27) {
  	   max_depth = atoi(argv[2]);
	   double longitude = atof(argv[3]);
	   double latitude = atof(argv[4]);
	   int order = atoi(argv[5]);
	   int step = atoi(argv[6]);
       TestMRHTransformersAndInfo(max_depth,longitude,latitude,order,step);
	}

	// Test HPX index to Base 0 HPX Normalization
	if(testNum == 28 ) {
  	   order = atoi(argv[2]);
       int numPoints = atoi(argv[3]);
	   TestHPXIndexToBase0Normalization(order,numPoints);
	}



}

