#include "RangeQueryTests.h"

using namespace std;

void DiscQueryTest(int order,float thetaDeg,float phiDeg,float radius ) {

	printf("\n######################################");
	printf("\n########## DISC QUERY TEST  ##########");
	printf("\n######################################\n");

	MultiResHpx mMRH2(NEST);
	Healpix_Map<int> hpxNEST(order,NEST);

	pointing testPt;
	vector<Point2D> queryDisc;
	rangeset<int> queryDiscHpx;
	vector<int> queryDiscList;
	Point2D nextCell;
	int rStartHpx,rEndHpx;

	testPt.theta = thetaDeg*degr2rad; 
	testPt.phi = phiDeg*degr2rad;	


	printf("Query: CoLat %5.2f Long %5.2f Radius %5.2f (z = %5.2f phi = %5.2f)\n",
		thetaDeg,
		phiDeg,
		radius,
		cos(thetaDeg*degr2rad),
		(phiDeg/180.0) );

	hpxNEST.query_disc(testPt,radius,queryDiscHpx);

	printf("\nTest #1: HEALPix\n\n");

	printf("%d Cells found in disc (range form):\n\n",queryDiscHpx.size());
	for(unsigned int i = 0; i < queryDiscHpx.size(); i++ ) {
		rStartHpx = queryDiscHpx.ivbegin(i);
		rEndHpx = queryDiscHpx.ivend(i);
		printf("From: hpxId = %d To: hpxId = %d\n",rStartHpx,rEndHpx);
	}

	queryDiscHpx.toVector(queryDiscList);
	printf("%d Cells found in disc (list form):\n\n",queryDiscHpx.size());
	for(unsigned int i = 0; i < queryDiscList.size(); i++ ) {
		printf("hpxId = %d\n",queryDiscList[i]);
	}

	queryDisc = mMRH2.query_disc(testPt,radius);

	printf("\nTest #2: MultiResHpx\n\n");
	printf("%d Cells found in disc:\n\n",queryDisc.size());
	for(unsigned int i = 0; i < queryDisc.size(); i++ ) {
		nextCell = queryDisc[i];
		printf("ix = %d iy = %d face_num = %d mapidx = %d\n",
			nextCell.ix,
			nextCell.iy,
			nextCell.fn,
			nextCell.mapidx);
	}

}

void DiscIncQueryTest(int order,float thetaDeg,float phiDeg,float radius ) {

	printf("\n#################################################");
	printf("\n##########  INCLUSIVE DISC QUERY TEST  ##########");
	printf("\n#################################################\n");

	MultiResHpx mMRH2(NEST);
	Healpix_Map<int> hpxNEST(order,NEST);

	pointing testPt;
	vector<Point2D> queryDisc;
	rangeset<int> queryDiscHpx;
	vector<int> queryDiscList;
	Point2D nextCell;
	int rStartHpx,rEndHpx;

	testPt.theta = thetaDeg*degr2rad; 
	testPt.phi = phiDeg*degr2rad;	


	printf("Query: CoLat %5.2f Long %5.2f Radius %5.2f (z = %5.2f phi = %5.2f)\n",
		thetaDeg,
		phiDeg,
		radius,
		cos(thetaDeg*degr2rad),
		(phiDeg/180.0) );

	hpxNEST.query_disc_inclusive(testPt,radius,queryDiscHpx,4);

	printf("\nTest #1: HEALPix\n\n");

	printf("%d Cells found in disc (range form):\n\n",queryDiscHpx.size());
	for(unsigned int i = 0; i < queryDiscHpx.size(); i++ ) {
		rStartHpx = queryDiscHpx.ivbegin(i);
		rEndHpx = queryDiscHpx.ivend(i);
		printf("From: hpxId = %d To: hpxId = %d\n",rStartHpx,rEndHpx);
	}

	queryDiscHpx.toVector(queryDiscList);
	printf("%d Cells found in disc (list form):\n\n",queryDiscList.size());
	for(unsigned int i = 0; i < queryDiscList.size(); i++ ) {
		printf("hpxId = %d\n",queryDiscList[i]);
	}

	queryDisc = mMRH2.query_disc_inclusive(testPt,radius,4);

	printf("\nTest #2: MultiResHpx\n\n");
	printf("%d Cells found in disc:\n\n",queryDisc.size());
	for(unsigned int i = 0; i < queryDisc.size(); i++ ) {
		nextCell = queryDisc[i];
		printf("ix = %d iy = %d face_num = %d mapidx = %d\n",
			nextCell.ix,
			nextCell.iy,
			nextCell.iy,
			nextCell.mapidx);
	}
}

void TriangleQueryTest(int order,float thetaDeg,float phiDeg) {

	printf("\n#################################################");
	printf("\n##########  TRIANGLE REGION QUERY TEST  ##########");
	printf("\n#################################################\n");

	MultiResHpx mMRH(NEST);
	Healpix_Map<int> hpxNEST(order,NEST);

	vector<pointing> testPoly;
	pointing nextVertex;
	vector<Point2D> queryResults;
	rangeset<int> queryResultsHpx;
	vector<int> queryResultsList;
	int rStartHpx,rEndHpx;
	Point2D nextCell;

	printf("Query Center: CoLat %5.2f Long %5.2f (z = %5.2f phi = %5.2f)\n",
		thetaDeg, phiDeg, cos(thetaDeg*degr2rad), (phiDeg/180.0) );
	// Create list of test polygon vertices, let input theta,phi be "center"
	// of the test polygon. Test polygon for this test will be a triangle.

	// Corner #1
	nextVertex.theta = (thetaDeg-30.0)*degr2rad;
	nextVertex.phi = (phiDeg+20.0)*degr2rad;
	printf("Corner #1: z = %5.2f phi = %5.2f)\n",cos(nextVertex.theta), ((phiDeg+20.0)/180.0) );
	testPoly.push_back(nextVertex); 

	// Corner #2
	nextVertex.theta = (thetaDeg+40.0)*degr2rad;
	nextVertex.phi = (phiDeg+40.0)*degr2rad;
	printf("Corner #2: z = %5.2f phi = %5.2f)\n",cos(nextVertex.theta), ((phiDeg+40.0)/180.0) );
	testPoly.push_back(nextVertex); 

	// Corner #3
	nextVertex.theta = (thetaDeg+20.0)*degr2rad;
	nextVertex.phi = (phiDeg-60.0)*degr2rad;
	printf("Corner #3: z = %5.2f phi = %5.2f)\n",cos(nextVertex.theta), ((phiDeg-60.0)/180.0) );
	testPoly.push_back(nextVertex); 


	printf("\nTest #1: HEALPix\n\n");
	hpxNEST.query_polygon(testPoly,queryResultsHpx);
	printf("%d Cells found in polygon (range form):\n\n",queryResultsHpx.size());
	for(unsigned int i = 0; i < queryResultsHpx.size(); i++ ) {
		rStartHpx = queryResultsHpx.ivbegin(i);
		rEndHpx = queryResultsHpx.ivend(i);
		printf("From: hpxId = %d To: hpxId = %d\n",rStartHpx,rEndHpx);
	}

	queryResultsHpx.toVector(queryResultsList);
	printf("\n%d Cells found in polygon (list form):\n\n",queryResultsList.size());
	for(unsigned int i = 0; i < queryResultsList.size(); i++ ) {
		printf("hpxId = %d\n",queryResultsList[i]);
	}

	queryResults = mMRH.query_polygon(testPoly);

	printf("\nTest #2: MultiResHpx\n\n");
	printf("%d Cells found in polygon:\n\n",queryResults.size());
	for(unsigned int i = 0; i < queryResults.size(); i++ ) {
		nextCell = queryResults[i];
		printf("ix = %d iy = %d face_num = %d mapidx = %d\n",
			nextCell.ix,
			nextCell.iy,
			nextCell.fn,
			nextCell.mapidx);
	}

}

void TriangleIncQueryTest(int order,float thetaDeg,float phiDeg) {
	printf("\n###########################################################");
	printf("\n##########  TRIANGLE INCLUSIVE REGION QUERY TEST  ##########");
	printf("\n###########################################################\n");

	MultiResHpx mMRH(NEST);
	Healpix_Map<int> hpxNEST(order,NEST);

	vector<pointing> testPoly;
	pointing nextVertex;
	vector<Point2D> queryResults;
	rangeset<int> queryResultsHpx;
	vector<int> queryResultsList;
	int rStartHpx,rEndHpx;
	Point2D nextCell;

	printf("Query Center: CoLat %5.2f Long %5.2f (z = %5.2f phi = %5.2f)\n",
		thetaDeg, phiDeg, cos(thetaDeg*degr2rad), (phiDeg/180.0) );
	// Create list of test polygon vertices, let input theta,phi be "center"
	// of the test polygon. Test polygon for this test will be a triangle.

	// Corner #1
	nextVertex.theta = (thetaDeg-30.0)*degr2rad;
	nextVertex.phi = (phiDeg+20.0)*degr2rad;
	printf("Corner #1: z = %5.2f phi = %5.2f)\n",cos(nextVertex.theta), ((phiDeg+20.0)/180.0) );
	testPoly.push_back(nextVertex); 

	// Corner #2
	nextVertex.theta = (thetaDeg+40.0)*degr2rad;
	nextVertex.phi = (phiDeg+40.0)*degr2rad;
	printf("Corner #2: z = %5.2f phi = %5.2f)\n",cos(nextVertex.theta), ((phiDeg+40.0)/180.0) );
	testPoly.push_back(nextVertex); 

	// Corner #3
	nextVertex.theta = (thetaDeg+20.0)*degr2rad;
	nextVertex.phi = (phiDeg-60.0)*degr2rad;
	printf("Corner #3: z = %5.2f phi = %5.2f)\n",cos(nextVertex.theta), ((phiDeg-60.0)/180.0) );
	testPoly.push_back(nextVertex); 


	printf("\nTest #1: HEALPix\n\n");
	hpxNEST.query_polygon_inclusive(testPoly,queryResultsHpx,4);
	printf("%d Cells found in polygon (range form):\n\n",queryResultsHpx.size());
	for(unsigned int i = 0; i < queryResultsHpx.size(); i++ ) {
		rStartHpx = queryResultsHpx.ivbegin(i);
		rEndHpx = queryResultsHpx.ivend(i);
		printf("From: hpxId = %d To: hpxId = %d\n",rStartHpx,rEndHpx);
	}

	queryResultsHpx.toVector(queryResultsList);
	printf("\n%d Cells found in polygon (list form):\n\n",queryResultsList.size());
	for(unsigned int i = 0; i < queryResultsList.size(); i++ ) {
		printf("hpxId = %d\n",queryResultsList[i]);
	}

	queryResults = mMRH.query_polygon_inclusive(testPoly,4);

	printf("\nTest #2: MultiResHpx\n\n");
	printf("%d Cells found in polygon:\n\n",queryResults.size());
	for(unsigned int i = 0; i < queryResults.size(); i++ ) {
		nextCell = queryResults[i];
		printf("ix = %d iy = %d face_num = %d mapidx = %d\n",
			nextCell.ix,
			nextCell.iy,
			nextCell.fn,
			nextCell.mapidx);
	}
}

void StripQueryTest(int order,float thetaDeg,float phiDeg) {
	printf("\n#########################################################");
	printf("\n########## COLATITUDE STRIP QUERY TEST (RING)  ##########");
	printf("\n#########################################################\n");

	MultiResHpx mMRH(RING);
	Healpix_Map<int> hpxNEST(order,RING);

	pointing testPt;
	double theta1,theta2;
	vector<Point2D> queryStrip;
	rangeset<int> queryStripHpx;
	vector<int> queryStripList;
	Point2D nextCell;
	int rStartHpx,rEndHpx;

	// Query Center
	testPt.theta = thetaDeg*degr2rad; 
	testPt.phi = phiDeg*degr2rad;	
	printf("Query Center: CoLat %5.2f Long %5.2f (z = %5.2f phi = %5.2f)\n",
		thetaDeg,phiDeg,cos(thetaDeg*degr2rad),(phiDeg/180.0) );

	// Strip Query Range Theta1 -> Theta2
	theta1 = (thetaDeg-20.0);
	theta2 = (thetaDeg+20.0);
	printf("\nStrip CoLatitude Query Range: theta1 = %5.2f To theta2 = %5.2f\n",
		cos(theta1*degr2rad),cos(theta2*degr2rad));

	hpxNEST.query_strip(theta1*degr2rad,theta2*degr2rad,false,queryStripHpx);

	printf("\nTest #1: HEALPix\n\n");

	printf("%d Cells found in strip (range form):\n\n",queryStripHpx.size());
	for(unsigned int i = 0; i < queryStripHpx.size(); i++ ) {
		rStartHpx = queryStripHpx.ivbegin(i);
		rEndHpx = queryStripHpx.ivend(i);
		printf("From: hpxId = %d To: hpxId = %d\n",rStartHpx,rEndHpx);
	}

	queryStripHpx.toVector(queryStripList);
	printf("%d Cells found in strip (list form):\n\n",queryStripList.size());
	for(unsigned int i = 0; i < queryStripList.size(); i++ ) {
		printf("hpxId = %d\n",queryStripList[i]);
	}

	queryStrip = mMRH.query_strip(theta1*degr2rad,theta2*degr2rad,false);

	printf("\nTest #2: MultiResHpx\n\n");
	printf("%d Cells found in strip:\n\n",queryStrip.size());
	for(unsigned int i = 0; i < queryStrip.size(); i++ ) {
		nextCell = queryStrip[i];
		printf("ix = %d iy = %d face_num = %d mapidx = %d\n",
			nextCell.ix,
			nextCell.iy,
			nextCell.fn,
			nextCell.mapidx);
	}
}

void StripIncQueryTest(int order,float thetaDeg,float phiDeg) {
	printf("\n###################################################################");
	printf("\n########## INCLUSIVE COLATITUDE STRIP QUERY TEST (RING)  ##########");
	printf("\n###################################################################\n");

	MultiResHpx mMRH(RING);
	Healpix_Map<int> hpxNEST(order,RING);

	pointing testPt;
	double theta1,theta2;
	vector<Point2D> queryStrip;
	rangeset<int> queryStripHpx;
	vector<int> queryStripList;
	Point2D nextCell;
	int rStartHpx,rEndHpx;

	// Query Center
	testPt.theta = thetaDeg*degr2rad; 
	testPt.phi = phiDeg*degr2rad;	
	printf("Query Center: CoLat %5.2f Long %5.2f (z = %5.2f phi = %5.2f)\n",
		thetaDeg,phiDeg,cos(thetaDeg*degr2rad),(phiDeg/180.0) );

	// Strip Query Range Theta1 -> Theta2
	theta1 = (thetaDeg-20.0);
	theta2 = (thetaDeg+20.0);
	printf("\nStrip CoLatitude Query Range: theta1 = %5.2f To theta2 = %5.2f\n",
		cos(theta1*degr2rad),cos(theta2*degr2rad));

	hpxNEST.query_strip(theta1*degr2rad,theta2*degr2rad,true,queryStripHpx);

	printf("\nTest #1: HEALPix\n\n");

	printf("%d Cells found in strip (range form):\n\n",queryStripHpx.size());
	for(unsigned int i = 0; i < queryStripHpx.size(); i++ ) {
		rStartHpx = queryStripHpx.ivbegin(i);
		rEndHpx = queryStripHpx.ivend(i);
		printf("From: hpxId = %d To: hpxId = %d\n",rStartHpx,rEndHpx);
	}

	queryStripHpx.toVector(queryStripList);
	printf("%d Cells found in strip (list form):\n\n",queryStripList.size());
	for(unsigned int i = 0; i < queryStripList.size(); i++ ) {
		printf("hpxId = %d\n",queryStripList[i]);
	}

	queryStrip = mMRH.query_strip(theta1*degr2rad,theta2*degr2rad,true);

	printf("\nTest #2: MultiResHpx\n\n");
	printf("%d Cells found in strip:\n\n",queryStrip.size());
	for(unsigned int i = 0; i < queryStrip.size(); i++ ) {
		nextCell = queryStrip[i];
		printf("ix = %d iy = %d face_num = %d mapidx = %d\n",
			nextCell.ix,
			nextCell.iy,
			nextCell.fn,
			nextCell.mapidx);
	}
}

void NeighborsQueryTest(int order,float thetaDeg,float phiDeg) {
	long hpxID;
	pointing testPt;
	fix_arr<int,8> neighbors;
	vector<Point2D> neighborsMRH;
	vector<int> mrhID;

	MultiResHpx mMRH(NEST);
	Healpix_Map<int> hpxNEST(order,NEST);

	printf("\n#########################################");
	printf("\n########## NEIGHBOR QUERY TEST ##########");
	printf("\n#########################################\n");

	printf("\nQuery Inputs: Order = %d, Lat = %5.2f, Long = %5.2f\n",order,thetaDeg,phiDeg);
	printf("z: %5.2f, phi: %5.2f\n\n",cos(thetaDeg*degr2rad),(phiDeg/180.0));

	printf("\nTest #1 HEALPix\n");
	testPt.theta = thetaDeg*degr2rad; //convert Latitude to CoLatitude then to Radians
	testPt.phi = phiDeg*degr2rad;
	hpxID = hpxNEST.ang2pix(testPt);
	hpxNEST.neighbors(hpxID,neighbors);

	// Print out the neighbor HEALPix indices
	printf("SW = %d\n",neighbors[0]);
	printf("W = %d\n",neighbors[1]);
	printf("NW = %d\n",neighbors[2]);
	printf("N = %d\n",neighbors[3]);
	printf("NE = %d\n",neighbors[4]);
	printf("E = %d\n",neighbors[5]);
	printf("SE = %d\n",neighbors[6]);
	printf("S = %d\n",neighbors[7]);

	printf("\nTest #2 MultiResHpx\n");
	mrhID = mMRH.ang2mrh(testPt);
	neighborsMRH = mMRH.neighbors(mrhID);

	// Print out the neighbor HEALPix indices
	printf("SW = ix %d iy %d f %d mapidx %d\n",neighborsMRH[0].ix,neighborsMRH[0].iy,neighborsMRH[0].fn,neighborsMRH[0].mapidx);
	printf("W = ix %d iy %d f %d mapidx %d\n",neighborsMRH[1].ix,neighborsMRH[1].iy,neighborsMRH[1].fn,neighborsMRH[1].mapidx);
	printf("NW = ix %d iy %d f %d mapidx %d\n",neighborsMRH[1].ix,neighborsMRH[2].iy,neighborsMRH[2].fn,neighborsMRH[2].mapidx);
	printf("N = ix %d iy %d f %d mapidx %d\n",neighborsMRH[3].ix,neighborsMRH[3].iy,neighborsMRH[3].fn,neighborsMRH[3].mapidx);
	printf("NE = ix %d iy %d f %d mapidx %d\n",neighborsMRH[4].ix,neighborsMRH[4].iy,neighborsMRH[4].fn,neighborsMRH[4].mapidx);
	printf("E = ix %d iy %d f %d mapidx %d\n",neighborsMRH[5].ix,neighborsMRH[5].iy,neighborsMRH[5].fn,neighborsMRH[5].mapidx);
	printf("SE = ix %d iy %d f %d mapidx %d\n",neighborsMRH[6].ix,neighborsMRH[6].iy,neighborsMRH[6].fn,neighborsMRH[6].mapidx);
	printf("S = ix %d iy %d f %d mapidx %d\n",neighborsMRH[7].ix,neighborsMRH[7].iy,neighborsMRH[7].fn,neighborsMRH[7].mapidx);

}