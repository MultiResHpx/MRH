#include "PointQueryTests.h"


#include <iterator>


void PointQueryTest(int hpxOrder,float thetaDeg,float phiDeg) {

//### TEST POINT QUERY ###//
	long hpxID;
	pointing testPt;
	srand(time(NULL));

	printf("\n######################################");
	printf("\n########## POINT QUERY TEST ##########");
	printf("\n######################################\n");

	printf("\nQuery Inputs: Order = %d, Lat = %5.2f, Long = %5.2f\n",hpxOrder,thetaDeg,phiDeg);
	printf("z: %5.2f, phi: %5.2f\n\n",cos(thetaDeg*degr2rad),(phiDeg/180.0));

	printf("\n#### Baseline HEALPix Test  ####\n");
	
	//NEST TEST
	Healpix_Map<int> hpxNEST(hpxOrder,NEST);
	testPt.theta = thetaDeg*degr2rad; //convert Latitude to CoLatitude then to Radians
	testPt.phi = phiDeg*degr2rad;
	hpxID = hpxNEST.ang2pix(testPt);
	printf("\nNEST hpxID = %d\n",hpxID);
	
	//RING TEST
	Healpix_Map<int> hpxRING(hpxOrder,RING);
	testPt.theta = thetaDeg*degr2rad; //convert Latitude to CoLatitude then to Radians
	testPt.phi = phiDeg*degr2rad;
	hpxID = hpxRING.ang2pix(testPt);
	printf("\nRING hpxID = %d\n",hpxID);


	printf("\n#### Multi-Resolution HEALPix Test ####\n");
	printf("\n RING & NEST indices should be the same in MRH since they access");
	printf("\n the same physical base QuadTree and QuadTree leaf.\n");
	
	//NEST TEST
	std::vector<int> id(3);
	int randValue1;
	int randValue2;
	int nodeCount = 0;
	MultiResHpx mrhNEST(NEST); 
	std::vector<Point2D> pointQuery;
	testPt.theta = thetaDeg*degr2rad;
	testPt.phi = phiDeg*degr2rad;
	id = mrhNEST.ang2mrh(testPt);
	printf("\nNEST mpxID: ix = %d iy = %d face_num = %d\n",id[0],id[1],id[2]);
	randValue1 = rand();
	printf("Random value to set at index = %d\n",randValue1);
	mrhNEST.Insert(testPt,randValue1);
	pointQuery = mrhNEST.At(testPt);
	for( unsigned int i = 0; i < pointQuery.size(); i++ ) 
	{
		printf("Value at index = %d\n",pointQuery[i].mapidx);
	}	

	//RING TEST
	MultiResHpx mrhRING(RING); 
	testPt.theta = thetaDeg*degr2rad;
	testPt.phi = phiDeg*degr2rad;
	id = mrhRING.ang2mrh(testPt);
	printf("\nRING mpxID: ix = %d iy = %d face_num = %d\n",id[0],id[1],id[2]);
	randValue2 = rand();
	printf("Random value to set at index = %d\n",randValue2);
	mrhRING.Insert(testPt,randValue2);
	pointQuery = mrhRING.At(testPt);
	for( unsigned int i = 0; i < pointQuery.size(); i++ ) 
	{
		printf("Value at index = %d\n",pointQuery[i].mapidx);
	}	

	//Now test the query_node_with_index method by setting known MRH index to value then
	//calling the query to search the MRH data structure for MRH index that is holding that value
	nodeCount = mrhNEST.query_nodecount_with_index(randValue1);
	
	printf("Found %d in QuadTree at %d locations:\n",randValue1,nodeCount);
	//for( unsigned int i = 0; i < pointQuery.size(); i++ ) 
	//{	
 // 	   printf("\nFound %d at NEST mpxID: ix = %d iy = %d face_num = %d\n",pointQuery[i].mapidx,
	//																	  pointQuery[i].ix,
	//																	  pointQuery[i].iy,
	//																	  pointQuery[i].fn);
	//}

   //Test failure
	pointQuery = mrhNEST.query_nodes_with_index(randValue1+1);
	for( unsigned int i = 0; i < pointQuery.size(); i++ ) 
	{	
 	   printf("\nSearch for %d: ix = %d iy = %d face_num = %d\n",pointQuery[i].mapidx,
																		  pointQuery[i].ix,
																		  pointQuery[i].iy,
																		  pointQuery[i].fn);
	}


}



void DataPointQueryTest(int hpxOrder) {
	pointing testPt;
	std::vector<int> id(3);
	std::vector<Point2D> pointQuery;

	int idT1,idT2;
	int hpxID,Nside,NsideNew;


	printf("\n###########################################");
	printf("\n########## DATA POINT QUERY TEST ##########");
	printf("\n###########################################\n\n");


	printf("Test Order: %d, NSide: %d\n\n",hpxOrder,int(1)<<hpxOrder);             

	// Create dummy data at variety of latitudes, longitudes
	MultiResHpx_Map<Data> mMRH(NEST);
	Healpix_Map<int> hpxNEST(hpxOrder,NEST); //Use HEALPix map for comparision
	Healpix_Map<int> hpxTEST(hpxOrder,NEST); //Use HEALPix map for comparision



	for( int i = 0; i < 10; i++)
	{
		// Store the incoming data in std::vector<Measurement>
		Data next;
		next.latitude = 180.0*(float(rand())/float(RAND_MAX));
		next.longitude = 360*(float(rand())/float(RAND_MAX));
		next.meas1 = float (rand())/float(RAND_MAX);
		next.meas2 = float (rand())/float(RAND_MAX);
		next.meas3 = float (rand())/float(RAND_MAX);
		next.meas4 = float (rand())/float(RAND_MAX);

		// Some utility parameters to store as well, used to compare 
		// HEALPix ID to z vs longOverPI HEALPix plot of Goski.
		next.z = cos(next.latitude*degr2rad);
		next.phi = (next.longitude/180.0);

		// Compute the Multi-Res HEALPix indices
		testPt.theta = next.latitude*degr2rad;
		testPt.phi = next.longitude*degr2rad;

		// Compute corresponding HEALPix index to Validate
		hpxID = hpxNEST.ang2pix(testPt);
		next.hpxId = hpxID;

		// Add data record to map and bind it's location in MRH Data structure

		// Want to test the re-building algorithm in MRH so let's figure out where this data point
		// should map to in MRH initially; could change if there is a data insertion collision
		id = mMRH.ang2mrh(testPt);
		printf("\nMapId %d with theta: %5.4f phi: %5.4f\n  Initially Maps to -> ix: %d iy: %d fn: %d hpix: %d\n",i,testPt.theta/degr2rad,testPt.phi/degr2rad,id[0],id[1],id[2],hpxID);

        mMRH.AddRecord(next,testPt);
	   
		// Now let's look at where the data insertion actually ended up mapping to.
	    pointQuery = mMRH.GetPointsWithMapIdx(i);
		for( unsigned int i = 0; i < pointQuery.size(); i++ ) 
		{
	       testPt = mMRH.mrh2ang(pointQuery[i].ix,pointQuery[i].iy,pointQuery[i].fn);
	       hpxID = hpxNEST.ang2pix(testPt);
	       printf("\n  Finally Mapped to -> ix: %d iy: %d fn: %d hpix: %d with center theta: %5.4f phi: %5.4f\n",id[0],id[1],id[2],hpxID,testPt.theta/degr2rad,testPt.phi/degr2rad);
		}
	}

	// Evaluate test data to map associations
	pointing dp;
	for(unsigned int i = 0; i < mMRH.Map().size(); i++) 
	{
		pointQuery = mMRH.GetPointsWithMapIdx(i);
		for( unsigned int j = 0; j < pointQuery.size(); j++ )
		{
		   testPt = mMRH.mrh2ang(pointQuery[j].ix,pointQuery[j].iy,pointQuery[j].fn);
		   hpxID = hpxNEST.ang2pix(testPt);
		   printf("\nMapId %d, theta: %5.4f phi: %5.4f, found at --> ix: %d iy: %d fn: %d with center theta: %5.4f phi: %5.4f hpix: %d\n",
		          i,dp.theta/degr2rad,dp.phi/degr2rad,pointQuery[j].ix,pointQuery[j].iy,pointQuery[j].fn,testPt.theta/degr2rad,testPt.phi/degr2rad,hpxID);
		}
	}
	//for( unsigned int i = 0; i < 10; i++ )
	//{
	//	if( mMRH.GetLocation(i,id2) ) {
	//	   testPt = mMRH.mrh2ang(id2[0],id2[1],id2[2]);
	//	   hpxID = hpxNEST.ang2pix(testPt);
	//	   printf("\nMapId %d found at --> ix: %d iy: %d fn: %d with center lat: %5.2f long: %5.2f hpix: %d\n",
	//	          i,id2[0],id2[1],id2[2],testPt.theta/degr2rad,testPt.phi/degr2rad,hpxID);
	//	}
	//	else {
 //          printf("MapId %d NOT FOUND!\n",i);
	//	}
	//}

	// Sanity check
	//Healpix_Map<int> tHPXNEST(3,NEST); //Use HEALPix map for comparision
 //   testPt.theta = 126.23*degr2rad;
 //   testPt.phi   = 328.64*degr2rad;
 //   idT1 = tHPXNEST.ang2pix(testPt);

 //   testPt.theta = 133.71*degr2rad;
 //   testPt.phi   = 345.86*degr2rad;
 //   idT2 = tHPXNEST.ang2pix(testPt);

	//printf("idT1: %d idT2: %d\n",idT1,idT2);

	// Now dump out the contents of the Multi-Res HEALPix structure
	//for(unsigned int face = 0; face < 12; face++) {
	//	for(unsigned int ix = 0; ix < mMRH2.getNsideAtFace(face); ix++) {
	//		for(unsigned int iy = 0; iy < mMRH2.getNsideAtFace(face); iy++ ) {
	//			if( mMRH2.at(ix,iy,face) != -1 ) {
	//				printf("Lat. = %5.2f, Long. = %5.2f, z = %5.2f, phi = %5.2f\n",
	//						vData[mMRH2.at(ix,iy,face)].latitude,
	//						vData[mMRH2.at(ix,iy,face)].longitude,
	//						vData[mMRH2.at(ix,iy,face)].z,
	//						vData[mMRH2.at(ix,iy,face)].phi);

	//				printf("mpxID: ix = %d iy = %d face_num = %d," 
	//						 "hpxID: %d = %5.2f %5.2f %5.2f %5.2f\n\n",
	//						ix,iy,face,vData[mMRH2.at(ix,iy,face)].hpxId,
	//						vData[mMRH2.at(ix,iy,face)].meas1,
	//						vData[mMRH2.at(ix,iy,face)].meas2,
	//						vData[mMRH2.at(ix,iy,face)].meas3,
	//						vData[mMRH2.at(ix,iy,face)].meas4 );
	//			}
	//		}
	//	}
	//}
}


void InterpolationQueryTest(int order,float thetaDeg,float phiDeg) {
	pointing testDir;
	fix_arr<int,4> pixels;
	fix_arr<double,4> weights;

	MultiResHpx_Map<int> mMRH(NEST);
	Healpix_Map<int> hpxNEST(order,NEST);

	printf("\n##############################################");
	printf("\n########## INTERPOLATION QUERY TEST ##########");
	printf("\n##############################################\n");

	printf("\nQuery Inputs: Order = %d, Theta Dir= %5.2f, Phi Dir= %5.2f\n",order,thetaDeg,phiDeg);

	printf("\nTest #1 HEALPix\n");
	testDir.theta = thetaDeg*degr2rad; //convert Latitude to CoLatitude then to Radians
	testDir.phi = phiDeg*degr2rad;
	hpxNEST.get_interpol(testDir,pixels,weights);

	// Print out the neighbor HEALPix indices
	printf("Pixel #1 = %d\n",pixels[0]);
	printf("Weight #1 = %5.2f\n",weights[0]);
	printf("Pixel #2 = %d\n",pixels[1]);
	printf("Weight #2 = %5.2f\n",weights[1]);
	printf("Pixel #3 = %d\n",pixels[2]);
	printf("Weight #3 = %5.2f\n",weights[2]);
	printf("Pixel #4 = %d\n",pixels[3]);
	printf("Weight #4 = %5.2f\n",weights[3]);

}