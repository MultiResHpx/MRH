#include "HEALPIX_v_MRH_ConstructionTests.h"

std::vector<std::vector<float>> GenRandData(int nPoints) {
	std::vector<float> nextPoint;
	std::vector<std::vector<float>> randPoints;
    
	float latitude,longitude,measurement;
	srand(time(NULL));

	for( int i = 0; i < nPoints; i++)
	{
	   nextPoint.clear();
	   latitude = 180.0*(float(rand())/float(RAND_MAX));
	   longitude = 360*(float(rand())/float(RAND_MAX));
	   measurement = float (rand())/float(RAND_MAX);	  
	   nextPoint.push_back(latitude);
	   nextPoint.push_back(longitude);
	   nextPoint.push_back(measurement);
       randPoints.push_back(nextPoint);
 	}
    return randPoints;
}

void LoadAndDump(MultiResHpx MRH, Healpix_Map<int> HPX, std::vector<std::vector<float>> data) {

	pointing Pt;
	std::vector<Point2D> pointQuery;
	int valHPX,dataVal;
	int hpxID; //HEALPix index
	std::vector<int> mrhID(3); //MRH index

	std::vector<int> hpxIDs;
	std::vector<std::vector<int>> mrhIDs;
	
	printf("#### Data INPUT to MRH and HPX data structures in the following order ####\n\n");

    printf("Lat, Long, HPXid, MRHid, Val,\n");
	// Load data into MRH and HPX data structures while noting respective tree node location and value
	for(unsigned int i = 0; i < data.size(); i++ ) 
	{
	   Pt.theta = data[i][0]*degr2rad; //convert Latitude to theta in Radians
	   Pt.phi = data[i][1]*degr2rad; //convert Longitude to phi in Radians

	   // Convert random float to int
	   dataVal = int( 1000.0 * data[i][2] );

	   //Calc HPX idx
       hpxID = HPX.ang2pix(Pt);

	   //Calc MRH idx
	   mrhID = MRH.ang2mrh(Pt);

	   //Assign data point to HPX idx and record that location
	   HPX[hpxID] = dataVal;
	   hpxIDs.push_back(hpxID);

	   //Assign data point to MRH idx and record that location
	   MRH.Insert(Pt,dataVal);
	   mrhIDs.push_back(mrhID);

	   //Output where data point i is assigned in HPX and MRH and the value
	   printf("%5.3f, %5.3f, %d, %d, %d, %d, %d, \n",data[i][0],data[i][1],hpxID,mrhID[0],mrhID[1],mrhID[2],dataVal);
   
	}

    // Dump the current MRH and HPX data structures noting index and data value at index

    printf("#### Data OUTPUT from MRH and HPX data structures after previous insertions ####\n\n");

    printf("HPXid, MRHid, ValHPX, ValMRH,\n");
	for(unsigned int i = 0; i < data.size(); i++ ) 
	{
	   //Get value stored at next HPX index
       valHPX = HPX[hpxIDs[i]];

	   //Get value stored at next MRH index
       pointQuery = MRH(mrhIDs[i]);
	   
	   for( unsigned int j = 0; j < pointQuery.size(); j++) 
	   {
	      //Output what data was stored (and where) in respective HPX and MRH data structures
	      printf("%d, %d, %d, %d, %d, %d, \n",hpxIDs[i],mrhIDs[i][0],mrhIDs[i][1],mrhIDs[i][2],valHPX,pointQuery[j].mapidx);
	   }
	}
}

void SpillageTest() {
	pointing testPt;
	std::vector<std::vector<float>> randPoints;

	printf("\n###############################################");
	printf("\n########## TEST QUADTREE SPILLAGE    ##########");
	printf("\n###############################################\n\n");

    // Generate batch of random data
	randPoints = GenRandData(10);

    // Compute HPX and MRH indices from random data given HPX and MRH data structures
	// of different sizes. Add the random data to HPX and MRH data structures and then dump trees
	// to analyze the structure and contents of the data structures.

	MultiResHpx_Map<int> mrh1(NEST);
	Healpix_Map<int> hpx1(3,NEST);
	printf("Nside: %d Npix: %d MRHSize: %d HPXSize: %d\n\n",hpx1.Nside(),hpx1.Npix(),mrh1.size(),hpx1.Map().size());

	LoadAndDump(mrh1,hpx1,randPoints);


	//MultiResHpx_Map<int> mMRH2(2,NEST);
	//Healpix_Map<int> hpx2(2,NEST);
	//printf("Nside: %d HPXNside: %d Size: %d HPXSize: %d\n",mMRH2.HPX_Handle().Nside(),hpx2.Nside(),mMRH2.size(),hpx2.Map().size());

	//MultiResHpx_Map<int> mMRH3(3,NEST);
	//Healpix_Map<int> hpx3(3,NEST);
	//printf("Nside: %d HPXNside: %d Size: %d HPXSize: %d\n",mMRH3.HPX_Handle().Nside(),hpx3.Nside(),mMRH3.size(),hpx3.Map().size());

	//MultiResHpx_Map<int> mMRH4(4,NEST);
	//Healpix_Map<int> hpx4(4,NEST);
	//printf("Nside: %d HPXNside: %d Size: %d HPXSize: %d\n",mMRH4.HPX_Handle().Nside(),hpx4.Nside(),mMRH4.size(),hpx4.Map().size());

	//MultiResHpx_Map<int> mMRH5(5,NEST);
	//Healpix_Map<int> hpx5(5,NEST);
	//printf("Nside: %d HPXNside: %d Size: %d HPXSize: %d\n",mMRH5.HPX_Handle().Nside(),hpx5.Nside(),mMRH5.size(),hpx5.Map().size());

	//MultiResHpx_Map<int> mMRH6(6,NEST);
	//Healpix_Map<int> hpx6(6,NEST);
	//printf("Nside: %d HPXNside: %d Size: %d HPXSize: %d\n",mMRH6.HPX_Handle().Nside(),hpx6.Nside(),mMRH6.size(),hpx6.Map().size());

	//MultiResHpx_Map<int> mMRH7(7,NEST);
	//Healpix_Map<int> hpx7(7,NEST);
	//printf("Nside: %d HPXNside: %d Size: %d HPXSize: %d\n",mMRH7.HPX_Handle().Nside(),hpx7.Nside(),mMRH7.size(),hpx7.Map().size());

	//MultiResHpx_Map<int> mMRH8(8,NEST);
	//Healpix_Map<int> hpx8(8,NEST);
	//printf("Nside: %d HPXNside: %d Size: %d HPXSize: %d\n",mMRH8.HPX_Handle().Nside(),hpx8.Nside(),mMRH8.size(),hpx8.Map().size());

	//MultiResHpx_Map<int> mMRH9(9,NEST);
	//Healpix_Map<int> hpx9(9,NEST);
	//printf("Nside: %d HPXNside: %d Size: %d HPXSize: %d\n",mMRH9.HPX_Handle().Nside(),hpx9.Nside(),mMRH9.size(),hpx9.Map().size());

	//MultiResHpx_Map<int> mMRH10(10,NEST);
	//Healpix_Map<int> hpx10(10,NEST);
	//printf("Nside: %d HPXNside: %d Size: %d HPXSize: %d\n",mMRH10.HPX_Handle().Nside(),hpx10.Nside(),mMRH10.size(),hpx10.Map().size());

	//MultiResHpx_Map<int> mMRH11(11,NEST);
	//Healpix_Map<int> hpx11(11,NEST);
	//printf("Nside: %d HPXNside: %d Size: %d HPXSize: %d\n",mMRH11.HPX_Handle().Nside(),hpx11.Nside(),mMRH11.size(),hpx11.Map().size());

	//MultiResHpx_Map<int> mMRH12(12,NEST);
	//Healpix_Map<int> hpx12(12,NEST);
	//printf("Nside: %d HPXNside: %d Size: %d HPXSize: %d\n",mMRH12.HPX_Handle().Nside(),hpx12.Nside(),mMRH12.size(),hpx12.Map().size());

	




}