#include "TestApplication.h"

#define TEST_POINTQUERY false
#define TEST_DISCQUERY true
#define TEST_INCLUSIVE_DISCQUERY false
#define TEST_TRIANGLE_QUERY false
#define TEST_TRIANGLE_INCLUSIVE_QUERY false
#define TEST_STRIP_QUERY false
#define TEST_STRIP_INCLUSIVE_QUERY false
#define TEST_NEIGHBORS_QUERY false
#define TEST_INTERPOLATION_QUERY true

//typedef unsigned long long timestamp_t;
//
// static timestamp_t
// get_timestamp ()
// {
//   struct timeval now;
//   gettimeofday (&now, NULL);
//   return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
// }

#ifdef WIN32
#include <Windows.h>
#else
#include <sys/time.h>
#include <ctime>
#endif

/* Returns the amount of milliseconds elapsed since the UNIX epoch. Works on both
 * windows and linux. */

uint64 GetTimeMs64()
{
#ifdef WIN32
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

void DoQueryAndOutputResults(Healpix_Map<int> hpxNEST, MultiResHpx mrh_db,std::vector<Data> db,std::string queryType, std::vector<float> queryParams) {
 //  ofstream out;
	//pointing latlong;
	//float calc_z,calc_phi;
	//out.open("dumpDataLocs.csv");
	//int db_index;
	//bool foundHPX = false;
	//bool foundMRH = false;
	//vector<int> queryResultListHPX;
	//vector<int> queryResultListMRH;
	//rangeset<int> queryResultRange;
	//vector<int> nextCell;

	//// Get the Multi-Res HEALPix indices of populated leaf nodes
	//std::vector< std::vector<int> > leafLocs = mrh_db.getForest();

	//// Output header that details MRH setup
	//for(unsigned int face = 0; face < 12; face++) {
	//	out << "FACE," << face << ",NSIDE," << mrh_db.getNsideAtFace(face) << endl;
	//}
	//out << endl;

	//// Next output the query type and parameters
	//out << "QUERY_TYPE," << queryType << ",";
	//for( unsigned int i=0; i < queryParams.size(); i++ ) {
	//	out << queryParams[i] << ",";
	//}
	//out << endl;

 //  // Next determine the query type and call the respective query method
 //  
	//// DISC QUERY
	//if( strcmp("DISC",queryType.c_str()) == 0 )
	//{
 // 	   pointing testPt;
	//	float radius = queryParams[2];
	//   testPt.phi = queryParams[0];	
	//   testPt.theta = queryParams[1]; 

	//	// Perform HEALPix DISC QUERY
	//	hpxNEST.query_disc(testPt,radius,queryResultRange);
	//	
	//	// Convert disc query results from range form to list form
	//	queryResultRange.toVector(queryResultListHPX);

	//	// Perform MRH DISC QUERY
	//	mrh_db.query_disc_map(testPt,radius,queryResultListMRH);

	//}

	//// Now dump out all the leaf indices to file, marking those indices that fell within the respective
	//// HPX or MRH query.
	//out << "IX,IY,FACE,CELL_LAT,CELL_LONG,CELL_Z,CELL_PHI,HPX_ID,TRUE_LAT,TRUE_LONG,TRUE_Z,TRUE_PHI,IN_HPX_Q,IN_MRH_Q" << endl;
	//for(unsigned int node = 0; node < leafLocs.size(); node++) {
 //     latlong = mrh_db.mrh2ang(leafLocs[node][0],leafLocs[node][1],leafLocs[node][2]);
	//	db_index = mrh_db.at(leafLocs[node][0],leafLocs[node][1],leafLocs[node][2]);
 //     calc_z = cos(latlong.theta);
	//	calc_phi = (latlong.phi/pi);
	//	out << leafLocs[node][0] << ","
	//		 << leafLocs[node][1] << ","
	//		 << leafLocs[node][2] << "," 
	//		 << latlong.theta*rad2degr << "," 
	//		 << latlong.phi*rad2degr << ","
	//		 << calc_z << ","
	//		 << calc_phi << ","
	//		 << db[db_index].hpxId << ","
	//		 << db[db_index].latitude << "," 
	//		 << db[db_index].longitude << "," 
	//		 << db[db_index].z << ","
	//		 << db[db_index].phi;
	//   
 //		// Determine which (if any) data point HPX/MRH indices match up with 
	//	// query HPX/MRH indices and mark them.
	//	foundHPX = false;
 // 	   for(unsigned int i = 0; i < queryResultListHPX.size(); i++ ) {

 //        // Was this HPX node found in the query?
	//		if( db[db_index].hpxId == queryResultListHPX[i] ) {
 //           out << "1" << ",";
	//			foundHPX = true;
	//		}
	//	}
	//	if(foundHPX == false) {
 //        out << "0" << ",";
	//	}
 //     
	//	foundMRH = false;
 // 	   for(unsigned int j = 0; j < queryResultListMRH.size(); j++ ) 
	//	{
	//		nextCell = queryResultListMRH[j];
	//	    // Was this MRH node found in the query?
	//		if( (nextCell[0] == leafLocs[node][0]) && 
	//			 (nextCell[1] == leafLocs[node][1]) && 
	//			 (nextCell[2] == leafLocs[node][2]) ) 
	//		{
	//			out << 1 << ",";
	//			out << "1" << ",";
	//			foundMRH = true;
	//		}
	//	}
	//	if(foundMRH == false) {
	//		out << 0 << ",";
 //        out << "0" << ",";
	//	}

	//	out << endl;

	//}

	//out.close();
}

void SimpleTest1() {
	pointing ptLoc;
	std::vector<int> id(3);
	std::vector<float> queryParams;
	float thetaDeg;
	float phiDeg;
	float radius;
	int hpxID;
	std::vector<Point2D> query;

	printf("\n#######################################################");
	printf("\n### TEST #1 MRH: INSERT KNOWN SPATIAL DATA INTO MRH ###");
	printf("\n#######################################################\n\n");

    // Create simple test set for V&V
	MultiResHpx_Map<Data> mrhNEST(NEST);
	std::vector<Data> vData;
	int x,y;
	srand(time(NULL));
	for( int i = 0; i < 100; i++)
	{
		// Store the incoming data in std::vector<Data>
		Data next;
		next.latitude = 90.0 * (-0.5 + (float(rand())/float(RAND_MAX)) );
		next.longitude = 360.0 * (float(rand())/float(RAND_MAX));
		next.meas1 = float (rand())/float(RAND_MAX);
		next.meas2 = float (rand())/float(RAND_MAX);
		next.meas3 = float (rand())/float(RAND_MAX);
		next.meas4 = float (rand())/float(RAND_MAX);

		// Compute the Multi-Res HEALPix indices
		ptLoc.phi = next.longitude*degr2rad;
		ptLoc.theta = next.latitude*degr2rad;

		// Add data record to MRH
		mrhNEST.AddRecord(next,ptLoc);

        //////////////TEST TEST REMOVE REMOVE
        query = mrhNEST.GetPointsWithMapIdx(i);
		//////////////TEST TEST REMOVE REMOVE


	}

	//// Now perform various queries on the stored data set.

	//// POINT QUERY
	//if( TEST_POINTQUERY ) {

	//}
	//
	//// INTERPOLATION QUERY
	//if( TEST_INTERPOLATION_QUERY ) {

	//}

	//// DISC QUERY
	//if( TEST_DISCQUERY ) {
	//	phiDeg = 180.0; // Phi/PI = 1.0
	//	thetaDeg = 90.0; // cos(theta) = 0.0
	//	radius = 1.0; // in radians
	//	queryParams.clear();
	//	queryParams.push_back(phiDeg*degr2rad); //Longitude or X value of query disc center in radians
	//	queryParams.push_back(thetaDeg*degr2rad); //Latitude or Y value of query disc center in radians
	//	queryParams.push_back(radius); //Radius of query in radians.
 //     DoQueryAndOutputResults(hpxNEST,mMRH2,vData,"DISC",queryParams);  
	//}

	//if( TEST_INCLUSIVE_DISCQUERY ) {

	//}

	//// TRIANGLE QUERY
	//if( TEST_TRIANGLE_QUERY ) {

	//}

	//if( TEST_TRIANGLE_INCLUSIVE_QUERY ) {

	//}


	//// STRIP QUERY
	//if( TEST_STRIP_QUERY ) {

	//}

	//if( TEST_STRIP_INCLUSIVE_QUERY ) {

	//}

	//// NEIGHBORS QUERY
	//if( TEST_NEIGHBORS_QUERY ) {

	//}

}