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

#define EMPTY -99999
#define TWOPI 6.2831853071796
#ifdef _WIN32
#include <Windows.h>
#else
#include <sys/time.h>
#include <ctime>
#endif

/* Remove if already defined */
typedef long long int64; typedef unsigned long long uint64;


class Measurement
{
public:
	Measurement() { rec = EMPTY; };

	~Measurement() {};

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

inline bool Measurement::equals(const Measurement &other) {
	if (this->rec == other.rec) { return true; }
	return false;
}

void GenRandomData(int numPoints, std::vector<Measurement>& ms, double MIN_PHI, double MAX_PHI, double MIN_THETA, double MAX_THETA)
{
	// GIS Longitude (Phi) Range: [-180.0,+180.0] degrees, relative to Prime Meridian.
	// GIS Latitude (Theta) Range: [-90.0,+90.0] degrees, relative to Equator.
	// HPX Longitude (Phi) Range: [0.0,360.0] degrees, relative to Prime Meridian and progressing East.
	// HPX Colatitude (Theta) Range: [0.0,180.0] degrees, relative to North Pole and progressing South to
	// South Pole.
	int MIN_INT = 0; int MAX_INT = 9999;
	double MIN_DOUBLE = -9999.0; double MAX_DOUBLE = 9999.0;
	ms.clear();
	int64 recNum = 0;
	for (int i = 0; i < numPoints; i++) {
		Measurement m;
		m.rec = recNum;

		// Create random point location
		m.pt.phi = MIN_PHI + static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / (MAX_PHI - MIN_PHI)));
		m.pt.theta = MIN_THETA + static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / (MAX_THETA - MIN_THETA)));

		// Random measurement data to demonstrate use of user defined Measurement structure.
		m.data1 = MIN_INT + static_cast <int> (rand()) / (static_cast <double> (RAND_MAX / (MAX_INT - MIN_INT)));
		m.data2 = MIN_INT + static_cast <int> (rand()) / (static_cast <double> (RAND_MAX / (MAX_INT - MIN_INT)));
		m.data3 = MIN_INT + static_cast <int> (rand()) / (static_cast <double> (RAND_MAX / (MAX_INT - MIN_INT)));
		m.data4 = MIN_INT + static_cast <int> (rand()) / (static_cast <double> (RAND_MAX / (MAX_INT - MIN_INT)));

		m.data1 = MIN_DOUBLE + static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / (MAX_DOUBLE - MIN_DOUBLE)));
		m.data2 = MIN_DOUBLE + static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / (MAX_DOUBLE - MIN_DOUBLE)));
		m.data3 = MIN_DOUBLE + static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / (MAX_DOUBLE - MIN_DOUBLE)));
		m.data4 = MIN_DOUBLE + static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / (MAX_DOUBLE - MIN_DOUBLE)));

		ms.push_back(m);


		recNum += 1;
	}
}



int main(int64 argc, char* argv[])

{
	int i;
	double nPhi, nTheta, nRadius;
	std::vector<Measurement> myMeasurements;
	std::vector<Measurement> foundMeasurements;
	MortonNode m;
	int NUMPOINTS, MAXDEPTH, NUMQUERIES;
	double MIN_PHI, MAX_PHI, MIN_THETA, MAX_THETA;

	srand(time(NULL));

	// Parse command line arguments: num data points, max tree depth, num queries
	if (argc != 8)
	{
		// Print out command line instructions
		std::cout << "\n*** Example.exe Usage ***\n\n";
		std::cout << "Arg 1: Number of random data points to generate\n";
		std::cout << "Arg 2,3: Min. and Max. HPX Longitude (Phi) Range: [0.0,360.0] degrees, relative to Prime Meridian and progressing East.\n";
		std::cout << "Arg 4,5: Min. and Max. HPX Colatitude (Theta) Range: [0.0,180.0] degrees, relative to North Pole and progressing South to the South Pole.\n";
		std::cout << "Arg 6: Max. MRH Tree Depth (1-29)\n";
		std::cout << "Arg 7: Number of Random Nearest Neighbor Queries\n\n";
		std::cout << "To run Example.exe with 1000 random data points,\n";
		std::cout << "in a range of 30-40 degrees HPX Phi,\n";
		std::cout << "and in a range of 90-100 degrees HPX Theta,\n";
		std::cout << "at a maximum tree depth of 4 and query the\n";
		std::cout << "MRH data structure with 5 nearest neighbor queries the user would type:\n\n";
		std::cout << "Example.exe 1000 30.0 40.0 90.0 100.0 0.1 20.0 4 5\n";

		exit(1);
	}
	else
	{
		NUMPOINTS = atoi(argv[1]);
		MIN_PHI = atof(argv[2])*D2R;
		MAX_PHI = atof(argv[3])*D2R;
		MIN_THETA = atof(argv[4])*D2R;
		MAX_THETA = atof(argv[5])*D2R;
		MAXDEPTH = atoi(argv[6]);
		NUMQUERIES = atoi(argv[7]);
	}

	MultiResHpx_Map<Measurement> mMRH(MAXDEPTH, NEST);


	//Create some random "Measurement" data with random Lat. Long. point positions for MRH data structure
	std::cout << "\nCreate some random ""Measurement"" data with random Lat. Long. point positions for MRH data structure.\n";
	GenRandomData(NUMPOINTS, myMeasurements, MIN_PHI, MAX_PHI, MIN_THETA, MAX_THETA);

	//Insert the Measurements into MRH data structure.
	std::cout << "\nInsert the Measurements into MRH data structure.\n";
	for (i = 0; i < NUMPOINTS; i++)
	{
		// Insert next GIS longitude,latitude, data index tuple into MRH
		mMRH.AddRecord(myMeasurements[i], myMeasurements[i].pt);
		if (i % 100 == 0)
		{
			cout << "Inserting Point " << i + 1 << " of " << NUMPOINTS << endl;
		}
	}

	//Print out the MRH data structure's quad trees
	std::cout << "\nPrint out the MRH data structure's quad trees.\n";
	for (i = 0; i < 12; i++)
	{
		cout << "\n\n######################\n";
		cout << "#### BASE CELL " << i << " ####\n";
		cout << "######################" << endl;
		mMRH.PrintTreeAtIndex(i);
	}

	//Create list of random Nearest Neighbor queries
	std::cout << "\nCreate list of random Nearest Neighbor Queries\n";

	// Run the Nearest Neighbor Queries on the MRH data structure
	std::cout << "\nRun the Nearest Neighbor Queries on the MRH data structure.\n";
	int64 index;
	for (i = 0; i < NUMQUERIES; i++)
	{
		// Random index draw
		index = rand() % (NUMPOINTS + 1);

		mMRH.GetMortonNodeAtDataIndex(index, m);

		std::cout << "\n\nNearest Neighbor Query #" << i + 1 << " definition:\n";
		std::cout << "Query Morton,Map Index: ";PrintMorton(m.m); std::cout << "," << index << "\n";

		// Report Neighbor Query Results
		foundMeasurements = mMRH.NearNeighbors(index);

		std::cout << "Found Measurements: \n";
		for (unsigned int j = 0; j < foundMeasurements.size(); j++)
		{
			mMRH.GetMortonNodeAtDataIndex(foundMeasurements[j].rec, m);
			std::cout << "Morton,Phi,Theta,P,Z,Rec: "; PrintMorton(m.m);
			std::cout << "," << foundMeasurements[j].pt.phi << "," << foundMeasurements[j].pt.theta << ","
				<< foundMeasurements[j].pt.phi / pi << "," << cos(foundMeasurements[j].pt.theta) << ","
				<< foundMeasurements[j].rec << "\n";
		}
	}

	// Archive the MRH data structure
	std::cout << "\nArchive the MRH data structure.\n";
	mMRH.SaveMapToArchive("testArchivedMRH");

	// Restore the archived MRH data structure
	std::cout << "\nRestore the archived MRH data structure.\n";
	MultiResHpx_Map<Measurement> new_mMRH(MAXDEPTH, NEST);
	new_mMRH.LoadMapFromArchive("testArchiveMRH");
}

