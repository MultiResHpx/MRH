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


class TwoPointType
{
public:
	TwoPointType() { radius = 0.0; };
	TwoPointType( double r) {radius = r; };
	~TwoPointType() {};
	double radius;
};

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



std::vector<TwoPointType> CreateRandomTwoPointQueries(int numqueries, double MIN_RAD, double MAX_RAD)
{
	std::vector<TwoPointType> _2pts;
	TwoPointType next2pt;

	for (int i = 0; i < numqueries; i++)
	{
			// Generate random query radius
		next2pt.radius = MIN_RAD + static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / (MAX_RAD - MIN_RAD)));
		_2pts.push_back(next2pt);
	}
	return _2pts;
}


int main(int64 argc, char* argv[])

{
	int i;
	double nPhi, nTheta, nRadius;
	std::vector<Measurement> myMeasurements;
	std::vector<pair<Measurement,Measurement>> foundMeasurements;
	std::vector< TwoPointType > TwoPts;
	MortonNode m,m2;
	int NUMPOINTS, MAXDEPTH, NUMQUERIES;
	double MIN_PHI, MAX_PHI, MIN_THETA, MAX_THETA, MIN_RAD, MAX_RAD;

	srand(time(NULL));

	// Parse command line arguments: num data points, max tree depth, num queries
	if (argc != 10)
	{
		// Print out command line instructions
		std::cout << "\n*** Example.exe Usage ***\n\n";
		std::cout << "Arg 1: Number of random data points to generate\n";
		std::cout << "Arg 2,3: Min. and Max. HPX Longitude (Phi) Range: [0.0,360.0] degrees, relative to Prime Meridian and progressing East.\n";
		std::cout << "Arg 4,5: Min. and Max. HPX Colatitude (Theta) Range: [0.0,180.0] degrees, relative to North Pole and progressing South to the South Pole.\n";
		std::cout << "Arg 6,7: Min. and Max. Query Radii for Two Point Correlation Queries in degrees.\n";
		std::cout << "Arg 8: Max. MRH Tree Depth (1-29)\n";
		std::cout << "Arg 9: Number of Random Two Point Corr. Queries\n\n";
		std::cout << "To run Example.exe with 1000 random data points,\n";
		std::cout << "in a range of 30-40 degrees HPX Phi,\n";
		std::cout << "and in a range of 90-100 degrees HPX Theta,\n";
		std::cout << "with two point corr. query radius between 0.1 and 20.0 degrees,\n";
		std::cout << "at a maximum tree depth of 4 and query the\n";
		std::cout << "MRH data structure with 5 random two point corr. queries the user would type:\n\n";
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
		MIN_RAD = atof(argv[6])*D2R;
		MAX_RAD = atof(argv[7])*D2R;
		MAXDEPTH = atoi(argv[8]);
		NUMQUERIES = atoi(argv[9]);
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

	//Create list of random two point correlation queries
	std::cout << "\nCreate list of random Two Point Correlation queries\n";
	TwoPts = CreateRandomTwoPointQueries(NUMQUERIES, MIN_RAD, MAX_RAD);

	//Run the Two Point Queries the MRH data structure
	std::cout << "\nRun the Two Point Queries on the MRH data structure.\n";
	for (i = 0; i < NUMQUERIES; i++)
	{
		std::cout << "\n\nTwo Point Corr. Query #" << i + 1 << " definition:\n";
		std::cout << "Radius: " << TwoPts[i].radius << "\n\n";

		// Report Disc Query Results
		foundMeasurements = mMRH.TwoPointCorrBin(TwoPts[i].radius);

		std::cout << "Found Measurement Pairs: \n";
		for (unsigned int j = 0; j < foundMeasurements.size(); j++)
		{
			mMRH.GetMortonNodeAtDataIndex(foundMeasurements[j].first.rec, m);
			mMRH.GetMortonNodeAtDataIndex(foundMeasurements[j].second.rec, m2);
			std::cout << "\t#1 Morton,Phi,Theta,P,Z,Rec: "; PrintMorton(m.m);
			std::cout << "," << foundMeasurements[j].first.pt.phi << "," << foundMeasurements[j].first.pt.theta << ","
				<< foundMeasurements[j].first.pt.phi / pi << "," << cos(foundMeasurements[j].first.pt.theta) << ","
				<< foundMeasurements[j].first.rec << "\n";
			std::cout << "\t#2 Morton,Phi,Theta,P,Z,Rec: "; PrintMorton(m2.m);
			std::cout << "," << foundMeasurements[j].second.pt.phi << "," << foundMeasurements[j].second.pt.theta << ","
				<< foundMeasurements[j].second.pt.phi / pi << "," << cos(foundMeasurements[j].second.pt.theta) << ","
				<< foundMeasurements[j].second.rec << "\n";
			std::cout << " Radial Distance Between: " <<
				RadialDist(pointing(foundMeasurements[j].first.pt.theta, foundMeasurements[j].first.pt.phi),
					pointing(foundMeasurements[j].second.pt.theta, foundMeasurements[j].second.pt.phi))
					<< "\n\n\n";
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

