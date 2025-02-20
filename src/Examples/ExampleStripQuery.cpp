/*
 * Copyright (C) 2017  Robert Youngren, robert.youngren@gmail.com
 *
 * This file is part of MultiResHpx.
 *
 * MultiResHpx is free software : you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MultiResHpx is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MultiResHpx.If not, see < https://www.gnu.org/licenses/>.
*/

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


class StripType
{
public:
	StripType() { theta1 = 0.0; theta2 = 0.0; };
	StripType(double t1, double t2) { theta1 = t1;theta2 = t2; };
	~StripType() {};
	double theta1;
	double theta2;
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


vector<StripType> CreateRandomStripQueries(int numqueries, double MIN_THETA, double MAX_THETA, double MIN_STRIP, double MAX_STRIP)
{
	int i;
	ofstream fp;
	std::vector<StripType> _strips;
	StripType s;
	bool done = false;
	double thetaWidth;

	// Generate random strip width;
	double strip_width = MIN_STRIP + static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / (MAX_STRIP - MIN_STRIP)));

	for (i = 0; i < numqueries; i++) {
		s.theta1 = MIN_THETA + static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / (MAX_THETA - MIN_THETA)));
		s.theta2 = s.theta1 + strip_width;
		_strips.push_back(s);
	}
	return _strips;
}

int main(int64 argc, char* argv[])

{
	int i;
	double nPhi, nTheta, nRadius;
	std::vector<Measurement> myMeasurements;
	std::vector<Measurement> foundMeasurements;
	std::vector< StripType > Strips;
	MortonNode m;
	int NUMPOINTS, MAXDEPTH, NUMQUERIES, NEIGHBOR_QUERY_RESOLUTION;
	double MIN_PHI, MAX_PHI, MIN_THETA, MAX_THETA, MIN_RAD, MAX_RAD, MIN_STRIP, MAX_STRIP;

	srand(time(NULL));

	// Parse command line arguments: num data points, max tree depth, num queries
	if (argc != 10)
	{
		// Print out command line instructions
		std::cout << "\n*** ExampleStripQuery.exe Usage ***\n\n";
		std::cout << "Arg 1: Number of random data points to generate\n";
		std::cout << "Arg 2,3: Min. and Max. HPX Longitude (Phi) Range: [0.0,360.0] degrees, relative to Prime Meridian and progressing East.\n";
		std::cout << "Arg 4,5: Min. and Max. HPX Colatitude (Theta) Range: [0.0,180.0] degrees, relative to North Pole and progressing South to the South Pole.\n";
		std::cout << "Arg 6,7: Min. and Max. Strip Widths for Strip Queries in degrees\n";
		std::cout << "Arg 6: Max. MRH Tree Depth (1-29)\n";
		std::cout << "Arg 7: Number of Random Strip Queries\n\n";
		std::cout << "To run Example.exe with 1000 random data points,\n";
		std::cout << "in a range of 30-40 degrees HPX Phi,\n";
		std::cout << "and in a range of 90-100 degrees HPX Theta,\n";
		std::cout << "and strip query widths ranging from 1-5 degrees,\n";
		std::cout << "at a maximum tree depth of 4 and query the\n";
		std::cout << "MRH data structure with 5 random strip queries the user would type:\n\n";
		std::cout << "Example.exe 1000 30.0 40.0 90.0 100.0 1 5 4 5\n";

		exit(1);
	}
	else
	{
		NUMPOINTS = atoi(argv[1]);
		MIN_PHI = atof(argv[2])*D2R;
		MAX_PHI = atof(argv[3])*D2R;
		MIN_THETA = atof(argv[4])*D2R;
		MAX_THETA = atof(argv[5])*D2R;
		MIN_STRIP = atof(argv[6])* D2R;
		MAX_STRIP = atof(argv[7])* D2R;
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

	//Create list of random Strip queries
	std::cout << "\nCreate list of random Strip Queries\n";
	Strips = CreateRandomStripQueries(NUMQUERIES, MIN_THETA, MAX_THETA, MIN_STRIP, MAX_STRIP);

	//Run the Strip Queries on the MRH data structure
	std::cout << "\nRun the Strip Queries on the MRH data structure.\n";
	for (i = 0; i < NUMQUERIES; i++)
	{
		// Report Strip Query Results
		foundMeasurements = mMRH.QueryStrip(Strips[i].theta1, Strips[i].theta2);
		std::cout << "\n\nStrip Query #" << i + 1 << " definition:\n";
		std::cout << "Theta1,Theta2,Z1,Z2,...\n";
		std::cout << Strips[i].theta1 << "," << Strips[i].theta2 << ","
			<< cos(Strips[i].theta1) << "," << cos(Strips[i].theta2) << "\n";
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

