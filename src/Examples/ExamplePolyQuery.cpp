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


class DiscType
{
public:
	DiscType() { pt.phi = 0.0;pt.theta = 0.0;radius = 0.0; };
	DiscType(pointing p, double r) { pt = p;radius = r; };
	~DiscType() {};
	pointing pt;
	double radius;
};

class PolyType
{
public:
	PolyType() { pts.clear(); };
	PolyType(std::vector<pointing> p) { pts = p; };
	~PolyType() {};
	std::vector<pointing> pts;
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



std::vector<DiscType> CreateRandomDiscQueries(int numqueries, double MIN_PHI, double MAX_PHI, double MIN_THETA, double MAX_THETA, double MIN_RAD, double MAX_RAD)
{
	// GIS Longitude (Phi) Range: [-180.0,+180.0] degrees, relative to Prime Meridian.
	// GIS Latitude (Theta) Range: [-90.0,+90.0] degrees, relative to Equator.
	// HPX Longitude (Phi) Range: [0.0,360.0] degrees, relative to Prime Meridian and progressing East.
	// HPX Colatitude (Theta) Range: [0.0,180.0] degrees, relative to North Pole and progressing South to
	// South Pole.
	std::vector<DiscType> _discs;
	DiscType nextDisc;
	pointing pt;
	bool validDisc = false;

	for (int i = 0; i < numqueries; i++)
	{
		validDisc = false;
		while (validDisc == false) {
			// Generate random query center location in HPX Coordinate System
			nextDisc.pt.phi = MIN_PHI + static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / (MAX_PHI - MIN_PHI)));
			nextDisc.pt.theta = MIN_THETA + static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / (MAX_THETA - MIN_THETA)));

			// Generate random query radius, making sure radius doesn't sweep
			// off of allowed latitude,longitude limits!
			nextDisc.radius = MIN_RAD + static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / (MAX_RAD - MIN_RAD)));

			if (nextDisc.pt.phi + nextDisc.radius < MAX_PHI  &&
				nextDisc.pt.phi - nextDisc.radius > MIN_PHI  &&
				nextDisc.pt.theta + nextDisc.radius < MAX_THETA  &&
				nextDisc.pt.theta - nextDisc.radius > MIN_THETA) {
				validDisc = true;
				_discs.push_back(nextDisc);
			}
		}
	}
	return _discs;
}

vector<PolyType>  CreateRandomConvexPolyQueries(int numqueries, double MIN_PHI, double MAX_PHI, double MIN_THETA, double MAX_THETA, double MIN_RAD, double MAX_RAD)
{
	// GIS Longitude (Phi) Range: [-180.0,+180.0] degrees, relative to Prime Meridian.
	// GIS Latitude (Theta) Range: [-90.0,+90.0] degrees, relative to Equator.
	// HPX Longitude (Phi) Range: [0.0,360.0] degrees, relative to Prime Meridian and progressing East.
	// HPX Colatitude (Theta) Range: [0.0,180.0] degrees, relative to North Pole and progressing South to
	// South Pole.
	int MIN_PTS = 3;
	double MIN_ARC = 5.0*D2R; double MAX_ARC, nArc;
	double ttlArc;
	pointing pt;
	vector<DiscType> _discs;
	vector<PolyType> _polys;
	vector<pointing> poly;

	// First get valid query discs (fits in the bounds of HPX coordinate system).
	_discs = CreateRandomDiscQueries(numqueries, MIN_PHI, MAX_PHI, MIN_THETA, MAX_THETA, MIN_RAD, MAX_RAD);

	for (int i = 0; i < numqueries; i++)
	{

		// Generate random maximum arc
		MAX_ARC = MIN_ARC + static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / (180.0*D2R - MIN_ARC)));

		// Generate an N sided, clockwise, around query disc
		// convex polyton. Basically keep adding sides to the polyton until
		// arc excedes 2PI (full circle). Coordinates of poly are in HPX.
		ttlArc = 0.0;
		poly.clear();
		while (ttlArc < TWOPI) {
			// Compute points on disc
			pt.phi = _discs[i].pt.phi + _discs[i].radius * cos(ttlArc);
			pt.theta = _discs[i].pt.theta + _discs[i].radius * sin(ttlArc);
			poly.push_back(pt);

			// Compute next arc length
			nArc = MIN_ARC + static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / (MAX_ARC - MIN_ARC)));
			ttlArc += nArc;
		}
		_polys.push_back(poly);
	}

	return _polys;
}



int main(int64 argc, char* argv[])

{
	int i;
	double nPhi, nTheta, nRadius;
	std::vector<Measurement> myMeasurements;
	std::vector<Measurement> foundMeasurements;
	std::vector< PolyType > Polys;
	MortonNode m;
	int NUMPOINTS, MAXDEPTH, NUMQUERIES;
	double MIN_PHI, MAX_PHI, MIN_THETA, MAX_THETA, MIN_RAD, MAX_RAD;

	srand(time(NULL));

	// Parse command line arguments: num data points, max tree depth, num queries
	if (argc != 10)
	{
		// Print out command line instructions
		std::cout << "\n*** ExamplePolyQuery.exe Usage ***\n\n";
		std::cout << "Arg 1: Number of random data points to generate\n";
		std::cout << "Arg 2,3: Min. and Max. HPX Longitude (Phi) Range: [0.0,360.0] degrees, relative to Prime Meridian and progressing East.\n";
		std::cout << "Arg 4,5: Min. and Max. HPX Colatitude (Theta) Range: [0.0,180.0] degrees, relative to North Pole and progressing South to the South Pole.\n";
		std::cout << "Arg 6,7: Min. and Max. Query Radii for Disc Queries in degrees.\n";
		std::cout << "Arg 8: Max. MRH Tree Depth (1-29)\n";
		std::cout << "Arg 9: Number of Random Polygon Queries\n\n";
		std::cout << "To run Example.exe with 1000 random data points,\n";
		std::cout << "in a range of 30-40 degrees HPX Phi,\n";
		std::cout << "and in a range of 90-100 degrees HPX Theta,\n";
		std::cout << "with disc query radius between 0.1 and 20.0 degrees,\n";
		std::cout << "at a maximum tree depth of 4 and query the\n";
		std::cout << "MRH data structure with 5 random polygon queries the user would type:\n\n";
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

	//Create list of random Convex Polygon queries
	std::cout << "\nCreate list of random Polygon Queries\n";
	Polys = CreateRandomConvexPolyQueries(NUMQUERIES, MIN_PHI, MAX_PHI, MIN_THETA, MAX_THETA, MIN_RAD, MAX_RAD);

	//Run the Polygon Queries the MRH data structure
	std::cout << "\nRun the Polygon Queries on the MRH data structure.\n";
	for (i = 0; i < NUMQUERIES; i++)
	{
		std::cout << "\n\nPoly Query #" << i + 1 << " definition:\n";
		std::cout << "Phi,Theta,P,Z,...\n";
		for (unsigned j = 0; j < Polys[i].pts.size(); j++) {
			std::cout << Polys[i].pts[j].phi << "," << Polys[i].pts[j].theta << ","
				<< Polys[i].pts[j].phi / pi << "," << cos(Polys[i].pts[j].theta) << "\n";
		}
		std::cout << Polys[i].pts[0].phi << "," << Polys[i].pts[0].theta << ","
			<< Polys[i].pts[0].phi / pi << "," << cos(Polys[i].pts[0].theta) << "\n";

		// Report Polygon Query Results
		foundMeasurements = mMRH.QueryPolygon(Polys[i].pts);
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

