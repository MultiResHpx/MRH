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

class StripType
{
public:
	StripType() { theta1 = 0.0; theta2 = 0.0; };
	StripType(double t1, double t2) { theta1 = t1;theta2 = t2; };
	~StripType() {};
	double theta1;
	double theta2;
};

class NeighborType
{
public:
	NeighborType() { pt.phi = 0.0; pt.theta = 0.0; };
	NeighborType(pointing _pt) { pt = _pt; };
	~NeighborType() {};
	pointing pt;
};

class Measurement 
{
public:
	Measurement(){rec = EMPTY;};

	~Measurement(){};

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

inline bool Measurement::equals(const Measurement &other){
	if ( this->rec == other.rec ) { return true; }
	return false;
}

void GenRandomData(int numPoints,std::vector<Measurement>& ms,double LO_THETA, double HI_THETA, double LO_PHI, double HI_PHI)
{
	// GIS Longitude (Phi) Range: [-180.0,+180.0] degrees, relative to Prime Meridian.
	// GIS Latitude (Theta) Range: [-90.0,+90.0] degrees, relative to Equator.
	// HPX Longitude (Phi) Range: [0.0,360.0] degrees, relative to Prime Meridian and progressing East.
	// HPX Colatitude (Theta) Range: [0.0,180.0] degrees, relative to North Pole and progressing South to
	// South Pole.
	LO_THETA *= D2R; 
	LO_PHI *= D2R;
	HI_THETA *= D2R; 
	HI_PHI *= D2R;
	int MIN_INT = 0; int MAX_INT = 9999;
	double MIN_DOUBLE = -9999.0; double MAX_DOUBLE = 9999.0;
	ms.clear();
	int64 recNum = 0;
	for(int i = 0; i < numPoints; i++) {
		Measurement m;
		m.rec = recNum;

		// Create random point location
		m.pt.phi = LO_PHI + static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / (HI_PHI - LO_PHI)));
		m.pt.theta = LO_THETA + static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / (HI_THETA - LO_THETA)));
	
		// Random measurement data to demonstrate use of user defined Measurement structure.
		m.data1 = MIN_INT + static_cast <int> (rand()) /( static_cast <double> (RAND_MAX/(MAX_INT-MIN_INT)));
		m.data2 = MIN_INT + static_cast <int> (rand()) /( static_cast <double> (RAND_MAX/(MAX_INT-MIN_INT)));
		m.data3 = MIN_INT + static_cast <int> (rand()) /( static_cast <double> (RAND_MAX/(MAX_INT-MIN_INT)));
		m.data4 = MIN_INT + static_cast <int> (rand()) /( static_cast <double> (RAND_MAX/(MAX_INT-MIN_INT)));

		m.data1 = MIN_DOUBLE + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(MAX_DOUBLE-MIN_DOUBLE)));
		m.data2 = MIN_DOUBLE + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(MAX_DOUBLE-MIN_DOUBLE)));
		m.data3 = MIN_DOUBLE + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(MAX_DOUBLE-MIN_DOUBLE)));
		m.data4 = MIN_DOUBLE + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(MAX_DOUBLE-MIN_DOUBLE)));

		ms.push_back(m);


		recNum += 1;
	}
}



std::vector<DiscType> CreateRandomDiscQueries(int numqueries)
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
	double LO_THETA = 0.01*D2R; double LO_PHI = 0.01*D2R;
	double HI_THETA = 179.99*D2R; double HI_PHI = 359.99*D2R;
	double MIN_RAD = 1.0*D2R; double MAX_RAD = 45.0*D2R;
	double THETA_BUFF = 20.0*D2R;
	double PHI_BUFF = 20.0*D2R;

	for (int i = 0; i < numqueries; i++)
	{
		validDisc = false;
		while (validDisc == false) {
			// Generate random query center location in HPX Coordinate System
			nextDisc.pt.phi = LO_PHI + static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / (HI_PHI - LO_PHI)));
			nextDisc.pt.theta = LO_THETA + static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / (HI_THETA - LO_THETA)));

			// Generate random query radius, making sure radius doesn't sweep
			// off of allowed latitude,longitude limits!
			nextDisc.radius = MIN_RAD + static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / (MAX_RAD - MIN_RAD)));

			if (nextDisc.pt.phi + nextDisc.radius < HI_PHI - PHI_BUFF &&
				nextDisc.pt.phi - nextDisc.radius > LO_PHI + PHI_BUFF &&
				nextDisc.pt.theta + nextDisc.radius < HI_THETA - THETA_BUFF &&
				nextDisc.pt.theta - nextDisc.radius > LO_THETA + THETA_BUFF) {
				validDisc = true;
				_discs.push_back(nextDisc);
			}
		}
	}
	return _discs;
}

vector<PolyType>  CreateRandomConvexPolyQueries(int numqueries)
{
	// GIS Longitude (Phi) Range: [-180.0,+180.0] degrees, relative to Prime Meridian.
	// GIS Latitude (Theta) Range: [-90.0,+90.0] degrees, relative to Equator.
	// HPX Longitude (Phi) Range: [0.0,360.0] degrees, relative to Prime Meridian and progressing East.
	// HPX Colatitude (Theta) Range: [0.0,180.0] degrees, relative to North Pole and progressing South to
	// South Pole.
	int MIN_PTS = 3; 
	double MIN_ARC = 5.0*D2R; double MAX_ARC,nArc;
	double ttlArc;
	pointing pt;
	vector<DiscType> _discs;
	vector<PolyType> _polys;
	vector<pointing> poly;

	// First get valid query discs (fits in the bounds of HPX coordinate system).
	_discs = CreateRandomDiscQueries(numqueries);
	
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

vector<StripType> CreateRandomStripQueries(int numqueries)
{
	int i;
	double LO_THETA = 0.01*D2R; double LO_PHI = 0.01*D2R;
	double HI_THETA = 179.99*D2R; double HI_PHI = 359.99*D2R;
	ofstream fp;
	std::vector<StripType> _strips;
	StripType s;
	bool done = false;
	double thetaWidth;
	for (i = 0; i < numqueries; i++) {
		s.theta1 = LO_THETA + static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / (HI_THETA - LO_THETA)));
		s.theta2 = LO_THETA + static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / (HI_THETA - LO_THETA)));
		_strips.push_back(s);
	}
	return _strips;
}

vector<NeighborType> CreateRandomNeighborQueries(int numqueries, std::vector<Measurement> m)
{
	vector<NeighborType> _neighbors;

	// Now randomly select numQueries from point data list to use
	// for neighbor queries
	bool QUERYPASS = false;
	bool DONE = false;
	int MIN_INT = 0;
	int MAX_INT = m.size() - 1;
	int numAttempts = 0;
	int numGenQueries = 0;
	int index = 0;
	NeighborType n;
	for (int i = 0; i < numqueries; i++)
	{
		// Random index draw
		index = MIN_INT + rand() % (MAX_INT - MIN_INT + 1);

		// Use drawn Measurement's point location as source of Neighbor Query
		// in hopes that there is a neighboring Measurement
		n.pt.phi = m[index].pt.phi;
		n.pt.theta = m[index].pt.theta;
		_neighbors.push_back(n);
	}


	return _neighbors;
}




int main(int64 argc, char* argv[])

{
	int i;
	double nPhi, nTheta, nRadius;
	std::vector<Measurement> myMeasurements;
	std::vector<Measurement> foundMeasurements;
	std::vector< PolyType > Polys;
	std::vector< DiscType > Discs;
	std::vector< StripType > Strips;
	std::vector< NeighborType > Neighbors;
	int NUMPOINTS, MAXDEPTH, NUMQUERIES, NEIGHBOR_QUERY_RESOLUTION;
	double MIN_PHI, MAX_PHI, MIN_THETA, MAX_THETA;

	// Parse command line arguments: num data points, max tree depth, num queries
	if (argc != 8)
	{
		// Print out command line instructions
		std::cout << "\n*** Example.exe Usage ***\n\n";
		std::cout << "Arg 1: Number of random data points to generate\n";
		std::cout << "Arg 2,3: Min. and Max. HPX Longitude (Phi) Range: [0.0,360.0] degrees, relative to Prime Meridian and progressing East.\n";
		std::cout << "Arg 4,5: Min. and Max. HPX Colatitude (Theta) Range: [0.0,180.0] degrees, relative to North Pole and progressing South to the South Pole.\n";
		std::cout << "Arg 6: Max. MRH Tree Depth (1-29)\n";
		std::cout << "Arg 7: Number of Random Disc & Polygon Queries\n\n";
		std::cout << "To run Example.exe with 1000 random data points,\n";
		std::cout << "in a range of 30-40 degrees HPX Phi,\n";
		std::cout << "and in a range of 90-100 degrees HPX Theta,";
		std::cout << "at a maximum tree depth of 4 and query the\n";
		std::cout << "MRH data structure with 5 random disc and polygon\n";
		std::cout << "queries the user would type:\n\n";
		std::cout << "Example.exe 1000 30 40 90 100 4 5\n";
		 
		exit(1);
	}
	else
	{
		NUMPOINTS = atoi(argv[1]);
		MIN_PHI = atof(argv[2]);
		MAX_PHI = atof(argv[3]);
		MIN_THETA = atof(argv[4]);
		MAX_THETA = atof(argv[5]);
		MAXDEPTH = atoi(argv[6]);
		NUMQUERIES = atoi(argv[7]);
	}
	NEIGHBOR_QUERY_RESOLUTION = int(MAXDEPTH / 2);

	MultiResHpx_Map<Measurement> mMRH(MAXDEPTH, NEST);


	// Step 1: Create some random "Measurement" data with random Lat. Long. point positions for MRH data structure
	std::cout << "\nStep 1: Create some random ""Measurement"" data with random Lat. Long. point positions for MRH data structure.\n";
	GenRandomData(NUMPOINTS, myMeasurements,MIN_PHI,MAX_PHI,MIN_THETA,MAX_THETA);

	// Step 2: Insert the Measurements into MRH data structure.
	std::cout << "\nStep 2: Insert the Measurements into MRH data structure.\n";
	for (i = 0; i < NUMPOINTS; i++)
	{
		// Insert next GIS longitude,latitude, data index tuple into MRH
		mMRH.AddRecord(myMeasurements[i], myMeasurements[i].pt);
		cout << "Inserting Point " << i + 1 << " of " << NUMPOINTS << endl;
	}

	// Step 3: Print out the MRH data structure's quad trees
	std::cout << "\nStep 3: Print out the MRH data structure's quad trees.\n";
	for (i = 0; i < 12; i++)
	{
		cout << "\n\n######################\n";
		cout << "#### BASE CELL " << i << " ####\n";
		cout << "######################" << endl;
		mMRH.PrintTreeAtIndex(i);
	}

	// Step 4: Create list of random Disc, Convex Polygon, Strip and Neighbor queries
	std::cout << "\nStep 4: Create list of random Disc and Polygon queries\n";
	Discs = CreateRandomDiscQueries(NUMQUERIES);
	Polys = CreateRandomConvexPolyQueries(NUMQUERIES);
	Strips = CreateRandomStripQueries(NUMQUERIES);
	Neighbors = CreateRandomNeighborQueries(NUMQUERIES, myMeasurements);
	
	// Step 5: Run the Disc Queries the MRH data structure
	std::cout << "\nStep 5: Run the Disc Queries the MRH data structure.\n";
	for (i = 0; i < NUMQUERIES; i++)
	{
		// Report Disc Query Results
		foundMeasurements = mMRH.QueryDisc(Discs[i].pt, Discs[i].radius);
		std::cout << "\n\nDisc Query #" << i + 1 << " definition:\n";
		std::cout << "Phi,Theta,P,Z,Radius: " << Discs[i].pt.phi << "," << Discs[i].pt.theta << ","
			<< Discs[i].pt.phi / pi << "," << cos(Discs[i].pt.theta) << "," << Discs[i].radius << "\n\n";
		std::cout << "Found Measurements: \n";
		for (unsigned int j = 0; j < foundMeasurements.size(); j++)
		{
			std::cout << "Phi,Theta,P,Z,Rec: " << foundMeasurements[j].pt.phi << "," << foundMeasurements[j].pt.theta << ","
				<< foundMeasurements[j].pt.phi / pi << "," << cos(foundMeasurements[j].pt.theta) << ","
				<< foundMeasurements[j].rec << "\n";
		}
	}

	// Step 6: Run the Polygon Queries the MRH data structure
	std::cout << "\nStep 6: Run the Polygon Queries on the MRH data structure.\n";
	for (i = 0; i < NUMQUERIES; i++)
	{
		// Report Polygon Query Results
		foundMeasurements = mMRH.QueryPolygon(Polys[i].pts);
		std::cout << "\n\nPoly Query #" << i + 1 << " definition:\n";
		std::cout << "Phi,Theta,P,Z,...\n";
		for (unsigned j = 0; j < Polys[i].pts.size(); j++) {
			std::cout << Polys[i].pts[j].phi << "," << Polys[i].pts[j].theta << ","
				<< Polys[i].pts[j].phi / pi << "," << cos(Polys[i].pts[j].theta) << "\n";
		}
		std::cout << Polys[i].pts[0].phi << "," << Polys[i].pts[0].theta << ","
			<< Polys[i].pts[0].phi / pi << "," << cos(Polys[i].pts[0].theta) << "\n";
		std::cout << "Found Measurements: \n";
		for (unsigned int j = 0; j < foundMeasurements.size(); j++)
		{
			std::cout << "Phi,Theta,P,Z,Rec: " << foundMeasurements[j].pt.phi << "," << foundMeasurements[j].pt.theta << ","
				<< foundMeasurements[j].pt.phi / pi << "," << cos(foundMeasurements[j].pt.theta) << ","
				<< foundMeasurements[j].rec << "\n";
		}
	}

	// Step 7: Run the Strip Queries on the MRH data structure
	std::cout << "\nStep 7: Run the Strip Queries on the MRH data structure.\n";
	for (i = 0; i < NUMQUERIES; i++)
	{
		// Report Strip Query Results
		foundMeasurements = mMRH.QueryStrip(Strips[i].theta1,Strips[i].theta2);
		std::cout << "\n\nStrip Query #" << i + 1 << " definition:\n";
		std::cout << "Theta1,Theta2,Z1,Z2,...\n";
		std::cout << Strips[i].theta1 << "," << Strips[i].theta2 << ","
	  			  << cos(Strips[i].theta1) << "," << cos(Strips[i].theta2) << "\n";
		std::cout << "Found Measurements: \n";
		for (unsigned int j = 0; j < foundMeasurements.size(); j++)
		{
			std::cout << "Phi,Theta,P,Z,Rec: " << foundMeasurements[j].pt.phi << "," << foundMeasurements[j].pt.theta << ","
				<< foundMeasurements[j].pt.phi / pi << "," << cos(foundMeasurements[j].pt.theta) << ","
				<< foundMeasurements[j].rec << "\n";
		}
	}

	// Step 8: Run the Neighbor Queries on the MRH data structure
	std::cout << "\nStep 8: Run the Neighbor Queries on the MRH data structure.\n";
	for (i = 0; i < NUMQUERIES; i++)
	{
		// Report Neighbor Query Results
		foundMeasurements = mMRH.Neighbors(Neighbors[i].pt, NEIGHBOR_QUERY_RESOLUTION);
		std::cout << "\n\nNeighbor Query #" << i + 1 << " definition:\n";
		std::cout << "Phi,Theta,Depth\n";
		std::cout << Neighbors[i].pt.phi << "," << Neighbors[i].pt.theta << ","
			<< NEIGHBOR_QUERY_RESOLUTION << "\n";
		std::cout << "Found Measurements: \n";
		for (unsigned int j = 0; j < foundMeasurements.size(); j++)
		{
			std::cout << "Phi,Theta,P,Z,Rec: " << foundMeasurements[j].pt.phi << "," << foundMeasurements[j].pt.theta << ","
				<< foundMeasurements[j].pt.phi / pi << "," << cos(foundMeasurements[j].pt.theta) << ","
				<< foundMeasurements[j].rec << "\n";
		}
	}

	// Step 9: Archive the MRH data structure
	std::cout << "\nStep 9: Archive the MRH data structure.\n";
	mMRH.SaveMapToArchive("testArchivedMRH");

	// Step 10: Restore the archived MRH data structure
	std::cout << "\nStep 10: Restore the archived MRH data structure.\n";
	MultiResHpx_Map<Measurement> new_mMRH(MAXDEPTH, NEST);
	new_mMRH.LoadMapFromArchive("testArchiveMRH");

	// Step 11: Re-run Disc queries on restored MRH data structure for V&V purposes
	std::cout << "\nStep 11: Re-run Disc Queries on restored MRH data structure for V&V purposes.\n";
	for (i = 0; i < NUMQUERIES; i++)
	{
		foundMeasurements = mMRH.QueryDisc(Discs[i].pt, Discs[i].radius);
		std::cout << "\n\nDisc Query #" << i + 1 << " definition:\n";
		std::cout << "Phi,Theta,P,Z,Radius: " << Discs[i].pt.phi << "," << Discs[i].pt.theta << ","
			<< Discs[i].pt.phi / pi << "," << cos(Discs[i].pt.theta) << "," << Discs[i].radius << "\n\n";
		std::cout << "Found Measurements: \n";
		for (unsigned int j = 0; j < foundMeasurements.size(); j++)
		{
			std::cout << "Phi,Theta,P,Z,Rec: " << foundMeasurements[j].pt.phi << "," << foundMeasurements[j].pt.theta << ","
				<< foundMeasurements[j].pt.phi / pi << "," << cos(foundMeasurements[j].pt.theta) << ","
				<< foundMeasurements[j].rec << "\n";
		}
	}

	// Step 12: Re-run the Polygon Queries the MRH data structure
	std::cout << "\nStep 12: Re-run Polygon Queries on restored MRH data structure for V&V purposes.\n";
	for (i = 0; i < NUMQUERIES; i++)
	{
		foundMeasurements = mMRH.QueryPolygon(Polys[i].pts);
		std::cout << "\n\nPoly Query #" << i + 1 << " definition:\n";
		std::cout << "Phi,Theta,P,Z,...\n";
		for (unsigned j = 0; j < Polys[i].pts.size(); j++) {
			std::cout << Polys[i].pts[j].phi << "," << Polys[i].pts[j].theta << ","
				<< Polys[i].pts[j].phi / pi << "," << cos(Polys[i].pts[j].theta) << "\n";
		}
		std::cout << Polys[i].pts[0].phi << "," << Polys[i].pts[0].theta << ","
			<< Polys[i].pts[0].phi / pi << "," << cos(Polys[i].pts[0].theta) << "\n";
		std::cout << "Found Measurements: \n";
		for (unsigned int j = 0; j < foundMeasurements.size(); j++)
		{
			std::cout << "Phi,Theta,P,Z,Rec: " << foundMeasurements[j].pt.phi << "," << foundMeasurements[j].pt.theta << ","
				<< foundMeasurements[j].pt.phi / pi << "," << cos(foundMeasurements[j].pt.theta) << ","
				<< foundMeasurements[j].rec << "\n";
		}
	}

	// Step 13: Re-run the Strip Queries on the MRH data structure
	std::cout << "\nStep 13: Re-run the Strip Queries on the MRH data structure for V&V purposes.\n";
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
			std::cout << "Phi,Theta,P,Z,Rec: " << foundMeasurements[j].pt.phi << "," << foundMeasurements[j].pt.theta << ","
				<< foundMeasurements[j].pt.phi / pi << "," << cos(foundMeasurements[j].pt.theta) << ","
				<< foundMeasurements[j].rec << "\n";
		}
	}

	// Step 14: Re-run the Neighbor Queries on the MRH data structure
	std::cout << "\nStep 14: Re-run the Neighbor Queries on the MRH data structure.\n";
	for (i = 0; i < NUMQUERIES; i++)
	{
		// Report Neighbor Query Results
		foundMeasurements = mMRH.Neighbors(Neighbors[i].pt, 1);
		std::cout << "\n\nNeighbor Query #" << i + 1 << " definition:\n";
		std::cout << "Phi,Theta,Depth\n";
		std::cout << Neighbors[i].pt.phi << "," << Neighbors[i].pt.theta << ","
			<< NEIGHBOR_QUERY_RESOLUTION << "\n";
		std::cout << "Found Measurements: \n";
		for (unsigned int j = 0; j < foundMeasurements.size(); j++)
		{
			std::cout << "Phi,Theta,P,Z,Rec: " << foundMeasurements[j].pt.phi << "," << foundMeasurements[j].pt.theta << ","
				<< foundMeasurements[j].pt.phi / pi << "," << cos(foundMeasurements[j].pt.theta) << ","
				<< foundMeasurements[j].rec << "\n";
		}
	}
}

