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
#define MAXDEPTH 4
#define NUMPOINTS 1000
#define NUMQUERIES 10

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

void GenRandomData(int numPoints,std::vector<Measurement>& ms)
{
	// GIS Longitude (Phi) Range: [-180.0,+180.0] degrees, relative to Prime Meridian.
	// GIS Latitude (Theta) Range: [-90.0,+90.0] degrees, relative to Equator.
	// HPX Longitude (Phi) Range: [0.0,360.0] degrees, relative to Prime Meridian and progressing East.
	// HPX Colatitude (Theta) Range: [0.0,180.0] degrees, relative to North Pole and progressing South to
	// South Pole.
	double LO_THETA = 0.001*D2R; double LO_PHI = 0.001*D2R;
	double HI_THETA = 179.999*D2R; double HI_PHI = 359.999*D2R;
	int MIN_INT = 0; int MAX_INT = 9999;
	double MIN_DOUBLE = -9999.0; double MAX_DOUBLE = 9999.0;
	Measurement m;
	ms.clear();
	int64 recNum = 0;
	for(int i = 0; i < numPoints; i++) {
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



std::vector<DiscType> RandomQueryDiscs(int numqueries)
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

vector<PolyType>  RandomConvexPolys(int numqueries)
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
	_discs = RandomQueryDiscs(numqueries);
	
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
	std::vector< DiscType > Discs;
	MultiResHpx_Map<Measurement> mMRH(MAXDEPTH, NEST);

	// Step 1: Create some random "Measurement" data with random Lat. Long. point positions for MRH data structure
	std::cout << "\nStep 1: Create some random ""Measurement"" data with random Lat. Long. point positions for MRH data structure.\n";
	GenRandomData(NUMPOINTS, myMeasurements);

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

	// Step 4: Create list of random Disc and Polygon queries
	std::cout << "\nStep 4: Create list of random Disc and Polygon queries\n";
	Discs = RandomQueryDiscs(NUMQUERIES);
	Polys = RandomConvexPolys(NUMQUERIES);
		
	// Step 5: Run the Disc Queries the MRH data structure
	std::cout << "\nStep 5: Run the Disc Queries the MRH data structure.\n";
	for (i = 0; i < NUMQUERIES; i++)
	{
		foundMeasurements = mMRH.QueryDisc(Discs[i].pt, Discs[i].radius);
		std::cout << "\n\nQuery #" << i + 1 << " definition:\n";
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
	std::cout << "\nStep 6: Run the Polygon Queries the MRH data structures.\n";
	for (i = 0; i < NUMQUERIES; i++)
	{
		foundMeasurements = mMRH.QueryPolygon(Polys[i].pts);
		std::cout << "\n\nQuery #" << i + 1 << " definition:\n";
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

	// Step 7: Archive the MRH data structure
	std::cout << "\nStep 7: Archive the MRH data structure.\n";
	mMRH.SaveMapToArchive("testArchivedMRH");

	// Step 8: Restore the archived MRH data structure
	std::cout << "\nStep 8: Restore the archived MRH data structure.\n";
	MultiResHpx_Map<Measurement> new_mMRH(MAXDEPTH, NEST);
	new_mMRH.LoadMapFromArchive("testArchiveMRH");

	// Step 9: Re-run Disc queries on restored MRH data structure for V&V purposes
	std::cout << "\nStep 9: Re-run Disc Queries on restored MRH data structure for V&V purposes.\n";
	for (i = 0; i < NUMQUERIES; i++)
	{
		foundMeasurements = mMRH.QueryDisc(Discs[i].pt, Discs[i].radius);
		std::cout << "\n\nQuery #" << i + 1 << " definition:\n";
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

	// Step 10: Run the Polygon Queries the MRH data structure
	std::cout << "\nStep 10: Re-run Polygon Queries on restored MRH data structure for V&V purposes.\n";
	for (i = 0; i < NUMQUERIES; i++)
	{
		foundMeasurements = mMRH.QueryPolygon(Polys[i].pts);
		std::cout << "\n\nQuery #" << i + 1 << " definition:\n";
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
}

