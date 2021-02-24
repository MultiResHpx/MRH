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

#include "MultiResHpx.h"

#include <windows.h>
#include <stdio.h>
#include <tchar.h>

using namespace std;

MultiResHpx::MultiResHpx
( 
 int Max_Depth,
 Healpix_Ordering_Scheme scheme
)
{
	// Setup template MortonLQT
	// Set maximum depth of MortonLQT w/o having to 
	// share node address
	MortonLQT _tree(Max_Depth); 
    forest_.clear();
	forest_.push_back(_tree); //Base Cell: 0
    forest_[0].SetFaceNum(0);

	forest_.push_back(_tree); //Base Cell: 1 
    forest_[1].SetFaceNum(1);

	forest_.push_back(_tree); //Base Cell: 2 
    forest_[2].SetFaceNum(2);

	forest_.push_back(_tree); //Base Cell: 3 
    forest_[3].SetFaceNum(3);

	forest_.push_back(_tree); //Base Cell: 4 
    forest_[4].SetFaceNum(4);

	forest_.push_back(_tree); //Base Cell: 5 
    forest_[5].SetFaceNum(5);

	forest_.push_back(_tree); //Base Cell: 6 
    forest_[6].SetFaceNum(6);

	forest_.push_back(_tree); //Base Cell: 7 
    forest_[7].SetFaceNum(7);

	forest_.push_back(_tree); //Base Cell: 8 
    forest_[8].SetFaceNum(8);

	forest_.push_back(_tree); //Base Cell: 9 
    forest_[9].SetFaceNum(9);

	forest_.push_back(_tree); //Base Cell: 10 
    forest_[10].SetFaceNum(10);

	forest_.push_back(_tree); //Base Cell: 11 
    forest_[11].SetFaceNum(11);

	hpx_scheme = scheme;
	hpxQ.Set(MIN_ORDER,scheme);
	lowHPX.Set(2,hpx_scheme);
	max_tree_depth = Max_Depth;
	num_rec = 0;
	search_sentinel = 0;
}


MultiResHpx::~MultiResHpx()
{

}

//####
//#### UTILITY METHODS
//####

// Sort ascending by HPX order
bool SortFunctionPixset( pair<int64,int> a, pair<int64,int> b ) {
	if( a.second < b.second )
		return true;
	if( a.second == b.second )
		return a.first < b.first;
	return false;
}

double RadialDist(pointing pt1,pointing pt2)
{
	return acos(fabs(cosdist_zphi(cos(pt1.theta),pt1.phi,cos(pt2.theta),pt2.phi)));
}

void HpxPZCellBoundaries(pair<int64,int> hpxid,int step)
{
	std::vector<vec3> out;
	unsigned int i;
	double p;
    Healpix_Custom hpxQuery(hpxid.second,NEST);	
	hpxQuery.boundaries(hpxid,step,out);
    
	// Convert 3-vector boundary to zphi
    for( i = 0; i < out.size(); i++)
	{
        p = (atan2(out[i].y,out[i].x))/pi;
		if( p < 0.0 ){
			p += 2.0;
		}
		cout << p*pi << "," << acos(out[i].z) << "," << p << "," << out[i].z << endl; 
	}
    p = (atan2(out[0].y,out[0].x))/pi;
	if( p < 0.0 ){
		p += 2.0;
	}
	cout << p*pi << "," << acos(out[0].z) << "," << p << "," << out[0].z << endl; 
}

bool IsPointInPoly2(pointing pt,std::vector<pointing> poly)
{
	double X_MAX = twopi;
	double X_MIN = 0.0; 
	double Y_MAX = pi; //south pole
	double Y_MIN = 0.0; //north pole
	double X_RANGE = twopi;
	double dx,dy,a,xi,xj,xmin,xmax;
	int oddNodes=0;
	int i,j=poly.size()-1;

	for(i=0; i < poly.size(); i++) {
		if( (poly[i].theta < pt.theta && poly[j].theta >= pt.theta) ||
			(poly[j].theta < pt.theta && poly[i].theta >= pt.theta) )
		{
			dx = poly[j].phi - poly[i].phi;
			xi = pt.phi - poly[i].phi;
			xj = pt.phi - poly[j].phi;
			xmax = pt.phi-X_MAX;
			xmin = pt.phi-X_MIN;

			if( 0 < xi && 0 < xj )
			{
				oddNodes = 1-oddNodes;
			}
			else if( !(0 > xi && 0 > xj ) )
			{
				dy = poly[j].theta - poly[i].theta;
				a = pt.theta - poly[i].theta;
				if( dy > 0 ) {
					if( a*dx < xi*dy )
						oddNodes = 1-oddNodes;
				}
				else if( a*dx > xi*dy )
					oddNodes = 1-oddNodes;
			}
		}
		j=i;
	}

	if(oddNodes == 1)
		return true;
	else
		return false;
}


bool IsPointInStrip(double theta1, double theta2, pointing pt) {

	if (theta1<theta2) {
		if( pt.theta >= theta1 && pt.theta <= theta2 ) {
			return true;
		}
	} else {
		if( (pt.theta >= 0.0 && pt.theta <= theta2) ||
			(pt.theta >= theta1 && pt.theta <= pi) ) {
			return true;
		}
	}
	return false;
}


bool IsPointInDisc(pointing center,double radius,pointing pt) {
	double dist = RadialDist(center,pt);
	if(dist <= radius){ return true; }
	return false;
}



bool IsMortonNodePairInList
(
 std::vector<pair<MortonNode,MortonNode>> list,
 pair<MortonNode,MortonNode> nPair
)
{
	for(unsigned int i = 0; i < list.size(); i++ ) 
	{
		// Check for:
		// (x,y) == (x,y)
		if ( Equals(list[i].first.m, nPair.first.m) && Equals(list[i].second.m, nPair.second.m))
		{
			return true;
		}
		// Check for reverse pair order:
		// (x,y) == (y,x)
		if( Equals(list[i].first.m,nPair.second.m ) && Equals(list[i].second.m,nPair.first.m) ) 
		{
			return true;
		}
	}
	return false;
}


//Input: pointing angle = HPX phi(phi)(radians), HPX theta (theta) (radians)
pointing HPXtoGIS(pointing pt)
{
	// GIS Longitude Range: [-180.0,+180.0] degrees, relative to Prime Meridian.
	// GIS Latitude Range: [-90.0,+90.0] degrees, relative to Equator.
	// HPX Longitude Range: [0.0,360.0] degrees, relative to Prime Meridian and progressing East.
	// HPX Colatitude Range: [0.0,180.0] degrees, relative to North Pole and progressing South to
	// South Pole. Input units must be RADIANS.
	pt.phi = pt.phi - pi;
	pt.theta = (pi/2.0)-pt.theta;	
	return pt;
}

//Input: pointing angle = GIS phi(phi)(radians), GIS theta (theta) (radians)
pointing GIStoHPX(pointing pt)
{
	// GIS Longitude Range: [-180.0,+180.0] degrees, relative to Prime Meridian.
	// GIS Latitude Range: [-90.0,+90.0] degrees, relative to Equator.
	// HPX Longitude Range: [0.0,360.0] degrees, relative to Prime Meridian and progressing East.
	// HPX Colatitude Range: [0.0,180.0] degrees, relative to North Pole and progressing South to
	// South Pole. Input units must be RADIANS.
	pt.phi = pt.phi + pi;
	pt.theta = -pt.theta + (pi/2.0);
	return pt;
}

pointing RADECtoHPX(pointing pt)
{
	// (RA) Right Ascension Range: [0.0,360.00] degrees, relative to Prime Meridian and progressing East.
	// (DEC) Declination: [-90.0,+90] degrees, relative to Celestial Equator.
	// HPX Longitude Range: [0.0,360.0] degrees, relative to Prime Meridian and progressing East.
	// HPX Colatitude Range: [0.0,180.0] degrees, relative to North Pole and progressing South to 
	// South Pole. Input units must be RADIANS.
	pt.theta = -pt.theta + (pi/2.0);

	return pt;
}

pointing HPXtoRADEC(pointing pt)
{
	// (RA) Right Ascension Range: [0.0,360.00] degrees, relative to Prime Meridian and progressing East.
	// (DEC) Declination: [-90.0,+90] degrees, relative to Celestial Equator.
	// HPX Longitude Range: [0.0,360.0] degrees, relative to Prime Meridian and progressing East.
	// HPX Colatitude Range: [0.0,180.0] degrees, relative to North Pole and progressing South to 
	// South Pole. Input units must be RADIANS.
	pt.theta = (pi/2.0)-pt.theta;	

	return pt;
}

void Vec2ang(vec3 vec,pointing &pt)
{
	// Convert 3-vector to phi,theta angles
    pt.phi = (atan2(vec.y,vec.x))/pi;
	if( pt.phi < 0.0 ){
		pt.phi += 2.0;
	}
	pt.phi *= pi;
	pt.theta = acos(vec.z);
}

void Ang2vec(pointing &pt, vec3 vec)
{
	// Convert phi,theta angles to 3-vector
	vec.x = cos(pt.phi)*sin(pt.theta);
	vec.y = sin(pt.phi)*sin(pt.theta);
	vec.z = cos(pt.theta);
}

///
/// MultiResHPX Method Definitions
///


std::vector<int> MultiResHpx::GetNodeSizeHistogram()
{
	std::vector<int> node_size_histo;
	std::vector<int> temp;
	node_size_histo.clear();
	node_size_histo.push_back(0);
	unsigned int i,j;

	// For each MortonLQT
	for(i = 0; i < 12; i++) 
	{
		cout << "#### Processing MLQ " << i+1 << " ####\n";
		
		if( forest_[i].GetNumMortonNodes() > 0 ) 
		{
			temp.clear();
			//Get its node size histogram
			temp = forest_[i].GetMortonNodeSizeHistogram();

			// Update the master node size histogram, expanding it if need be...
			if( temp.size() >= node_size_histo.size() )
			{
				while( temp.size() >= node_size_histo.size() )
				{
					node_size_histo.push_back(0);
				}
			}

			for(  j = 0; j < temp.size(); j++)
			{
				node_size_histo[j] += temp[j];
			}
		}
		cout << "\tDONE!\n";
	}
	return node_size_histo;
}


int64 MultiResHpx::NumNodes()
{
	int64 count = 0;
	//Count number of MortonNodes in the MortonLQT Forest
	for(int i = 0; i < 12; i++) {
		count += forest_[i].GetNumMortonNodes();
	}
	return count;
}

int64 MultiResHpx::NumRec()
{
	return num_rec;
}

int64 MultiResHpx::MemSizeKb()
{
	return (NumNodes()*sizeof(MortonNode))/1024; //Convert to Kb
}

int MultiResHpx::MaxDepth()
{
	int maxdepth = 0;
	for(int i = 0; i < 12; i++) {
		if( forest_[i].GetMortonTreeDepth() > maxdepth )
		{
			maxdepth = forest_[i].GetMortonTreeDepth();
		}
	}
	return maxdepth;
}

int MultiResHpx::MinDepth()
{
	int mindepth = 100;
	for(int i = 0; i < 12; i++) {
		if( forest_[i].GetMortonTreeDepth() < mindepth )
		{
			mindepth = forest_[i].GetMortonTreeDepth();
		}
	}
	return mindepth;
}

double MultiResHpx::AvgDepth()
{
	double sumDepth = 0;
	for(int i = 0; i < 12; i++) {
		sumDepth += double(forest_[i].GetMortonTreeDepth());
	}
	return sumDepth/12.0;
}

int64 MultiResHpx::NumNodesAtMLQ(int idx)
{
	return forest_[idx].GetNumMortonNodes();
}

int MultiResHpx::DepthAtMLQ(int idx)
{
	return forest_[idx].GetMortonTreeDepth();
}

double MultiResHpx::AvgDepthAtMLQ(int idx)
{
	return forest_[idx].GetAvgMortonTreeDepth();
}


int MultiResHpx::GetCoverMapSize()
{
	return MRH_COVER_MAP_SIZE;
}

int MultiResHpx::GetCoverMapCellRes()
{
	return MRH_COVER_MAP_CELL_RES;
}



//Inputs: pt, pointing angle = HPX phi(phi)(radians), HPX theta (theta) (radians), reference index to data item
//        data_idx, external reference index to data location.
//        If used in conjunction with MultiResHpx_Map, data_idx is the
//        std::vector<T> index of the element.
void MultiResHpx::Insert(pointing pt, int64 data_idx)
{
	// GIS Longitude Range: [-180.0,+180.0] degrees, relative to Prime Meridian.
	// GIS Latitude Range: [-90.0,+90.0] degrees, relative to Equator.
	// HPX Longitude Range: [0.0,360.0] degrees, relative to Prime Meridian and progressing East.
	// HPX Colatitude Range: [0.0,180.0] degrees, relative to North Pole and progressing South to
	// South Pole. Input units must be RADIANS.
	int fn;
	std::vector<int64> data_indices;
	data_indices.push_back(data_idx);

	// Determine which face number (fn) HPX Longitude,Colatitude pair map to.
	fn = hpxQ.FaceNum(pt);

	// Insert data index into appropriate tree using normalized
	// HEALPix indices
	forest_[fn].InsertMortonNodeAtPhiTheta(pt,data_indices);
	
	// Keep count of total number of records inserted into MRH.
	num_rec++;

}


//Inputs: pointing angle = HPX phi(phi)(radians), HPX theta (theta) (radians)
void MultiResHpx::Delete(pointing pt)
{
	// GIS Longitude Range: [-180.0,+180.0] degrees, relative to Prime Meridian.
	// GIS Latitude Range: [-90.0,+90.0] degrees, relative to Equator.
	// HPX Longitude Range: [0.0,360.0] degrees, relative to Prime Meridian and progressing East.
	// HPX Colatitude Range: [0.0,180.0] degrees, relative to North Pole and progressing South to
	// South Pole. Input units must be RADIANS.
	int fn;
	
	// Determine which face number (fn) HPX Longitude,Colatitude pair map to.
	fn = hpxQ.FaceNum(pt);

	// Delete node from appropriate tree 
	forest_[fn].DeleteMortonNodeAtPhiTheta(pt);
	num_rec--;
}


////#####################################################
////##### PUBLIC TRANSLATION METHODS ####################
////#####################################################

void MultiResHpx::Pix2ring 
(
 int64 hpxid, 
 int order,
 int64& ring
) 
{
	hpxQ.Set(order,hpx_scheme);
	ring = hpxQ.pix2ring (hpxid);
}

void MultiResHpx::Xyf2pix
(
 int ix, 
 int iy, 
 int fn,
 int order,
 int64& hpxid 
) 
{
	hpxQ.Set(order,hpx_scheme);
	hpxid = hpxQ.xyf2pix(ix,iy,fn);
}

void MultiResHpx::Pix2xyf
(
 int64 hpxid,
 int order, 
 int &ix, 
 int &iy, 
 int &fn
) 
{
  hpxQ.Set(order,hpx_scheme);
  hpxQ.pix2xyf(hpxid,ix,iy,fn);
}

    /*! Translates a pixel number from NEST to RING. */
void MultiResHpx::Nest2ring 
(
 int64 hpxid, 
 int order, 
 int64& Rhpxid
) 
{
  hpxQ.Set(order,hpx_scheme);
  Rhpxid = hpxQ.nest2ring(hpxid);
}
   
	//*! Translates a pixel number from RING to NEST. */
void MultiResHpx::Ring2nest 
(
 int64 Rhpxid, 
 int order, 
 int64& Nhpxid
) 
{
  hpxQ.Set(order,hpx_scheme);
  Nhpxid = hpxQ.ring2nest(Rhpxid);
}

	//*! Translates a pixel number from NEST to its Peano index. */
void MultiResHpx::Nest2peano 
(
 int64 hpxid, 
 int order, 
 int64& pid
) 
{
  hpxQ.Set(order,hpx_scheme);
  pid = hpxQ.nest2peano(hpxid);
}

	//*! Translates a pixel number from its Peano index to NEST. */
void MultiResHpx::Peano2nest 
(
 int64 pid, 
 int order, 
 int64& Nhpxid
) 
{
  hpxQ.Set(order,hpx_scheme);
  Nhpxid = hpxQ.peano2nest(pid);
}

    /*! Returns the number of the pixel which contains the angular coordinates
        (\a z:=cos(theta), \a phi).
        \note This method is inaccurate near the poles at high resolutions. */
void MultiResHpx::Zphi2pix 
(
 double z, 
 double phi, 
 int order, 
 int64& hpxid
) 
{
  hpxQ.Set(order,hpx_scheme);
  hpxid = hpxQ.zphi2pix(z,phi);
}

    /*! Returns the number of the pixel which contains the angular coordinates
        \a ang. */
void MultiResHpx::Ang2pix 
(
 pointing &ang, 
 int order, 
 int64& hpxid
) 
{
  //ShiftAngleGIStoHPX(ang,1.0);
  hpxQ.Set(order,hpx_scheme);
  hpxid = hpxQ.ang2pix(ang);
}

    /*! Returns the number of the pixel which contains the vector \a vec
        (\a vec is normalized if necessary). */
void MultiResHpx::Vec2pix 
(
 vec3 &v,
 int order, 
 int64& hpxid
) 
{
  hpxQ.Set(order,hpx_scheme);
  hpxid = hpxQ.ang2pix(v);
}

	/*! Returns the angular coordinates (\a z:=cos(theta), \a phi) of the center
        of the pixel with number \a pix.
        \note This method is inaccurate near the poles at high resolutions. */
void MultiResHpx::Pix2zphi 
(
 int64 hpxid, 
 int order, 
 double &z, 
 double &phi
) 
{
  hpxQ.Set(order,hpx_scheme);
  hpxQ.pix2zphi(hpxid,z,phi);
}

    /*! Returns the angular coordinates of the center of the pixel with
        number \a pix. */
void MultiResHpx::Pix2ang 
(
 int64 hpxid, 
 int order, 
 pointing& pt
) 
{
  hpxQ.Set(order,hpx_scheme);
  pt = hpxQ.pix2ang(hpxid);
}

	/*! Returns the vector to the center of the pixel with number \a pix. */
void MultiResHpx::Pix2vec 
(
 int64 hpxid, 
 int order, 
 vec3& v
) 
{
  hpxQ.Set(order,hpx_scheme);
  v = hpxQ.pix2vec(hpxid);
}

////#####################################################
////##### PUBLIC INFORMATION METHODS ####################
////#####################################################

    /*! Returns useful information about a given ring of the map.
        \param ring the ring number (the number of the first ring is 1)
		\param order the HEALPix order (level of resolution) for the query
        \param startpix the number of the first pixel in the ring
               (NOTE: this is always given in the RING numbering scheme!)
        \param ringpix the number of pixels in the ring
        \param costheta the cosine of the colatitude of the ring
        \param sintheta the sine of the colatitude of the ring
        \param shifted if \a true, the center of the first pixel is not at
               \a phi=0 */
void MultiResHpx::Get_Ring_Info 
(
 int ring, 
 int order, 
 int64 &startpix, 
 int64 &ringpix,
 double &costheta, 
 double &sintheta, 
 bool &shifted
)
{
  hpxQ.Set(order,hpx_scheme);
  hpxQ.get_ring_info(ring,startpix,ringpix,costheta,sintheta,shifted);
}

	/*! Returns useful information about a given ring of the map.
        \param ring the ring number (the number of the first ring is 1)
		\param order the HEALPix order (level of resolution) for the query
        \param startpix the number of the first pixel in the ring
               (NOTE: this is always given in the RING numbering scheme!)
        \param ringpix the number of pixels in the ring
        \param theta the colatitude (in radians) of the ring
        \param shifted if \a true, the center of the first pixel is not at
               \a phi=0 */
void MultiResHpx::Get_Ring_Info2 
(
 int ring, 
 int order, 
 int64 &startpix, 
 int64 &ringpix,
 double &theta, 
 bool &shifted
)
{
  hpxQ.Set(order,hpx_scheme);
  hpxQ.get_ring_info2(ring,startpix,ringpix,theta,shifted);
}

    /*! Returns useful information about a given ring of the map.
        \param ring the ring number (the number of the first ring is 1)
		\param order the HEALPix order (level of resolution) for the query
        \param startpix the number of the first pixel in the ring
               (NOTE: this is always given in the RING numbering scheme!)
        \param ringpix the number of pixels in the ring
        \param shifted if \a true, the center of the first pixel is not at
               \a phi=0 */
void MultiResHpx::Get_Ring_Info_Small 
(
 int ring, 
 int order, 
 int64 &startpix, 
 int64 &ringpix,
 bool &shifted
)
{
  hpxQ.Set(order,hpx_scheme);
  hpxQ.get_ring_info_small(ring,startpix,ringpix,shifted);
}

    /*! Returns the maximum angular distance (in radian) between any pixel
        center and its corners. 
		\param order the HEALPix order (level of resolution) for the query*/
double MultiResHpx::Max_pixrad(int order)
{
  hpxQ.Set(order,hpx_scheme);
  return hpxQ.max_pixrad();
}

    /*! Returns the maximum angular distance (in radian) between any pixel
        center and its corners in a given ring. 
		\param order the HEALPix order (level of resolution) for the query*/
double MultiResHpx::Max_pixrad(int ring,int order)
{
  hpxQ.Set(order,hpx_scheme);
  return hpxQ.max_pixrad(ring);
}

    /*! Returns a set of points along the boundary of the given pixel.
        \a step=1 gives 4 points on the corners. The first point corresponds
        to the northernmost corner, the subsequent points follow the pixel
        boundary through west, south and east corners.
        \param pix pixel index number
		\param order the HEALPix order (level of resolution) for the query
        \param step the number of returned points is 4*step. */
void MultiResHpx::Boundaries 
(
 int64 hpxid, 
 int order, 
 tsize step, 
 std::vector<vec3> &out
)
{
  pair<int64,int> hpxidx;
  hpxidx.first = hpxid;
  hpxidx.second = order;
  hpxQ.Set(order,hpx_scheme);
  hpxQ.boundaries(hpxidx,step,out);
}


////###############################################
////##### PUBLIC QUERY METHODS ####################
////###############################################

//Inputs: pointing angle = HPX phi(phi)(radians), HPX theta (theta) (radians)
std::vector<MortonNode> MultiResHpx::Search(pointing pt)
{
	// GIS Longitude Range: [-180.0,+180.0] degrees, relative to Prime Meridian.
	// GIS Latitude Range: [-90.0,+90.0] degrees, relative to Equator.
	// HPX Longitude Range: [0.0,360.0] degrees, relative to Prime Meridian and progressing East.
	// HPX Colatitude Range: [0.0,180.0] degrees, relative to North Pole and progressing South to
	// South Pole. Input units must be RADIANS.
	int fn;
	std::vector<MortonNode> found;

	// Check the hpxQ indexing scheme setting.
	// If RING, change to NEST for purpose of query.
	// Then reset hpxQ indexing schme back to RING.
	if( hpx_scheme == RING ) {
		hpxQ.Set(hpxQ.Order(),NEST);
	}

	// Determine which face number (fn) HPX Longitude,Colatitude pair map to.
	fn = hpxQ.FaceNum(pt);

    found = forest_[fn].SearchMortonNodeAtPhiTheta(pt,search_sentinel);

	// Set hpxQ ordering scheme back to whatever it was before.
	if( hpx_scheme == RING ) {
		hpxQ.Set(hpxQ.Order(),RING);
	}

	return found;
}

std::vector<MortonNode> MultiResHpx::Search(int64 hpxid, int order, bool upsearch)
{
	bool Done = false;
	pointing pt;
	int fn;
	std::vector<MortonNode> found;
	std::vector<MortonNode> temp;
	int stop_order;
	Morton m;

	// Get min leaf depth for this tree
	//stop_order = forest_[fn].GetMortonTreeMinDepth();
	//stop_order = MinDepth(); //int(order/3);
	stop_order = 0;

	if(hpxid == -1) {
		return found;
	}
    
	// Check the hpxQ indexing scheme setting.
	// If RING, change to NEST for purpose of query.
	// Then reset hpxQ indexing schme back to RING.
	if( hpx_scheme == RING ) {
		hpxQ.Set(hpxQ.Order(),NEST);
	}

	// First determine the face number where hpxid,order is located this determines
	// which MortonLQT to search.
	fn = hpxQ.FaceNum(hpxid,order);
	
	// Next shift hpxid,order from base N to base 0 range.
	hpxid = hpxid - fn*order_to_npface[order];

	// Now search for query in respective Morton LQT
	found = forest_[fn].SearchMortonNodeHpx(hpxid,order,search_sentinel);

	// Possibility of search hpxid,order pair is of higher resolution than data point
	// in MLQT is mapped to so will search UP the MLQT until either data point is found 
	// or reach specified stop order.
	if( found.size() == 0 && upsearch == true ) {
	
		while( Done == false ) {

			// Compute parent of hpxid,order
			hpxid = hpxid >> 2;
			order -= 1;
			m = HpxToMorton(hpxid,order);

			found = forest_[fn].GetNodeAtMorton(m,-1);

			// If parent hpxid is of order = 0 (root), we're done
			if( found.size() > 0 || order < stop_order ) {

				for(unsigned int i = 0; i < found.size(); i++) {
					// Check if MortonNode has been visited already
					if( found[i].sentinel == search_sentinel ) {
						found.clear();
						Done = true;
					}
					else
					{
						// Update MortonNode as being visited
						if( found[i].data.size() > 0 ) {
							MortonNode update_node;
							update_node = found[i];
							update_node.sentinel = search_sentinel;
							forest_[fn].UpdateMortonNode(update_node);
						}
					}
				}
				Done = true;
			}
		}
	}

	// Set hpxQ ordering scheme back to whatever it was before.
	if( hpx_scheme == RING ) {
		hpxQ.Set(hpxQ.Order(),RING);
	}
	return found;
}

std::vector<MortonNode> MultiResHpx::Search(MortonNode qm, bool upsearch)
{
	bool Done = false;
	pointing pt;
	int fn;
	std::vector<MortonNode> found;
	std::vector<MortonNode> temp;
	int stop_order;
	Morton m;
	int64 hpxid;
	int order;

	// Get min leaf depth for this tree
	stop_order = 0;

	// Check the hpxQ indexing scheme setting.
	// If RING, change to NEST for purpose of query.
	// Then reset hpxQ indexing schme back to RING.
	if (hpx_scheme == RING) {
		hpxQ.Set(hpxQ.Order(), NEST);
	}

	MortonToHpx(qm.m, hpxid, order);

	// Now search for query in respective Morton LQT
	found = forest_[qm.facenum].SearchMortonNodeHpx(hpxid,order, search_sentinel);

	// Possibility of search hpxid,order pair is of higher resolution than data point
	// in MLQT is mapped to so will search UP the MLQT until either data point is found 
	// or reach specified stop order.
	if (found.size() == 0 && upsearch == true) {

		while (Done == false) {

			// Compute parent of hpxid,order
			hpxid = hpxid >> 2;
			order -= 1;
			m = HpxToMorton(hpxid, order);

			found = forest_[fn].GetNodeAtMorton(m, -1);

			// If parent hpxid is of order = 0 (root), we're done
			if (found.size() > 0 || order < stop_order) {

				for (unsigned int i = 0; i < found.size(); i++) {
					// Check if MortonNode has been visited already
					if (found[i].sentinel == search_sentinel) {
						found.clear();
						Done = true;
					}
					else
					{
						// Update MortonNode as being visited
						if (found[i].data.size() > 0) {
							MortonNode update_node;
							update_node = found[i];
							update_node.sentinel = search_sentinel;
							forest_[fn].UpdateMortonNode(update_node);
						}
					}
				}
				Done = true;
			}
		}
	}

	// Set hpxQ ordering scheme back to whatever it was before.
	if (hpx_scheme == RING) {
		hpxQ.Set(hpxQ.Order(), RING);
	}
	return found;
}

//Inputs: center, pointing angle = HPX phi(phi)(radians), HPX theta (theta) (radians)
//        radius, radial distance (radians)
// Using GeographicLib:
// http://geographiclib.sourceforge.net/html/index.html
std::vector<MortonNode> MultiResHpx::QueryDisc(pointing center,double radius)
{
	pointing pt;
	rangeset<int64> covermapHPX;
	std::vector<pair<int64,int>> covermapMRH;
    std::vector<MortonNode> candidate;
	std::vector<int64> covermapHPX2;
	std::vector<MortonNode> found;
	pair<int64,int> nPair;
	int i,j,order=2;
	bool order_set = false;
	double area = radius*radius;

    for( i = 0; i <= 29; i++ ){
		if(order_to_cellres[i] < area && order_set == false){
			order = i-1;
			if(order < 1){ order = 1; }
			order_set = true; i = 30;
			break;
		}
    }
	// Want to cap the maximum HPX disc query resolution to
	// never be greater than maximum resolution of MRH.
	if( order >= MaxDepth() )order = MaxDepth()-2;
	
	// Must make sure order doesn't go less than 1.
	if( order < 1 )order = 1;

	lowHPX.Set(order,hpx_scheme);

	// Update search sentinel value. During search if sentinel value matches discovered
	// the sentinel value in potential found Morton LQT node (data point) we won't return
	// that data point as it's already in the list of found data points.
	search_sentinel += 1;

	// Do low-res hpx disc query to get blanket coverage of data points in and 
	// around query disc which helps limit number of searches
	lowHPX.query_disc_inclusive(center,radius,covermapHPX,1);
	covermapHPX.toVector(covermapHPX2);
	for( i = 0; i < covermapHPX2.size(); i++ ) {
		nPair.first = covermapHPX2[i];
		nPair.second = order;
		covermapMRH.push_back(nPair);
	}

	// Now search for hpxQ results in Morton LQT forest
	found.clear();
	for( i = 0; i < covermapMRH.size(); i++) {
	   candidate.clear();
	   candidate = Search(covermapMRH[i].first,covermapMRH[i].second,true);
	   for( j = 0; j < candidate.size(); j++) {
			// Attempt to filter out false matches by computing distance
			// from query center and compare to query radius.
			// First convert found result GIS long/lat to HPX z/phi
			pt.phi = candidate[j].phi;
			pt.theta = candidate[j].theta;
			if( IsPointInDisc(center,radius,pt) == true ) {
				found.push_back(candidate[j]);
			}
	   }
	}
	return found;
}


// Polygon Query, returns data points that fall within user specified convex polygon.
// The defined polygon must wind counter-clockwise.
// Using GeographicLib:
// http://geographiclib.sourceforge.net/html/index.html
std::vector<MortonNode> MultiResHpx::QueryPolygon
( 
 std::vector<pointing> poly
)  
{
	pointing pt,pt0,p0proj;
	std::vector<pointing> p;
	std::vector<pair<int64,int>> covermapMRH;
    std::vector<MortonNode> candidate,found;
	pair<int64,int> nPair;
	std::vector<int64> covermapHPX2;
	rangeset<int64> covermapHPX;
	int i,j,k,order=2;
	bool over_horizon = false,order_set = false;
	double x,y,azi,rk,min_phi,max_phi,min_theta,max_theta,poly_area;
	GeographicLib::Gnomonic* g = new GeographicLib::Gnomonic(GeographicLib::Geodesic(1.0,0.0));
	min_phi = 9999.0; min_theta = 9999.0; max_phi = -9999.0; max_theta = -9999.0;
	
	// Compute extents of query polygon in order to determine the lowest resolution HPX cell size
	// that will produce a HPX coverage map with as few cells as possible while at same time not
	// generating too many unnecessary MRH searches.
	for(unsigned int n = 0; n < poly.size(); n++) {
		if( poly[n].phi < min_phi ) { min_phi = poly[n].phi; }
		if( poly[n].phi > max_phi ) { max_phi = poly[n].phi; }
		if( poly[n].theta < min_theta ) { min_theta = poly[n].theta; }
		if( poly[n].theta > max_theta ) { max_theta = poly[n].theta; }
	}
	poly_area = RadialDist(pointing(min_theta,min_phi),pointing(max_theta,max_phi));
	poly_area *= poly_area;
    for( k = 0; k <= 29; k++ ) 
    {
		if(order_to_cellres[k] < poly_area && order_set == false)
		{
			order = k-1;
			if(order < 1){ order = 1; }
			order_set = true; k = 30;
			break;
		}
    }
	// Want to cap the maximum HPX disc query resolution to
	// never be greater than maximum resolution of MRH.
	if( order >= MaxDepth() ) {
		order = MaxDepth()-2;
	}
	// Must make sure order isn't less than 1 either
	if( order < 1 ) {
		order = 1;
	}
	lowHPX.Set(order,hpx_scheme);

	// Update search sentinel value. During search if sentinel value matches discovered
	// the sentinel value in potential found Morton LQT node (data point) we won't return
	// that data point as it's already in the list of found data points.
	search_sentinel += 1;

	// Do low-res hpx polygon query to get blanket coverage of data points in and 
	// around query polygon which helps limit number of searches
	lowHPX.query_polygon_inclusive(poly,covermapHPX,1);
	covermapHPX.toVector(covermapHPX2);
	covermapMRH.clear();
	for( i = 0; i < covermapHPX2.size(); i++ ) {
		nPair.first = covermapHPX2[i];
		nPair.second = order;
		covermapMRH.push_back(nPair);
	}

	// Now Search through blanket coverage to discover data points within
	// query polygon.
	found.clear();
	for( i = 0; i < covermapMRH.size(); i++) {

		candidate = Search(covermapMRH[i].first,covermapMRH[i].second,true);

		for( j = 0; j < candidate.size(); j++) {
			//pt.phi = temp[j].phi;
			//pt.theta = temp[j].theta;

			// Convert point to GIS degrees longitude,latitude
			pt0 = HPXtoGIS(pointing(candidate[j].theta,candidate[j].phi));
			pt0.phi *= rad2degr;
			pt0.theta *= rad2degr;

			// Do ellipsoidal gnomonic projection of found data point onto 2D flat surface
			g->Forward(pt0.theta,pt0.phi,pt0.theta,pt0.phi,x,y,azi,rk);
			p0proj.theta = x;
			p0proj.phi = y;

			// Use GeographicLib to do ellipsoidal gnomonic projection of query
			// polygon so can then use standard point-in-polygon test. Basically
			// project each polygon point given query point as the origin.
			p.clear();
			over_horizon = false;
			for( k = 0; k < poly.size(); k++ ) {
				pt = HPXtoGIS(poly[k]);
				pt.phi *= rad2degr;
				pt.theta *= rad2degr;
				g->Forward(pt0.theta,pt0.phi,pt.theta,pt.phi,x,y,azi,rk);
				p.push_back(pointing(y,x));

				// Flag over the horizon point
				if(rk < 0 ) {
					over_horizon = true;
					k = poly.size(); // exit the loop
				}
			}

			// If found point is flagged as over the horizon then reject, else
			// apply point-in-polygon test to filter out data points outside
			// polygon query
			if( over_horizon == false ) {
				if( IsPointInPoly2(p0proj,p) ) {
					found.push_back(candidate[j]);
				}
			}
		}
		candidate.clear();
	}
	return found;
}


std::vector<MortonNode> MultiResHpx::QueryStrip 
( 
 double theta1, 
 double theta2 
)
{
	pointing pt;
    std::vector<MortonNode> found;
	MortonNode m;
	int i,j=0;
	found.clear();

	// Do linear search through all Twelve MortonLQTs and pick out
	// all data points whose theta values fall within the strip query
	// range.
	for( i = 0; i < 12; i++ ) {
		for( j = 0; j < forest_[i].GetNumMortonNodes(); j++) {
			m = forest_[i].GetNodeAtIndex(j);
			if( m.data.size() > 0 ) {
				pt.phi = m.phi;
				pt.theta = m.theta;
				if( IsPointInStrip(theta1,theta2,pt) ) {
					found.push_back(m);
				}				
			}
		}
	}

	return found;
}




std::vector<MortonNode> MultiResHpx::Neighbors( pointing pt, int64 order )
{
	pointing pt0,p0proj;
	std::vector<pointing> p;
	bool over_horizon = false;
	bool SKIP = false;
	GeographicLib::Gnomonic* g = new GeographicLib::Gnomonic(GeographicLib::Geodesic(1.0,0.0));
	double x,y,azi,rk;
	int i,j,k=0;
	pair<int64,int> hpxid,Qhpxid;
	pointing newpt;
	std::vector<vec3> out;
	fix_arr<pair<int64,int>,8> covermapHPX;
	std::vector<pointing> poly;
    std::vector<MortonNode> QNode,found,candidate;
	int64 save_order;
	pointing min_phi,max_phi,min_theta,max_theta;

	min_phi.phi = 999999.0; min_phi.theta = 0.0;
	max_phi.phi = -999999.0; max_phi.theta = 0.0;
	min_theta.theta = 999999.0; min_theta.phi = 0.0;
	max_theta.theta = -999999.0; max_theta.phi = 0.0;

	// Update search sentinel value. During search if sentinel value matches discovered
	// the sentinel value in potential found Morton LQT node (data point) we won't return
	// that data point as it's already in the list of found data points.
	search_sentinel += 1;

	// Check the hpxQ indexing scheme setting.
	// If RING, change to NEST for purpose of query.
	// Then reset hpxQ indexing schme back to RING.
	if( hpx_scheme == RING ) {
		hpxQ.Set(hpxQ.Order(),NEST);
	}

	// First search MRH for query point. If exists continue
	// neighbor query at user specified level of resolution (order).
    // Otherwise return empty vector. Save original hpxQ order.
	QNode = Search(pt);
	if( QNode.size() > 0 ) {
		save_order = hpxQ.Order();
		hpxQ.Set(order,NEST);
	} else {
			return found;
	}

	// Next compute query point's hpxid for neighbor query
	Qhpxid.first = hpxQ.ang2pix(pt);
    Qhpxid.second = order;

	//Call healpix_custom's neighbors
    hpxQ.neighbors( Qhpxid, covermapHPX);

	// Compute the neighbor net using the boundaries of
    // neighboring cells of query cell
	for( i = 0; i < covermapHPX.size(); i++ ) {
		hpxid = covermapHPX[i];
		hpxQ.boundaries(hpxid,1,out);

		// Convert vec3 to phi/theta pairs
		for(j = out.size()-1; j >= 0; j--) {
			Vec2ang(out[j],newpt);

			//Track min and max phi theta pairs to determine the
			//boundary of the neighbor net

			//New min theta?
			if( newpt.theta < min_theta.theta ) {
				min_theta.theta = newpt.theta;  min_theta.phi = newpt.phi;
			}
			//New max theta?
			if( newpt.theta > max_theta.theta ) {
				max_theta.theta = newpt.theta;  max_theta.phi = newpt.phi;
			}
			//New min phi?
			if( newpt.phi < min_phi.phi ) {
				min_phi.theta = newpt.theta;  min_phi.phi = newpt.phi;
			}
			//New max phi?
			if( newpt.phi > max_phi.phi ) {
				max_phi.theta = newpt.theta;  max_phi.phi = newpt.phi;
			}
		}
	}

	// Will use neighbor net polygon definition to call 
	// define point-in-polygon test to filtering out
	// any discovered neighboring data points.
	poly.clear();
	poly.push_back(min_theta);  poly.push_back(max_phi);
	poly.push_back(max_theta);	poly.push_back(min_phi);  

	// Now Search through blanket coverage to discover data points within
	// query polygon.
	found.clear();
	for( i = 0; i < covermapHPX.size(); i++) {
        
		// Search for each HPXid in MRH
		candidate.clear();
		hpxid = covermapHPX[i];
		candidate = Search(hpxid.first,hpxid.second,true);

		for( j = 0; j < candidate.size(); j++) {
			SKIP = false;

			// Skip if candidate is query point
			if( approx(candidate[j].phi,0.0) && approx(candidate[j].theta,0.0) ) {
				SKIP = true;
			}

			if( approx(QNode[0].phi,candidate[j].phi) && approx(QNode[0].theta,candidate[j].theta) ) {
				SKIP = true;      
			}

			// Check if candidate already in found list
			for(k = 0; k < found.size(); k++) {
				if( approx(found[k].phi,candidate[j].phi) &&
					approx(found[k].theta,candidate[j].theta) )
				{
					SKIP = true;
				}
			}

			if( SKIP == false ) {
				// Convert point to GIS degrees longitude,latitude
				pt0 = HPXtoGIS(pointing(candidate[j].theta,candidate[j].phi));
				pt0.phi *= rad2degr;
				pt0.theta *= rad2degr;

				// Do ellipsoidal gnomonic projection of found data point onto 2D flat surface
				g->Forward(pt0.theta,pt0.phi,pt0.theta,pt0.phi,x,y,azi,rk);
				p0proj.theta = x;
				p0proj.phi = y;

				// Use GeographicLib to do ellipsoidal gnomonic projection of query
				// polygon so can then use standard point-in-polygon test. Basically
				// project each polygon point given query point as the origin.
				p.clear();
				over_horizon = false;
				for( k = 0; k < poly.size(); k++ ) {
					pt = HPXtoGIS(poly[k]);
					pt.phi *= rad2degr;
					pt.theta *= rad2degr;
					g->Forward(pt0.theta,pt0.phi,pt.theta,pt.phi,x,y,azi,rk);
					p.push_back(pointing(y,x));

					// Flag over the horizon point
					if(rk < 0 ) {
						over_horizon = true;
						k = poly.size(); // exit the loop
					}
				}

				// If found point is flagged as over the horizon then reject, else
				// apply point-in-polygon test to filter out data points outside
				// polygon query
				if( over_horizon == false ) {
					if( IsPointInPoly2(p0proj,p) ) {
						found.push_back(candidate[j]);
					}
				}
			}
		}

	}

	// Set hpxQ order and scheme back to whatever it was before.
	hpxQ.Set(save_order,hpx_scheme);
	if( hpx_scheme == RING ) {
		hpxQ.Set(hpxQ.Order(),RING);
	}
	return found;
}


std::vector<MortonNode> MultiResHpx::NearNeighbors(MortonNode m)
{
	int i,j = 0;
	pair<int64, int> hpxid, Qhpxid;
	fix_arr<pair<int64, int>, 8> covermapHPX;
	std::vector<pointing> poly;
	std::vector<MortonNode> found,temp;
	int64 save_order;
	MortonNode nextMN;
	Morton nextM;

	// Update search sentinel value. During search if sentinel value matches discovered
	// the sentinel value in potential found Morton LQT node (data point) we won't return
	// that data point as it's already in the list of found data points.
	search_sentinel += 1;

	// Check the hpxQ indexing scheme setting.
	// If RING, change to NEST for purpose of query.
	// Then reset hpxQ indexing schme back to RING.
	if (hpx_scheme == RING) {
		hpxQ.Set(hpxQ.Order(), NEST);
	}

	// First verify the query MortonNode is valid. If so continue
	// neighbor query at user specified level of resolution (order).
	// Otherwise return empty vector. Save original hpxQ order.
	if (m.sentinel == -1 && m.data[0] == -1)
	{
		return found;
	}

	// First search MRH for query point. If exists continue
	// neighbor query at user specified level of resolution (order).
	// Otherwise return empty vector. Save original hpxQ order.
	save_order = hpxQ.Order();
	hpxQ.Set(m.m.LEVEL, NEST);

	// Next compute query point's hpxid for neighbor query
	MortonToHpx(m.m, Qhpxid.first, Qhpxid.second);

	//Call healpix_custom's neighbors
	hpxQ.neighbors(Qhpxid, covermapHPX);

	// Now process the HPX covermap to convert each hpx id/order to Morton
	for (i = 0; i < covermapHPX.size(); i++) 
	{
		hpxid = covermapHPX[i];

		// Convert HPXID/ORDER to equivalent Morton number
		nextM = HpxToMorton(hpxid.first, hpxid.second);

		// Next retrieve MortonNode (if it exists)
		if( GetMortonNodeAtMorton(nextM, nextMN) )
		{
			// If this MortonNode has children we need to search them!
			if (nextMN.childrenYN == 1)
			{
				// Now do a Search for this MortonNode to possibly uncover populated
				// child nodes of this node (if it has children)
				temp.clear();
				temp = Search(nextMN, false);

				for (j = 0; j < temp.size(); j++)
				{
					found.push_back(temp[j]);
				}
			}
			else
			{
				found.push_back(nextMN);
			}
		}
	}
	// Set hpxQ order and scheme back to whatever it was before.
	hpxQ.Set(save_order, hpx_scheme);
	if (hpx_scheme == RING) {
		hpxQ.Set(hpxQ.Order(), RING);
	}

	return found;
}

bool MultiResHpx::IsPointInCoverageMap(pointing pt,std::vector<pair<int64,int>> pixset)
{
	int64 hpxid_pt;
	//See if pt maps to ANY of the pixset coverage map cells
	for(unsigned int i = 0; i < pixset.size(); i++) {
		Ang2pix(pt,pixset[i].second,hpxid_pt);
		if(hpxid_pt == pixset[i].first) {
			return true;
		}
	}
	return false;
}


std::vector<std::pair<MortonNode,MortonNode>> MultiResHpx::TwoPointCorrBin( double radius )
{
	return _twopointcorrbin_internal(radius);
}


// Two-Point Correlation Binning Algorithm
// Essentially need to discover all unique pairs of points that lie
// distance "radius" between each other and return this list of pairs.
// Brute force method will do QueryDisc search from EACH data point
// in the data structure and record each unique pair discovered, where
// a pair is the query point and data point found within distance "radius"
// from the query point. If pair is unique it's added to the return list.
std::vector<std::pair<MortonNode,MortonNode>> MultiResHpx::_twopointcorrbin_internal(double radius )
{
	std::vector<std::pair<MortonNode,MortonNode>> pairs;
	std::vector<MortonNode> found;
	pair<MortonNode,MortonNode> foundPair;
	found.clear();
	MortonNode nextNode;
	pointing pt;
	unsigned int i, j, k;
	pairs.clear();

	// Loop through each MortonLQT
	for( i = 0; i < 12; i++) {

		// Loop through each MortonLQT Node, checking for data points
		for( j = 0; j < forest_[i].GetNumMortonNodes(); j++) {
			nextNode = forest_[i].GetNodeAtIndex(j);

			// Check for data point
			if( nextNode.data.size() > 0 ) {

				// Do QueryDisc from this data point
				pt.phi = nextNode.phi;
				pt.theta = nextNode.theta;

				found = QueryDisc(pt,radius);

				// Create found pairs
				for( k = 0; k < found.size(); k++) {

					// Check for data point
					if (found[k].data.size() > 0)
					{
						// Skip the query point matching with itself in found
						// disc query points.
						if (Equals(nextNode.m, found[k].m))
						{
							continue;
						}
						else
						{
							foundPair.first = nextNode;
							foundPair.second = found[k];

							// Check that foundPair is unique!
							if (IsMortonNodePairInList(pairs, foundPair))
							{
								continue;
							}
							else
							{
								pairs.push_back(foundPair);
							}
						}
					}
				}
				found.clear();
			}
		}
	}
	return pairs;
}

// Given a data index from MultiResHpx_Map, find and return the stored MortonNode
bool MultiResHpx::GetMortonNodeAtDataIndex(int64 data_idx, MortonNode& foundM)
{
	// Do linear search through all Twelve MortonLQTs to find and return the 
	// MortonNode that contains data_idx.
	for (unsigned int i = 0; i < 12; i++) {
		for (unsigned int j = 0; j < forest_[i].GetNumMortonNodes(); j++) {
			foundM = forest_[i].GetNodeAtIndex(j);
			for (unsigned int k = 0; k < foundM.data.size(); k++)
			{
				if (foundM.data[k] == data_idx)
				{
					return true;
				}
			}
		}
	}
	// If made it here, then MortonNode with data_idx wasn't found. Return false.
	return false;
}

// Given a Morton number, find the MortonNode (if it exists)
bool MultiResHpx::GetMortonNodeAtMorton(Morton qm, MortonNode& foundM)
{
	// Do linear search through all Twelve MortonLQTs to find and return the 
	// MortonNode who's address is the given Morton number
	for (unsigned int i = 0; i < 12; i++) {
		for (unsigned int j = 0; j < forest_[i].GetNumMortonNodes(); j++) {
			foundM = forest_[i].GetNodeAtIndex(j);
			if (Equals(qm, foundM.m))
			{
				return true;
			}
		}
	}
	// If made it here, then MortonNode with data_idx wasn't found.
	// Return false

	return false;
}

vector< MortonLQT > MultiResHpx::GetForest()
{
	return forest_;
}

int MultiResHpx::SaveToFile(std::string filename)
{
	ofstream fp;
	fp.open(filename.c_str());

	if(fp.fail())
	{
	  cout << "Unable to open MRH output file: " << filename.c_str() << endl;
	  return false;     
	}	

	//// Write out header & data, base cell by base cell
	for(unsigned int i = 0; i < 12; i++) {
		forest_[i].SaveTreeToFile(fp);	 
	}

	fp.close();  	

	return true;
}


int MultiResHpx::LoadFromFile(std::string filename)
{
	ifstream fp;
	fp.open(filename.c_str());

	if(fp.fail())
	{
	  cout << "Unable to open MRH input file: " << filename.c_str() << endl;
	  return false;     
	}

	//// Read in MRH input file, base cell by base cell
	for(unsigned int i = 0; i < 12; i++) {
		forest_[i].LoadTreeFromFile(fp);	 
	}

	fp.close();  	


	return true;

}