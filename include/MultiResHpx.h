/*
 *  MultiResHpx
 *
 *  Author: Robert Youngren, robert.youngren@gmail.com
 *
 */

#ifndef MULTIRESHEALPIX_H
#define MULTIRESHEALPIX_H

#define _USE_MATH_DEFINES
#define _NUMFACES 12
#define MAX_HPX_ORDER32 13
#define MAX_HPX_ORDER64 29
#define MIN_ORDER 1
#define MIN_NSIDE 2
#define MAX_NSIDE_32 8192
#define MAX_NSIDE_64 536870912

//#define CRITCOUNT
//#define UPSEARCHCOUNT

//#define MRH_VERBOSE 

//#define HPX_COVERAGE_MAP
//#define HPX_COVERAGE_MAP_CELL_COUNT

#include "math.h"


//HEALPix References
#include "healpix_custom.h"
#include "healpix_base.h"
#include "geom_utils.h"
//#include "MortonLQT.h"
#include "MortonLQT64.h"
#include "MultiResHpx_tables.h"


//GeographicLib References
#include <GeographicLib/Gnomonic.hpp>

using namespace std;


pointing HPXtoGIS(pointing pt);
pointing GIStoHPX(pointing pt);
pointing RADECtoHPX(pointing pt);
pointing HPXtoRADEC(pointing pt);
bool my_equals(double a,double b);
bool pointing_equals(pointing a, pointing b);
bool IsPointInPoly(pointing pt,std::vector<pointing> poly); 
bool IsPointInPoly2(pointing pt,std::vector<pointing> poly);
bool IsPointInPoly3(pointing pt,std::vector<pointing> poly);
bool IsPointInStrip(double theta1, double theta2, pointing pt);
bool IsPointInDisc(pointing center,double radius,pointing pt);
double ComputeArcLengthBetween(pointing pt1, pointing pt2);
void HpxPZCellBoundaries(pair<int64,int> hpxid,int step);
double RadialDist(pointing pt1,pointing pt2);
void Vec2ang(vec3 vec,pointing &ang);
void Ang2vec(pointing &ang, vec3 vec);
bool SortFunctionPixset( pair<int64,int> a, pair<int64,int> b );

class MultiResHpx 
{
public:
	MultiResHpx(){}
	MultiResHpx( int Max_Depth, Healpix_Ordering_Scheme scheme );
	~MultiResHpx();

	// Accessors
	const int& emptyValue( int fn ) const;


	// Mutators

	// Utilities
	void NormGIStoHPX(double phi, double theta,pointing &ptg,int& fn);
	void DenormHPXtoGIS(double phi, double theta,pointing &ptg,int& fn);
	//void ShiftAngleGIStoHPX(pointing &ptg,double dir);
	//void ShiftAngleBase0toBaseN(int fn, pointing &ptg,double dir);
	void PrintTreeAtIndex(int fn);
	void PrintListAtIndex(int fn);
	void Delete(pointing pt);
	int  NodeCountTree(int fn);
	int SaveToFile(std::string filename);
	int LoadFromFile(std::string filename);
 

	// Indexing operators
	//std::vector<Point2D> operator() ( pointing pt );
	//std::vector<Point2D> operator() ( std::vector<int> mrh );
	//std::vector<Point2D> At( pointing pt );
	//std::vector<Point2D> At( std::vector<int> mrh )const;

	//// Setters
	void Insert (pointing pt, int64 mapidx);
	void Insert2 (pointing pt, int64 mapidx);
	void SetMaxDepth(int maxDepth);

	// Build the MortonLQT Forest
	void BuildForest();

	void BuildForestFromArchive(ifstream& fp);

	//// Getters
	bool GetMortonNodeAtDataIndex(int64 qidx, MortonNode& foundM);
	bool GetMortonNodeAtMorton(Morton qm, MortonNode& foundM);

	vector< MortonLQT > GetForest();

	//int getBytesAtFace(int fn) { return forest_[fn].bytes(); }
	//inline int size();

    //// Translators & Conversions 

	//*! Returns the number of the ring in which \a pix,order lies. */
    void Pix2ring (int64 hpxid, int order,int64& ring);

    void Xyf2pix(int ix, int iy, int fn, int order, int64& hpxid);
	
    void Pix2xyf(int64 hpxid,int order, int &ix, int &iy, int &fn);

    /*! Translates a pixel number from NEST to RING. */
    void Nest2ring (int64 hpxid, int order, int64& Rhpxid);
   
	//*! Translates a pixel number from RING to NEST. */
    void Ring2nest (int64 Rhpxid, int order, int64& Nhpxid);

	//*! Translates a pixel number from NEST to its Peano index. */
    void Nest2peano (int64 hpxid, int order, int64& pid);

	//*! Translates a pixel number from its Peano index to NEST. */
    void Peano2nest (int64 pid, int order, int64& Nhpxid);

    /*! Returns the number of the pixel which contains the angular coordinates
        (\a z:=cos(theta), \a phi).
        \note This method is inaccurate near the poles at high resolutions. */
    void Zphi2pix (double z, double phi, int order, int64& hpxid);

    /*! Returns the number of the pixel which contains the angular coordinates
        \a ang. */
    void Ang2pix (pointing &ang, int order, int64& hpxid);

    /*! Returns the number of the pixel which contains the vector \a vec
        (\a vec is normalized if necessary). */
    void Vec2pix (vec3 &v,int order, int64& hpxid);

    /*! Returns the angular coordinates (\a z:=cos(theta), \a phi) of the center
        of the pixel with number \a pix.
        \note This method is inaccurate near the poles at high resolutions. */
    void Pix2zphi (int64 hpxid, int order, double &z, double &phi);

    /*! Returns the angular coordinates of the center of the pixel with
        number \a pix. */
    void Pix2ang (int64 hpxid, int order, pointing& pt);

	/*! Returns the vector to the center of the pixel with number \a pix. */
    void Pix2vec (int64 hpxid, int order, vec3& v);

	//// Information

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
    void Get_Ring_Info (int ring, int order, int64 &startpix, int64 &ringpix,
      double &costheta, double &sintheta, bool &shifted);

	/*! Returns useful information about a given ring of the map.
        \param ring the ring number (the number of the first ring is 1)
		\param order the HEALPix order (level of resolution) for the query
        \param startpix the number of the first pixel in the ring
               (NOTE: this is always given in the RING numbering scheme!)
        \param ringpix the number of pixels in the ring
        \param theta the colatitude (in radians) of the ring
        \param shifted if \a true, the center of the first pixel is not at
               \a phi=0 */
	void Get_Ring_Info2 (int ring, int order, int64 &startpix, int64 &ringpix,
      double &theta, bool &shifted);

    /*! Returns useful information about a given ring of the map.
        \param ring the ring number (the number of the first ring is 1)
		\param order the HEALPix order (level of resolution) for the query
        \param startpix the number of the first pixel in the ring
               (NOTE: this is always given in the RING numbering scheme!)
        \param ringpix the number of pixels in the ring
        \param shifted if \a true, the center of the first pixel is not at
               \a phi=0 */
    void Get_Ring_Info_Small (int ring, int order, int64 &startpix, int64 &ringpix,
        bool &shifted);


    /*! Returns the maximum angular distance (in radian) between any pixel
        center and its corners. 
		\param order the HEALPix order (level of resolution) for the query*/
    double Max_pixrad(int order);

    /*! Returns the maximum angular distance (in radian) between any pixel
        center and its corners in a given ring. 
		\param order the HEALPix order (level of resolution) for the query*/
    double Max_pixrad(int ring,int order);

    /*! Returns a set of points along the boundary of the given pixel.
        \a step=1 gives 4 points on the corners. The first point corresponds
        to the northernmost corner, the subsequent points follow the pixel
        boundary through west, south and east corners.
        \param pix pixel index number
		\param order the HEALPix order (level of resolution) for the query
        \param step the number of returned points is 4*step. */
    void Boundaries (int64 hpxid, int order, tsize step, std::vector<vec3> &out);

	//// Queries
	std::vector<MortonNode> Search (pointing pt);
	std::vector<MortonNode> Search (int64 hpxid, int order, bool upsearch);
	std::vector<MortonNode> Search (MortonNode qm, bool upsearch);

    std::vector<MortonNode> QueryDisc (pointing pt, double radius);

	std::vector<MortonNode> QueryPolygon ( std::vector<pointing> poly);

    std::vector<MortonNode> QueryStrip ( double theta1, double theta2);

	std::vector<MortonNode> Neighbors( pointing pt, int64 order );

	std::vector<MortonNode> NearNeighbors(MortonNode m);

	std::vector<std::pair<MortonNode,MortonNode>> TwoPointCorrBin( double radius );

	bool IsPointInCoverageMap(pointing pt,std::vector<pair<int64,int>> pixset);


	//// Data dump


    /*! Returns a constant reference to the healpix_custom query handle. */
    const Healpix_Custom &HPX_Query_Handle() const { return hpxQ; }

	int64 NumNodes();
	int64 NumRec();
	int64 MemSizeKb();
	int MaxDepth();
	int MinDepth();
	double AvgDepth();
	int64 NumNodesAtMLQ(int);
	int DepthAtMLQ(int);
	double AvgDepthAtMLQ(int);
	int GetTreeDepthLimit();
	int GetCritCount();
	int GetCoverMapSize();
	int GetCoverMapCellRes();
	void ResetCritCount();
	std::vector<int> GetNodeSizeHistogram();


	int GetUpSearchCount();
	void ResetUpSearchCount();

protected:

	//### Internal Queries 
	std::vector<std::pair<MortonNode,MortonNode>> _twopointcorrbin_internal(double radius );


private:
    Healpix_Ordering_Scheme hpx_scheme;
	int max_tree_depth;
	int64 num_rec;
	vector< MortonLQT > forest_;
	Healpix_Custom hpxQ;
	Healpix_Base lowHPX;
	int search_sentinel; 
	std::vector<double> HpxCellRadiusTable;
	int MRH_CRIT_COUNT;
	int MRH_COVER_MAP_SIZE;
	int MRH_COVER_MAP_CELL_RES;
	int MRH_UPSEARCH_COUNT;
    int UPSEARCH_COUNT;
};


inline void MultiResHpx::PrintTreeAtIndex(int fn)
{
	forest_[fn].PrintMortonTree();   
}

inline void MultiResHpx::PrintListAtIndex(int fn)
{
	forest_[fn].PrintMortonList();
}

inline int MultiResHpx::NodeCountTree(int fn)
{
	return forest_[fn].GetNumMortonNodes();
}

inline int MultiResHpx::GetTreeDepthLimit()
{
	return max_tree_depth;
}

inline 	void MultiResHpx::SetMaxDepth(int maxDepth)
{
	max_tree_depth = maxDepth;
}



#endif
