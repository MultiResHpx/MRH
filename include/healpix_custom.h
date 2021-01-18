#ifndef HEALPIX_CUSTOM_H
#define HEALPIX_CUSTOM_H

#include "healpix_base.h"
#include "array.h"

#define DEBUG false

using namespace std;

class Healpix_Custom: public Healpix_Base
{

public:
	//*** CONSTRUCTORS ***//
	/*! Constructs an unallocated object. */
	Healpix_Custom () {};
	
	/*! Constructs an object with a given \a order and the ordering
	scheme \a scheme. */
	Healpix_Custom (int order, Healpix_Ordering_Scheme scheme)
	{ Set (order, scheme); }

	/*! Constructs an object with a given \a nside and the ordering
	scheme \a scheme. The \a nside_dummy parameter must be set to
	SET_NSIDE. */
	Healpix_Custom (int nside, Healpix_Ordering_Scheme scheme, const nside_dummy)
	{ SetNside (nside, scheme); }

   /*! Given Pointing Angle, returns the face number of query pixel. */
   const int FaceNum( const pointing ptg )const;

   const int FaceNum( const int64 hpxid, const int order )const;

   void boundaries(pair<int64,int> pix, tsize step,vector<vec3> &out)const;

   /*! Range Queries */

    /*! Returns the range set of all pixels whose centers lie within the disk
        defined by \a dir and \a radius.
        \param dir the angular coordinates of the disk center
        \param radius the radius (in radians) of the disk
        \param pixset a \a pair object containing the index and order of the lowest
		       level resolution pixel whose center lies inside the disk
        \note This method is more efficient in the RING scheme. */
    void query_disc (	pointing ptg, 
						double radius, 
						vector< pair<int64,int> > &pixset) const;
    
	/*! Returns the range set of all pixels which overlap with the disk
        defined by \a dir and \a radius.
        \param dir the angular coordinates of the disk center
        \param radius the radius (in radians) of the disk
        \param pixset a \a pair object containing the index and order of the lowest
		       level resolution pixel that overlaps with the disk.
        \param fact The overlapping test will be done at the resolution
           \a fact*nside. For NESTED ordering, \a fact must be a power of 2,
           else it can be any positive integer. A typical choice would be 4.
        \note This method may return some pixels which don't overlap with
           the disk at all. The higher \a fact is chosen, the fewer false
           positives are returned, at the cost of increased run time.
        \note This method is more efficient in the RING scheme. */
    void query_disc_inclusive (	pointing ptg, 
								double radius, 
								vector< pair<int64,int> > &pixset,
							    int fact=1) const;

    void query_multidisc ( const arr<vec3> &norm,
                           const arr<double> &rad, 
						   int fact, 
						   vector< pair<int64,int> > &pixset) const;

    void query_multidisc_general ( const arr<vec3> &norm, 
		                           const arr<double> &rad,
                                   bool inclusive, 
								   const std::vector<int> &cmds, 
								   vector< pair<int64,int> > &pixset) const;


    /*! Returns a range set of pixels whose centers lie within the convex
        polygon defined by the \a vertex array.
        \param vertex array containing the vertices of the polygon.
        \param pixset a \a pair object containing the index and order of the lowest
		       level resolution pixel whose centers lie inside the polygon
        \note This method is more efficient in the RING scheme. */
    void query_polygon ( const std::vector<pointing> &vertex,
                         vector< pair<int64,int> > &pixset) const;

    /*! Returns a range set of pixels which overlap with the convex
        polygon defined by the \a vertex array.
        \param vertex array containing the vertices of the polygon.
        \param pixset a \a pair object containing the index and order of the lowest
		       level resolution pixel overlapping with the polygon.
        \param fact The overlapping test will be done at the resolution
           \a fact*nside. For NESTED ordering, \a fact must be a power of 2,
           else it can be any positive integer. A typical choice would be 4.
        \note This method may return some pixels which don't overlap with
           the polygon at all. The higher \a fact is chosen, the fewer false
           positives are returned, at the cost of increased run time.
        \note This method is more efficient in the RING scheme. */
    void query_polygon_inclusive ( const std::vector<pointing> &vertex,
                                   vector< pair<int64,int> > &pixset, 
								   int fact=1) const;

    /*! Returns a \a pair object containing the index and order of the lowest
		       level resolution pixel whose centers lie within the colatitude
        range defined by \a theta1 and \a theta2 (if \a inclusive==false), or
        which overlap with this region (if \a inclusive==true). If
        \a theta1<theta2, the region between both angles is considered,
        otherwise the regions \a 0<theta<theta2 and \a theta1<theta<pi.
        \param theta1 first colatitude
        \param theta2 second colatitude
        \param inclusive if \a false, return the exact set of pixels whose
           pixels centers lie within the region; if \a true, return all pixels
           that overlap with the region. */
    void query_strip ( double theta1, 
		               double theta2, 
					   bool inclusive,
                       vector< pair<int64,int> > &pixset) const;

    /*! Returns the neighboring pixels of \a pix in \a result.
        On exit, \a result contains (in this order)
        the pixel pairs (index,order) of the SW, W, NW, N, NE, E, SE and S neighbor
        of \a pix. If a neighbor does not exist (this can only be the case
        for the W, N, E and S neighbors), its entry is set to -1.

        \note This method works in both RING and NEST schemes, but is
          considerably faster in the NEST scheme. */
	void neighbors( pair<int64,int> pix, fix_arr<pair<int64,int>,8> &result) const;

protected:


	/*! Range Query Helpers */
   void query_disc_internal (	pointing ptg, 
								double radius,
								int fact, 
								vector< pair<int64,int> > &pixset) const;

   void query_strip_internal ( double theta1, 
	                           double theta2, 
							   bool inclusive,
                               vector< pair<int64,int> > &pixset) const;

   void query_polygon_internal ( const std::vector<pointing> &vertex, 
	                             int fact,
                                 vector< pair<int64,int> > &pixset) const;


};

#endif