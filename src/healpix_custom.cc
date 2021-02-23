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

#include "healpix_custom.h"
#include "lsconstants.h"
#include "geom_utils.h"

using namespace std;

namespace {

vec3 locToVec3 (double z, double phi, double sth, bool have_sth)
  {
  if (have_sth)
    return vec3(sth*cos(phi),sth*sin(phi),z);
  else
    {
    vec3 res;
    res.set_z_phi (z, phi);
    return res;
    }
  }

} // unnamed namespace


//#########################################################
//#                   UTILITY METHODS                     #
//#########################################################

void Healpix_Custom::boundaries
(
 pair<int64,int> pix, 
 tsize step,
 vector<vec3> &out
) const
{
  out.resize(4*step);
  int ix, iy, face;
  pix2xyf(pix.first, ix, iy, face);
  double dc = 0.5 / nside_;
  double xc = (ix + 0.5)/nside_, yc = (iy + 0.5)/nside_;
  double d = 1.0/(step*nside_);
  for (tsize i=0; i<step; ++i){
    double z, phi, sth;
    bool have_sth;
    xyf2loc(xc+dc-i*d, yc+dc, face, z, phi, sth, have_sth);
	out[i] = locToVec3(z, phi, sth, have_sth);
    xyf2loc(xc-dc, yc+dc-i*d, face, z, phi, sth, have_sth);
    out[i+step] = locToVec3(z, phi, sth, have_sth);
    xyf2loc(xc-dc+i*d, yc-dc, face, z, phi, sth, have_sth);
    out[i+2*step] = locToVec3(z, phi, sth, have_sth);
    xyf2loc(xc+dc, yc-dc+i*d, face, z, phi, sth, have_sth);
    out[i+3*step] = locToVec3(z, phi, sth, have_sth);
  }
}


const int Healpix_Custom::FaceNum
( 
 const pointing ptg 
)const
{
  return( ang2pix(ptg)>>(2*order_) );
}

const int Healpix_Custom::FaceNum
(
const int64 hpxid,
const int order
)const
{
	return hpxid >> (2*order);
}

inline void append_pixel_range
(
 int64 start, 
 int64 end,
 int order, 
 vector<pair<int64,int>>& container)
{
	int64 i = 0;
	for(i=start;i<end;i++) {
		container.push_back(make_pair(i,order));
	}
}

/* Short note on the "zone":
   zone = 0: pixel lies completely outside the queried shape
          1: pixel may overlap with the shape, pixel center is outside
          2: pixel center is inside the shape, but maybe not the complete pixel
          3: pixel lies completely inside the shape */

template<typename int64> inline bool check_pixel_ring 
(
 const T_Healpix_Base<int64> &b1,
 const T_Healpix_Base<int64> &b2, 
 int64 pix, 
 int64 nr, 
 int64 ipix1, 
 int fct,
 double cz, 
 double cphi, 
 double cosrp2, 
 int64 cpix
)
{
  if (pix>=nr) pix-=nr;
  if (pix<0) pix+=nr;
  pix+=ipix1;
  if (pix==cpix) return false; // disk center in pixel => overlap
  int px,py,pf;
  b1.pix2xyf(pix,px,py,pf);
  for (int i=0; i<fct-1; ++i) // go along the 4 edges
    {
    int64 ox=fct*px, oy=fct*py;
    double pz,pphi;
    b2.pix2zphi(b2.xyf2pix(ox+i,oy,pf),pz,pphi);
    if (cosdist_zphi(pz,pphi,cz,cphi)>cosrp2) // overlap
      return false;
    b2.pix2zphi(b2.xyf2pix(ox+fct-1,oy+i,pf),pz,pphi);
    if (cosdist_zphi(pz,pphi,cz,cphi)>cosrp2) // overlap
      return false;
    b2.pix2zphi(b2.xyf2pix(ox+fct-1-i,oy+fct-1,pf),pz,pphi);
    if (cosdist_zphi(pz,pphi,cz,cphi)>cosrp2) // overlap
      return false;
    b2.pix2zphi(b2.xyf2pix(ox,oy+fct-1-i,pf),pz,pphi);
    if (cosdist_zphi(pz,pphi,cz,cphi)>cosrp2) // overlap
      return false;
    }
  return true;
}

template<typename int64> inline void check_pixel 
(
 int o, 
 int order_,
 int omax, 
 int zone, 
 vector<pair<int64,int> > &pixset, 
 int64 pix, 
 vector<pair<int64,int> > &stk,
 bool inclusive, 
 int &stacktop
)
{
/* Short note on the "zone":
   zone = 0: pixel lies completely outside the queried shape
          1: pixel may overlap with the shape, pixel center is outside
          2: pixel center is inside the shape, but maybe not the complete pixel
          3: pixel lies completely inside the shape */
  if (zone==0) return; //pixel lies completely outside the queried shape

  if (o<order_)
    {
    if (zone>=3)//pixel lies completely inside the shape
      {
      pixset.push_back(make_pair(pix,o)); // output lowest level pixel
      }
    else // (zone>=1) pixel may overlap with the shape, pixel center is outside
      for (int i=0; i<4; ++i)
        stk.push_back(make_pair(4*pix+3-i,o+1)); // add children
    }
  else if (o>order_) // this implies that inclusive==true
    {
    if (zone>=2) // pixel center is inside the shape, but maybe not the complete pixel
      {
      pixset.push_back(make_pair(pix>>(2*(o-order_)),o)); // ROB output the parent pixel at order_
      stk.resize(stacktop); // unwind the stack
      }
    else // (zone>=1): pixel center in safety range
      {
      if (o<omax) // check sublevels
        for (int i=0; i<4; ++i) // add children in reverse order
          stk.push_back(make_pair(4*pix+3-i,o+1));
      else // at resolution limit
        {
        pixset.push_back(make_pair(pix>>(2*(o-order_)),o)); // output the parent pixel at order_
        stk.resize(stacktop); // unwind the stack
        }
      }
    }
  else // o==order_
    {
    if (zone>=2) // pixel center is inside the shape, but maybe not the complete pixel
      pixset.push_back(make_pair(pix,o)); 
    else if (inclusive) // and (zone>=1) pixel may overlap with the shape, pixel center is outside
      {
      if (order_<omax) // check sublevels
        {
        stacktop=stk.size(); // remember current stack position
        for (int i=0; i<4; ++i) // add children in reverse order
          stk.push_back(make_pair(4*pix+3-i,o+1));
        }
      else // at resolution limit
        pixset.push_back(make_pair(pix,o)); // output the pixel
      }
    }
}


//#########################################################
//#                 RANGE POINT QUERIES                   #
//#########################################################


//###
//### QUERY DISC
//###
void Healpix_Custom::query_disc
(
 pointing ptg, 
 double radius, 
 vector< pair<int64,int> > &pixset
) const
{
  query_disc_internal (ptg, radius, 0, pixset);
}

//###
//### QUERY DISC INCLUSIVE
//###
void Healpix_Custom::query_disc_inclusive
(
 pointing ptg, 
 double radius, 
 vector< pair<int64,int> > &pixset, 
 int fact
) const
{
  planck_assert(fact>0,"fact must be a positive integer");
  if ((sizeof(int64)<8) && (((int64(1)<<order_max)/nside_)<fact))
    {
    Healpix_Custom base2(nside_,scheme_,SET_NSIDE);
    base2.query_disc_internal(ptg,radius,fact,pixset);
    return;
    }
  query_disc_internal (ptg, radius, fact, pixset);
}

//###
//### QUERY DISC INTERNAL
//###
void Healpix_Custom::query_disc_internal
(
 pointing ptg, 
 double radius, 
 int fact, 
 vector< pair<int64,int> > &pixset
) const
{
  bool inclusive = (fact!=0);
  pixset.clear();
  ptg.normalize();

  if (scheme_==RING)
    {
    int64 fct=1;
    if (inclusive)
      {
      planck_assert (((int64(1)<<order_max)/nside_)>=fact,
        "invalid oversampling factor");
      fct = fact;
      }
    T_Healpix_Base b2;
    double rsmall, rbig;
    if (fct>1)
      {
      b2.SetNside(fct*nside_,RING);
      rsmall = radius+b2.max_pixrad();
      rbig = radius+max_pixrad();
      }
    else
      rsmall = rbig = inclusive ? radius+max_pixrad() : radius;

    if (rsmall>=pi)
	{ append_pixel_range(0,npix_,order_,pixset); return; }

    rbig = min(pi,rbig);

    double cosrsmall = cos(rsmall);
    double cosrbig = cos(rbig);

    double z0 = cos(ptg.theta);
    double xa = 1./sqrt((1-z0)*(1+z0));

    int64 cpix=zphi2pix(z0,ptg.phi);

    double rlat1 = ptg.theta - rsmall;
    double zmax = cos(rlat1);
    int64 irmin = ring_above (zmax)+1;

    if ((rlat1<=0) && (irmin>1)) // north pole in the disk
      {
      int64 sp,rp; bool dummy;
      get_ring_info_small(irmin-1,sp,rp,dummy);
	  append_pixel_range(0,sp+rp,order_,pixset);
      }

    if ((fct>1) && (rlat1>0)) irmin=max(int64(1),irmin-1);

    double rlat2 = ptg.theta + rsmall;
    double zmin = cos(rlat2);
    int64 irmax = ring_above (zmin);

    if ((fct>1) && (rlat2<pi)) irmax=min(4*nside_-1,irmax+1);

    for (int64 iz=irmin; iz<=irmax; ++iz)
      {
      double z=ring2z(iz);
      double x = (cosrbig-z*z0)*xa;
      double ysq = 1-z*z-x*x;
      double dphi = (ysq<=0) ? pi-1e-15 : atan2(sqrt(ysq),x);
      int64 nr, ipix1;
      bool shifted;
      get_ring_info_small(iz,ipix1,nr,shifted);
      double shift = shifted ? 0.5 : 0.;

      int64 ipix2 = ipix1 + nr - 1; // highest pixel number in the ring

      int64 ip_lo = ifloor<int64>(nr*inv_twopi*(ptg.phi-dphi) - shift)+1;
      int64 ip_hi = ifloor<int64>(nr*inv_twopi*(ptg.phi+dphi) - shift);

      if (fct>1)
        {
        while ((ip_lo<=ip_hi) && check_pixel_ring
               (*this,b2,ip_lo,nr,ipix1,fct,z0,ptg.phi,cosrsmall,cpix))
          ++ip_lo;
        while ((ip_hi>ip_lo) && check_pixel_ring
               (*this,b2,ip_hi,nr,ipix1,fct,z0,ptg.phi,cosrsmall,cpix))
          --ip_hi;
        }

      if (ip_lo<=ip_hi)
        {
        if (ip_hi>=nr)
          { ip_lo-=nr; ip_hi-=nr; }
        if (ip_lo<0)
          {
          append_pixel_range(ipix1,ipix1+ip_hi+1,order_,pixset);
		  append_pixel_range(ipix1+ip_lo+nr,ipix2+1,order_,pixset);
  		  }
        else
          append_pixel_range(ipix1+ip_lo,ipix1+ip_hi+1,order_,pixset);
        }
      }
    if ((rlat2>=pi) && (irmax+1<4*nside_)) // south pole in the disk
      {
      int64 sp,rp; bool dummy;
      get_ring_info_small(irmax+1,sp,rp,dummy);
      append_pixel_range(sp,npix_,order_,pixset);
      }
    }
  else // scheme_==NEST
    {
    if (radius>=pi) // disk covers the whole sphere
      { append_pixel_range(0,npix_,order_,pixset); return; }

    int oplus = 0;
    if (inclusive)
      {
      planck_assert ((int64(1)<<(order_max-order_))>=fact,
        "invalid oversampling factor");
      planck_assert ((fact&(fact-1))==0,
        "oversampling factor must be a power of 2");
      oplus=ilog2(fact);
      }
    int omax=order_+oplus; // the order up to which we test

    vec3 vptg(ptg);
    arr<T_Healpix_Base<int64> > base(omax+1);
    arr<double> crpdr(omax+1), crmdr(omax+1);
    double cosrad=cos(radius);
    for (int o=0; o<=omax; ++o) // prepare data at the required orders
      {
      base[o].Set(o,NEST);
      double dr=base[o].max_pixrad(); // safety distance
      crpdr[o] = (radius+dr>pi) ? -1. : cos(radius+dr);
      crmdr[o] = (radius-dr<0.) ?  1. : cos(radius-dr);
      }
    vector<pair<int64,int> > stk; // stack for pixel numbers and their orders
    stk.reserve(12+3*omax); // reserve maximum size to avoid reallocation
    for (int i=0; i<12; ++i) // insert the 12 base pixels in reverse order
      stk.push_back(make_pair(int64(11-i),0));

    int stacktop=0; // a place to save a stack position

    while (!stk.empty()) // as long as there are pixels on the stack
      {
      // pop current pixel number and order from the stack
      int64 pix=stk.back().first;
      int o=stk.back().second;
      stk.pop_back();

      double z,phi;
      base[o].pix2zphi(pix,z,phi);
      // cosine of angular distance between pixel center and disk center
      double cangdist=cosdist_zphi(vptg.z,ptg.phi,z,phi);

      if (cangdist>crpdr[o])
        {
        int zone = (cangdist<cosrad) ? 1 : ((cangdist<=crmdr[o]) ? 2 : 3);

        check_pixel (o, order_, omax, zone, pixset, pix, stk, inclusive,
          stacktop);
        }
      }
    }
}

//###
//### QUERY MULTI-DISC 
//###
void Healpix_Custom::query_multidisc 
( 
 const arr<vec3> &norm,
 const arr<double> &rad, 
 int fact, 
 vector< pair<int64,int> > &pixset
) const
{
  bool inclusive = (fact!=0);
  tsize nv=norm.size();
  planck_assert(nv==rad.size(),"inconsistent input arrays");
  pixset.clear();

  if (scheme_==RING)
    {
    int64 fct=1;
    if (inclusive)
      {
      planck_assert (((int64(1)<<order_max)/nside_)>=fact,
        "invalid oversampling factor");
      fct = fact;
      }
    T_Healpix_Base b2;
    double rpsmall, rpbig;
    if (fct>1)
      {
      b2.SetNside(fct*nside_,RING);
      rpsmall = b2.max_pixrad();
      rpbig = max_pixrad();
      }
    else
      rpsmall = rpbig = inclusive ? max_pixrad() : 0;

    int64 irmin=1, irmax=4*nside_-1;
    vector<double> z0,xa,cosrsmall,cosrbig;
    vector<pointing> ptg;
    vector<int64> cpix;
    for (tsize i=0; i<nv; ++i)
      {
      double rsmall=rad[i]+rpsmall;
      if (rsmall<pi)
        {
        double rbig=min(pi,rad[i]+rpbig);
        pointing pnt=pointing(norm[i]);
        cosrsmall.push_back(cos(rsmall));
        cosrbig.push_back(cos(rbig));
        double cth=cos(pnt.theta);
        z0.push_back(cth);
        if (fct>1) cpix.push_back(zphi2pix(cth,pnt.phi));
        xa.push_back(1./sqrt((1-cth)*(1+cth)));
        ptg.push_back(pnt);

        double rlat1 = pnt.theta - rsmall;
        double zmax = cos(rlat1);
        int64 irmin_t = (rlat1<=0) ? 1 : ring_above (zmax)+1;

        if ((fct>1) && (rlat1>0)) irmin_t=max(int64(1),irmin_t-1);

        double rlat2 = pnt.theta + rsmall;
        double zmin = cos(rlat2);
        int64 irmax_t = (rlat2>=pi) ? 4*nside_-1 : ring_above (zmin);

        if ((fct>1) && (rlat2<pi)) irmax_t=min(4*nside_-1,irmax_t+1);

        if (irmax_t < irmax) irmax=irmax_t;
        if (irmin_t > irmin) irmin=irmin_t;
        }
      }

    for (int64 iz=irmin; iz<=irmax; ++iz)
      {
      double z=ring2z(iz);
      int64 ipix1,nr;
      bool shifted;
      get_ring_info_small(iz,ipix1,nr,shifted);
      double shift = shifted ? 0.5 : 0.;
      rangeset<int64> tr;
      tr.append(ipix1,ipix1+nr);
      for (tsize j=0; j<z0.size(); ++j)
        {
        double x = (cosrbig[j]-z*z0[j])*xa[j];
        double ysq = 1.-z*z-x*x;
        double dphi = (ysq<=0) ? pi-1e-15 : atan2(sqrt(ysq),x);
        int64 ip_lo = ifloor<int64>(nr*inv_twopi*(ptg[j].phi-dphi) - shift)+1;
        int64 ip_hi = ifloor<int64>(nr*inv_twopi*(ptg[j].phi+dphi) - shift);
        if (fct>1)
          {
          while ((ip_lo<=ip_hi) && check_pixel_ring
            (*this,b2,ip_lo,nr,ipix1,fct,z0[j],ptg[j].phi,cosrsmall[j],cpix[j]))
            ++ip_lo;
          while ((ip_hi>ip_lo) && check_pixel_ring
            (*this,b2,ip_hi,nr,ipix1,fct,z0[j],ptg[j].phi,cosrsmall[j],cpix[j]))
            --ip_hi;
          }
        if (ip_hi>=nr)
          { ip_lo-=nr; ip_hi-=nr;}
        if (ip_lo<0)
          tr.remove(ipix1+ip_hi+1,ipix1+ip_lo+nr);
        else
          tr.intersect(ipix1+ip_lo,ipix1+ip_hi+1);
        }
	  append_pixel_range(tr.ivbegin(0),tr.ivend(0),order_,pixset);
      }
    }
  else // scheme_ == NEST
    {
    int oplus = 0;
    if (inclusive)
      {
      planck_assert ((int64(1)<<(order_max-order_))>=fact,
        "invalid oversampling factor");
      planck_assert ((fact&(fact-1))==0,
        "oversampling factor must be a power of 2");
      oplus=ilog2(fact);
      }
    int omax=order_+oplus; // the order up to which we test

    // TODO: ignore all disks with radius>=pi

    arr<T_Healpix_Base<int64> > base(omax+1);
    arr3<double> crlimit(omax+1,nv,3);
    for (int o=0; o<=omax; ++o) // prepare data at the required orders
      {
      base[o].Set(o,NEST);
      double dr=base[o].max_pixrad(); // safety distance
      for (tsize i=0; i<nv; ++i)
        {
        crlimit(o,i,0) = (rad[i]+dr>pi) ? -1. : cos(rad[i]+dr);
        crlimit(o,i,1) = (o==0) ? cos(rad[i]) : crlimit(0,i,1);
        crlimit(o,i,2) = (rad[i]-dr<0.) ?  1. : cos(rad[i]-dr);
        }
      }

    vector<pair<int64,int> > stk; // stack for pixel numbers and their orders
    stk.reserve(12+3*omax); // reserve maximum size to avoid reallocation
    for (int i=0; i<12; ++i) // insert the 12 base pixels in reverse order
      stk.push_back(make_pair(int64(11-i),0));

    int stacktop=0; // a place to save a stack position

    while (!stk.empty()) // as long as there are pixels on the stack
      {
      // pop current pixel number and order from the stack
      int64 pix=stk.back().first;
      int o=stk.back().second;
      stk.pop_back();

      vec3 pv(base[o].pix2vec(pix));

      tsize zone=3;
      for (tsize i=0; i<nv; ++i)
        {
        double crad=dotprod(pv,norm[i]);
        for (tsize iz=0; iz<zone; ++iz)
          if (crad<crlimit(o,i,iz))
            if ((zone=iz)==0) goto bailout;
        }

      check_pixel (o, order_, omax, zone, pixset, pix, stk, inclusive,
        stacktop);
      bailout:;
      }
    }
}

//###
//### QUERY MULTI-DISC GENERAL
//###
void Healpix_Custom::query_multidisc_general 
( 
 const arr<vec3> &norm, 
 const arr<double> &rad,
 bool inclusive, 
 const std::vector<int> &cmds, 
 vector< pair<int64,int> > &pixset
) const
{
  tsize nv=norm.size();
  planck_assert(nv==rad.size(),"inconsistent input arrays");
  pixset.clear();

  if (scheme_==RING)
    {
    planck_fail ("not yet implemented");
    }
  else // scheme_ == NEST
    {
    int oplus=inclusive ? 2 : 0;
    int omax=min(order_max,order_+oplus); // the order up to which we test

    // TODO: ignore all disks with radius>=pi

    arr<T_Healpix_Base<int64> > base(omax+1);
    arr3<double> crlimit(omax+1,nv,3);
    for (int o=0; o<=omax; ++o) // prepare data at the required orders
      {
      base[o].Set(o,NEST);
      double dr=base[o].max_pixrad(); // safety distance
      for (tsize i=0; i<nv; ++i)
        {
        crlimit(o,i,0) = (rad[i]+dr>pi) ? -1. : cos(rad[i]+dr);
        crlimit(o,i,1) = (o==0) ? cos(rad[i]) : crlimit(0,i,1);
        crlimit(o,i,2) = (rad[i]-dr<0.) ?  1. : cos(rad[i]-dr);
        }
      }

    vector<pair<int64,int> > stk; // stack for pixel numbers and their orders
    stk.reserve(12+3*omax); // reserve maximum size to avoid reallocation
    for (int i=0; i<12; ++i) // insert the 12 base pixels in reverse order
      stk.push_back(make_pair(int64(11-i),0));

    int stacktop=0; // a place to save a stack position
    arr<tsize> zone(nv);

    vector<tsize> zstk; zstk.reserve(cmds.size());

    while (!stk.empty()) // as long as there are pixels on the stack
      {
      // pop current pixel number and order from the stack
      int64 pix=stk.back().first;
      int o=stk.back().second;
      stk.pop_back();

      vec3 pv(base[o].pix2vec(pix));

      for (tsize i=0; i<nv; ++i)
        {
        zone[i]=3;
        double crad=dotprod(pv,norm[i]);
        for (tsize iz=0; iz<zone[i]; ++iz)
          if (crad<crlimit(o,i,iz))
            zone[i]=iz;
        }

      for (tsize i=0; i<cmds.size(); ++i)
        {
        tsize tmp;
        switch (cmds[i])
          {
          case -1: // union
            tmp=zstk.back(); zstk.pop_back();
            zstk.back() = max(zstk.back(),tmp);
            break;
          case -2: // intersection
            tmp=zstk.back(); zstk.pop_back();
            zstk.back() = min(zstk.back(),tmp);
            break;
          default: // add value
            zstk.push_back(zone[cmds[i]]);
          }
        }
      planck_assert(zstk.size()==1,"inconsistent commands");
      tsize zn=zstk[0]; zstk.pop_back();

      check_pixel (o, order_, omax, zn, pixset, pix, stk, inclusive,
        stacktop);
      }
    }
}

//###
//### QUERY POLYGON INTERNAL
//###
void Healpix_Custom::query_polygon_internal 
( 
 const std::vector<pointing> &vertex, 
 int fact,
 vector< pair<int64,int> > &pixset
) const
{
  bool inclusive = (fact!=0);
  tsize nv=vertex.size();
  tsize ncirc = inclusive ? nv+1 : nv;
  planck_assert(nv>=3,"not enough vertices in polygon");
  arr<vec3> vv(nv);
  for (tsize i=0; i<nv; ++i)
    vv[i]=vertex[i].to_vec3();
  arr<vec3> normal(ncirc);
  int flip=0;
  for (tsize i=0; i<nv; ++i)
    {
    normal[i]=crossprod(vv[i],vv[(i+1)%nv]);
    double hnd=dotprod(normal[i],vv[(i+2)%nv]);
    planck_assert(abs(hnd)>1e-10,"degenerate corner");
    if (i==0)
      flip = (hnd<0.) ? -1 : 1;
    else
      planck_assert(flip*hnd>0,"polygon is not convex");
    normal[i]*=flip/normal[i].Length();
    }
  arr<double> rad(ncirc,halfpi);
  if (inclusive)
    {
    double cosrad;
    find_enclosing_circle (vv, normal[nv], cosrad);
    rad[nv]=acos(cosrad);
    }
  this->query_multidisc(normal,rad,fact,pixset);
}

//###
//### QUERY POLYGON 
//###
void Healpix_Custom::query_polygon 
( 
 const std::vector<pointing> &vertex,
 vector< pair<int64,int> > &pixset
) const
{
  query_polygon_internal(vertex, 0, pixset);
}

//###
//### QUERY POLYGON INCLUSIVE
//###
void Healpix_Custom::query_polygon_inclusive 
( 
 const std::vector<pointing> &vertex,
 vector< pair<int64,int> > &pixset, 
 int fact
) const
{
  planck_assert(fact>0,"fact must be a positive integer");
  if ((sizeof(int64)<8) && (((int64(1)<<order_max)/nside_)<fact))
    {
    Healpix_Custom base2(nside_,scheme_,SET_NSIDE);
    base2.query_polygon_internal(vertex,fact,pixset);
    return;
    }
  query_polygon_internal(vertex, fact, pixset);
}

//###
//### QUERY STRIP INTERNAL 
//###
void Healpix_Custom::query_strip_internal
(
 double theta1, 
 double theta2, 
 bool inclusive, 
 vector< pair<int64,int> > &pixset
) const
{
	if (scheme_==RING)
    {
    int64 ring1 = max(int64(1),1+ring_above(cos(theta1))),
      ring2 = min(4*nside_-1,ring_above(cos(theta2)));
    if (inclusive)
      {
      ring1 = max(int64(1),ring1-1);
      ring2 = min(4*nside_-1,ring2+1);
      }

    int64 sp1,rp1,sp2,rp2;
    bool dummy;
    get_ring_info_small(ring1,sp1,rp1,dummy);
    get_ring_info_small(ring2,sp2,rp2,dummy);
    int64 pix1 = sp1,
      pix2 = sp2+rp2;

    if (pix1<=pix2) append_pixel_range(pix1,pix2,order_,pixset);
    }
  else
    planck_fail("query_strip not yet implemented for NESTED");
}

//###
//### QUERY STRIP
//###
void Healpix_Custom::query_strip 
(
 double theta1,
 double theta2, 
 bool inclusive, 
 vector< pair<int64,int> > &pixset
) const
{
  pixset.clear();

  if (theta1<theta2)
    query_strip_internal(theta1,theta2,inclusive,pixset);
  else
  {
    query_strip_internal(0.,theta2,inclusive,pixset);
    vector<pair<int64,int>> ps2;
    query_strip_internal(theta1,pi,inclusive,ps2);
	for(unsigned int i = 0; i < ps2.size(); i++) {
		pixset.push_back(ps2[i]);
	}
  }
}

//###
//### QUERY NEIGHBORS
//###
void Healpix_Custom::neighbors 
(
 pair<int64,int> pix,
 fix_arr<pair<int64,int>,8> &result
) const
{
  int ix, iy, fn;
  (scheme_==RING) ? ring2xyf(pix.first,ix,iy,fn) : nest2xyf(pix.first,ix,iy,fn);

  const int64 nsm1 = nside_-1;
  if ((ix>0)&&(ix<nsm1)&&(iy>0)&&(iy<nsm1))
  {
    if (scheme_==RING)
	 {
      for (int m=0; m<8; ++m)
	   {
        result[m].first = xyf2ring(ix+nb_xoffset[m],iy+nb_yoffset[m],fn);
	     result[m].second = order_;
	   }
	 }
    else
    {
         int64 fpix = int64(fn)<<(2*order_),
         px0=spread_bits(ix  ), py0=spread_bits(iy  )<<1,
         pxp=spread_bits(ix+1), pyp=spread_bits(iy+1)<<1,
         pxm=spread_bits(ix-1), pym=spread_bits(iy-1)<<1;

         result[0].first = fpix+pxm+py0; result[0].second = order_;
		 result[1].first = fpix+pxm+pyp; result[1].second = order_;
		 result[2].first = fpix+px0+pyp; result[2].second = order_;
		 result[3].first = fpix+pxp+pyp; result[3].second = order_;
		 result[4].first = fpix+pxp+py0; result[4].second = order_; 
		 result[5].first = fpix+pxp+pym; result[5].second = order_;
		 result[6].first = fpix+px0+pym; result[6].second = order_; 
		 result[7].first = fpix+pxm+pym; result[7].second = order_;
    }
  }
  else
  {
    for (int i=0; i<8; ++i)
    {
      int x=ix+nb_xoffset[i], y=iy+nb_yoffset[i];
      int nbnum=4;
      if (x<0)
      { x+=nside_; nbnum-=1; }
      else if (x>=nside_)
      { x-=nside_; nbnum+=1; }
      
		if (y<0)
      { y+=nside_; nbnum-=3; }
      else if (y>=nside_)
      { y-=nside_; nbnum+=3; }

      int f = nb_facearray[nbnum][fn];
      if (f>=0)
      {
        int bits = nb_swaparray[nbnum][fn>>2];
        if (bits&1) x=nside_-x-1;
        if (bits&2) y=nside_-y-1;
        if (bits&4) std::swap(x,y);
        result[i].first = (scheme_==RING) ? xyf2ring(x,y,f) : xyf2nest(x,y,f);
		result[i].second = order_;
      }
      else
	   {
        result[i].first = -1;
        result[i].second = order_;
	   }
    }
  }
}

