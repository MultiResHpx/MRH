
#ifndef __POINT2D_H__
#define __POINT2D_H__

class Point2D
{
public:
   Point2D () {};
   Point2D( float nx, float ny, int fn, int mapidx );
   
   float nx;
   float ny;
   int fn;
   int mapidx;
};

#endif