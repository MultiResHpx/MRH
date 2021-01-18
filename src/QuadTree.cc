// Code adapted from original downloaded on 15-APRIL-2014 from:
// http://veendeta.wordpress.com/2012/10/26/algorithm-adaptive-quadtree/


#include "QuadTree.h"
#include "Point2D.h"

QuadTree::QuadTree() :
   LEFT( 0.0 ),
   RIGHT( 1.0 ),
   TOP( 1.0 ),
   BOTTOM( 0.0 ),
   MAX_LEAF_ELEMENTS( 1 ),
   NODES( 0 ),
   isLeaf( true )
{
}

QuadTree::QuadTree( float _left, float _right, float _bottom, float _top, unsigned int _MAX_LEAF_ELEMENTS ) :
   LEFT( _left ),
   RIGHT( _right ),
   TOP( _top ),
   BOTTOM( _bottom ),
   MAX_LEAF_ELEMENTS( _MAX_LEAF_ELEMENTS ),
   NODES( 0 ),
   isLeaf( true )
{
}

QuadTree::QuadTree( int _NSIDE, unsigned int _MAX_LEAF_ELEMENTS ) :
   LEFT( 0.0 ),
   RIGHT( 1.0 ),
   BOTTOM( 0.0 ),
   TOP( 1.0 ),
   NSIDE( _NSIDE ),
   MAX_LEAF_ELEMENTS( _MAX_LEAF_ELEMENTS ),
   NODES( 0 ),
   isLeaf( true )
{
}

QuadTree::~QuadTree()
{
   if ( !isLeaf )
      delete [] NODES;
}

int QuadTree::Insert( Point2D point )
{
   return _Insert(point,0);
}

int QuadTree::_Insert( Point2D point, int curDepth)
{
   if ( isLeaf ) 
   {
      POINTS.push_back( point );
      if ( POINTS.size() > MAX_LEAF_ELEMENTS ) 
      {
         createLeaves();
         moveObjectsToLeaves();
         isLeaf = false;
      }
      return curDepth;
   }

   NODES[whichQuad(point)].Insert(point);
   curDepth += 1;

   return curDepth;

}

int QuadTree::whichQuad(Point2D point)const
{
	return whichQuad(point.nx,point.ny);
}

int QuadTree::whichQuad( float x, float y)const
{
	// Test which Quadrant the point lies inside.
	// In the case of point lying on a vertical border
	// the point is defined to belong to the Quadrants LEFT 
	// of the vertical border. 
	// In the case of point lying on a horizontal border 
	// the point is defined to belong to the Quadrants BELOW
	// the horizontal border.

	// Test if point lies in NW Quad
	if( (x >= LEFT && x <= (LEFT+RIGHT)/2) &&
		(y <= TOP  && y > (TOP+BOTTOM)/2) )
	{
		return NW; 
	} 
	// Test if point lies in NE Quad
	else if( (x > (LEFT+RIGHT)/2 && x <= RIGHT) &&
		     (y <= TOP && y > (TOP+BOTTOM)/2) )
	{
		return NE;
	} 
	// Test if point lies in SW Quad
	if( (x >= LEFT && x <= (LEFT+RIGHT)/2) &&
		(y <= (TOP+BOTTOM)/2  && y >= BOTTOM) )
	{
		return SW;
	} 
	// Test if point lies in SE Quad
	else if( (x > (LEFT+RIGHT)/2 && x <= RIGHT) &&
		     (y <= (TOP+BOTTOM)/2  && y >= BOTTOM) )
	{
		return SE;
	}
}


void QuadTree::Clear()
{
   POINTS.clear();
   
   if ( !isLeaf ) 
   {
      for ( int n = 0; n < NodeCount; ++n ) 
      {
         NODES[n].Clear();
      }
   }
}

vector<Point2D> QuadTree::operator() (float nx, float ny)const
{
   return At(nx,ny);
}

vector<Point2D> QuadTree::At( float nx, float ny ) const
{
   if ( isLeaf ) 
   {
      return POINTS;
   }
   
   vector<Point2D> returnedObjects;
   vector<Point2D> childObjects;
   Point2D queryObject;
   queryObject.nx = nx;
   queryObject.ny = ny;
   
   if ( !POINTS.empty() )
      returnedObjects.insert( returnedObjects.end(), POINTS.begin(), POINTS.end() );
   
   childObjects = NODES[whichQuad(nx,ny)].At( nx, ny );
   for( unsigned int i = 0; i < childObjects.size(); i++) 
   {
     queryObject.mapidx = childObjects[i].mapidx;
     returnedObjects.push_back(queryObject);
   }
   
   return returnedObjects;
}

//
// Brute force search method to discover count of all the leaf nodes that contain "mapidx"
//
int QuadTree::CountLeavesWithMapidx(int mapidx) const
{
	std::vector<Point2D> fPoint2D;
	int count = 0;
	for( int ix = 0; ix < NSIDE; ix++ ) 
	{
		for( int iy = 0; iy < NSIDE; iy++ )
		{
			fPoint2D = this->At(ix,iy);
			if( fPoint2D.size() > 0 ) 
			{
               count++;
			   fPoint2D.clear();
			}
		}
	}
	return count;
}

//
// Brute force search method to discover all the leaf nodes that contain "mapidx"
//
std::vector<Point2D> QuadTree::LeavesWithMapidx(int mapidx)const
{
	float nx,ny;
	std::vector<Point2D> fLeaves;
	std::vector<Point2D> fPoint2D;
	fLeaves.clear();
	for( int ix = 0; ix < NSIDE; ix++ ) 
	{
		for( int iy = 0; iy < NSIDE; iy++ )
		{
			// Normalize index
            nx = float(ix)/float(NSIDE);
			ny = float(iy)/float(NSIDE);
			fPoint2D = this->At(nx,ny);
			if( fPoint2D.size() > 0 ) 
			{
				if( fPoint2D[0].mapidx == mapidx )
				{
					fLeaves.insert(fLeaves.end(),fPoint2D.begin(),fPoint2D.end());
				}
				fPoint2D.clear();
			}
		}
	}
	return fLeaves;
}


void QuadTree::createLeaves()
{
   NODES = new QuadTree[4];
   NODES[NW] = QuadTree( LEFT, (LEFT+RIGHT)/2, (TOP+BOTTOM)/2, TOP, MAX_LEAF_ELEMENTS );
   NODES[NE] = QuadTree( (LEFT+RIGHT)/2, RIGHT, (TOP+BOTTOM)/2, TOP, MAX_LEAF_ELEMENTS );
   NODES[SW] = QuadTree( LEFT, (LEFT+RIGHT)/2, BOTTOM, (TOP+BOTTOM)/2, MAX_LEAF_ELEMENTS );
   NODES[SE] = QuadTree( (LEFT+RIGHT)/2, RIGHT, BOTTOM, (TOP+BOTTOM)/2, MAX_LEAF_ELEMENTS );
}

void QuadTree::moveObjectsToLeaves()
{
   for ( int n = 0; n < NodeCount; ++n ) 
   {
      for ( unsigned int m = 0; m < POINTS.size(); ++m ) 
      {
	     if( whichQuad(POINTS[m]) == NW )
         {
            NODES[NW].Insert( POINTS[m] );
		 }
		 else if( whichQuad(POINTS[m]) == NE )
		 {
            NODES[NE].Insert( POINTS[m] );
		 }
		 else if( whichQuad(POINTS[m]) == SW )
		 {
            NODES[SW].Insert( POINTS[m] );
		 } 
		 else if( whichQuad(POINTS[m]) == SE )
		 {
            NODES[SE].Insert( POINTS[m] );
		 }
		 POINTS.erase( POINTS.begin() + m );
         --m;
      }
   }
}

int QuadTree::bytes()
{
  return sizeof(NODES);
}

void QuadTree::SetNside(int nside)
{
   NSIDE = nside;
}
