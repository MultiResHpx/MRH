

#ifndef __QUADTREE_H__
#define __QUADTREE_H__

#include <vector>
#include "QuadTree.h"
#include "Point2D.h"
using std::vector;


class QuadTree {

   enum Node {
   NW = 0,
   NE,
   SW,
   SE,
   NodeCount
   };

public:
	QuadTree();
	QuadTree( float LEFT, float RIGHT, float TOP, float BOTTOM, unsigned int MAX_LEAF_ELEMENTS = 1 );
	QuadTree( int nside, unsigned int MAX_LEAF_ELEMENTS = 1 );
	~QuadTree();
	int	Insert( Point2D point );   
	void	Clear();
	std::vector<Point2D>	At( float nx, float ny ) const;
	std::vector<Point2D> operator() (float nx, float ny)const;
	int CountLeavesWithMapidx(int mapidx)const;
	std::vector<Point2D> LeavesWithMapidx(int mapidx)const;
	int  bytes();

	// SETTERS
	void SetNside(int nside);

	// GETTERS

protected:
	int _Insert( Point2D point, int curDepth);

private:
   float LEFT;
   float RIGHT;
   float TOP;   
   float BOTTOM;
   int NSIDE;
   unsigned int MAX_LEAF_ELEMENTS;
	std::vector<Point2D>	POINTS;
   QuadTree * NODES;
   bool	isLeaf;
   int whichQuad(Point2D point)const;
   int whichQuad( float nx, float ny)const;
   void	createLeaves();
   void	moveObjectsToLeaves();
};

#endif 