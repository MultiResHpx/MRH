#include "TestQuadTree.h"
#include "Point2D.h"


void SimpleQuadA()
{
	srand(time(NULL));

	int rangeX,rangeY,qX,qY,cX,cY;
	//QuadTree tree( 0, MAX_NSIDE_32, 0, MAX_NSIDE_32 ) ;
	vector<Point2D> getObjects;
	vector<Point2D> fobjects;
    int objCenterX = 0;
	int objCenterY = 0;
	long int hpxID = 0;
	// ######### TEST #1 QUADTREE A: INSERT KNOWN SPATIAL DATA INTO QUADTREE #############
    // Create simple test set for V&V
	QuadTree treeA( 0, 15, 0, 15 ) ;
	vector<Point2D> objectsA;
	getObjects.clear();
	printf("\n##### QUADTREE A: #####\n");
	objectsA.push_back( Point2D(  7, 3, 0, 1) );
	objectsA.push_back( Point2D(  4, 6, 0, 2 ) );
	objectsA.push_back( Point2D( 10,12, 0, 3 ) );
	objectsA.push_back( Point2D(  3,10, 0, 4 ) );
	objectsA.push_back( Point2D(  5, 4, 0, 5 ) );
	objectsA.push_back( Point2D( 12, 3, 0, 6 ) );
	objectsA.push_back( Point2D(  5,10, 0, 7 ) );
	objectsA.push_back( Point2D(  2,10, 0, 8 ) );
	objectsA.push_back( Point2D(  9,10, 0, 9 ) );
	objectsA.push_back( Point2D(  8,10, 0, 10 ) );

	for ( unsigned int n = 0; n < objectsA.size(); n++ ) 
	{
		printf("Inserting Point2D: %d %d %d into QuadTreeA\n",objectsA[n].ix,objectsA[n].iy,objectsA[n].mapidx);
	    treeA.Insert( objectsA[n] );
	}

	// Retrieve objects from the quad tree
	printf("\n\n##### Now retrieve spatial data from QuadTreeA #####\n");
	for( int i = 0; i < objectsA.size(); i++)
	{
		getObjects = treeA.At( objectsA[i].ix, objectsA[i].iy );
		printf("Under %d %d found %d Objects:\n",objectsA[i].ix, objectsA[i].iy,getObjects.size());
		for( int j = 0; j < getObjects.size(); j++ )
		{
		   printf( "   Found Point2D %d %d %d\n", getObjects[j].ix, getObjects[j].iy, getObjects[j].mapidx );
		}
		getObjects.clear();
	}

	// Test for empty nodes
	//printf("\n##### Test for empty leaf nodes #####\n");
	//for( int i = 0; i < objectsA.size(); i++ )
	//{
	//    objCenterX = int( float(15) * (float(rand())/float(RAND_MAX)) );
 //	    objCenterY = int( float(15) * (float(rand())/float(RAND_MAX)) );
	//    getObjects = treeA.At( objCenterX, objCenterY );
	//	printf("Under %d %d found %d Objects:\n",objCenterX, objCenterY,getObjects.size());
	//	for( int j = 0; j < getObjects.size(); j++ )
	//	{
	//	   printf( "   Found Point2D %d %d %d\n", getObjects[j].x, getObjects[j].y, getObjects[j].hpxID );
	//	}
	//	getObjects.clear();
	//}

    // Do simple square region query
	getObjects.clear();
	rangeX = 8;
	rangeY = 8;
	cX = 7;
	cY = 7;
	printf("\n##### Simple Range Query on QuadTreeA#####\n");
	printf("  Centered at %d %d with dimensions %d x %d\n",cX,cY,rangeX,rangeY);
	for( int i = 0; i < rangeX; i++)
	{
		for( int j = 0; j < rangeY; j++)
		{
			qX = cX - int(0.5*rangeX) + i;
			qY = cY - int(0.5*rangeY) + j;
			fobjects = treeA.At( qX, qY );
			getObjects.insert(getObjects.end(),fobjects.begin(),fobjects.end());
		}
		fobjects.clear();
	}
	printf("Found %d Objects:\n",getObjects.size());
    for( int k = 0; k < getObjects.size(); k++ )
    {
      printf( "   %d %d %d\n", getObjects[k].ix, getObjects[k].iy, getObjects[k].mapidx );
    }

}

void SimpleQuadB()
{
	srand(time(NULL));

	int rangeX,rangeY,qX,qY,cX,cY;
	//QuadTree tree( 0, MAX_NSIDE_32, 0, MAX_NSIDE_32 ) ;
	vector<Point2D> getObjects;
	vector<Point2D> fobjects;
    int objCenterX = 0;
	int objCenterY = 0;
	long int hpxID = 0;
	// ######### TEST #2 QUADTREE B: INSERT KNOWN SPATIAL DATA INTO QUADTREE #############
    // Create simple test set for V&V
	QuadTree treeB( 0, 15, 0, 15 ) ;
	vector<Point2D> objectsB;
	getObjects.clear();
	objectsB.clear();
	printf("\n##### QUADTREE B: #####\n");
	objectsB.push_back( Point2D(  9, 11, 0, 1) );
	objectsB.push_back( Point2D(  7,  1, 0, 2 ) );
	objectsB.push_back( Point2D( 10,  3, 0, 3 ) );
	objectsB.push_back( Point2D(  5, 10, 0, 4 ) );
	objectsB.push_back( Point2D(  5,  7, 0, 5 ) );
	objectsB.push_back( Point2D(  7,  3, 0, 6 ) );
	objectsB.push_back( Point2D(  1,  1, 0, 7 ) );
	objectsB.push_back( Point2D(  7,  8, 0, 8 ) );
	objectsB.push_back( Point2D( 11, 11, 0, 9 ) );
	objectsB.push_back( Point2D( 11,  5, 0, 10 ) );
	
	for ( unsigned int n = 0; n < objectsB.size(); n++ ) 
	{
		printf("Inserting Point2D: %d %d %d into QuadTreeB\n",objectsB[n].ix,objectsB[n].iy,objectsB[n].mapidx);
	    treeB.Insert( objectsB[n] );
	}

	// Retrieve objects from the quad tree
	printf("\n\n##### Now retrieve spatial data from QuadTreeB #####\n");
	for( int i = 0; i < objectsB.size(); i++)
	{
		getObjects = treeB.At( objectsB[i].ix, objectsB[i].iy );
		printf("Under %d %d found %d Objects:\n",objectsB[i].ix, objectsB[i].iy,getObjects.size());
		for( int j = 0; j < getObjects.size(); j++ )
		{
		   printf( "   Found Point2D %d %d %d\n", getObjects[j].ix, getObjects[j].iy, getObjects[j].mapidx );
		}
		getObjects.clear();
	}

	// Test for empty nodes
	//printf("\n##### Test for empty leaf nodes #####\n");
	//for( int i = 0; i < objectsB.size(); i++ )
	//{
	//    objCenterX = int( float(19) * (float(rand())/float(RAND_MAX)) );
 //	    objCenterY = int( float(19) * (float(rand())/float(RAND_MAX)) );
	//    getObjects = treeB.At( objCenterX, objCenterY );
	//	printf("Under %d %d found %d Objects:\n",objCenterX, objCenterY,getObjects.size());
	//	for( int j = 0; j < getObjects.size(); j++ )
	//	{
	//	   printf( "   Found Point2D %d %d %d\n", getObjects[j].x, getObjects[j].y, getObjects[j].hpxID );
	//	}
	//	getObjects.clear();
	//}

    // Do simple square region query
	rangeX = 8;
	rangeY = 8;
	cX = 7;
	cY = 7;
	printf("\n##### Simple Range Query on QuadTreeB#####\n");
	printf("  Centered at %d %d with dimensions %d x %d\n",cX,cY,rangeX,rangeY);
	for( int i = 0; i < rangeX; i++)
	{
		for( int j = 0; j < rangeY; j++)
		{
			qX = cX - int(0.5*rangeX) + i;
			qY = cY - int(0.5*rangeY) + j;
			fobjects = treeB.At( qX, qY );
			getObjects.insert(getObjects.end(),fobjects.begin(),fobjects.end());
		}
	}
	printf("Found %d Objects:\n",getObjects.size());
    for( int k = 0; k < getObjects.size(); k++ )
    {
      printf( "   %d %d %d\n", getObjects[k].ix, getObjects[k].iy, getObjects[k].mapidx );
    }

}

void SimpleQuadC()
{
	srand(time(NULL));

	int rangeX,rangeY,qX,qY,cX,cY;
	//QuadTree tree( 0, MAX_NSIDE_32, 0, MAX_NSIDE_32 ) ;
	vector<Point2D> getObjects;
	vector<Point2D> fobjects;
    int objCenterX = 0;
	int objCenterY = 0;
	long int hpxID = 0;
	// ######### TEST #2 QUADTREE C: INSERT KNOWN SPATIAL DATA INTO QUADTREE #############
    // Create simple test set for V&V
	QuadTree treeC( 0, 15, 0, 15 ) ;
	vector<Point2D> objectsC;
	getObjects.clear();
	printf("\n##### QUADTREE C: #####\n");
	objectsC.push_back( Point2D( 3, 1, 0, 1) );
	objectsC.push_back( Point2D(10,12, 0, 2 ) );
	objectsC.push_back( Point2D( 1, 5, 0, 3 ) );
	objectsC.push_back( Point2D( 3, 7, 0, 4 ) );
	objectsC.push_back( Point2D( 6, 9, 0, 5 ) );
	objectsC.push_back( Point2D(10, 9, 0, 6 ) );
	objectsC.push_back( Point2D( 6, 6, 0, 7 ) );
	objectsC.push_back( Point2D( 3, 4, 0, 8 ) );
	objectsC.push_back( Point2D(11, 4, 0, 9 ) );
	objectsC.push_back( Point2D( 4, 6, 0, 10 ) );

	for ( unsigned int n = 0; n < objectsC.size(); n++ ) 
	{
		printf("Inserting Point2D: %d %d %d into QuadTreeC\n",objectsC[n].ix,objectsC[n].iy,objectsC[n].mapidx);
	    treeC.Insert( objectsC[n] );
	}

	// Retrieve objects from the quad tree
	printf("\n\n##### Now retrieve spatial data from QuadTreeC #####\n");
	for( int i = 0; i < objectsC.size(); i++)
	{
		getObjects = treeC.At( objectsC[i].ix, objectsC[i].iy );
		printf("Under %d %d found %d Objects:\n",objectsC[i].ix, objectsC[i].iy,getObjects.size());
		for( int j = 0; j < getObjects.size(); j++ )
		{
		   printf( "   Found Point2D %d %d %d\n", getObjects[j].ix, getObjects[j].iy, getObjects[j].mapidx );
		}
		getObjects.clear();
	}

	// Test for empty nodes
	//printf("\n##### Test for empty leaf nodes #####\n");
	//for( int i = 0; i < objectsC.size(); i++ )
	//{
	//    objCenterX = int( float(24) * (float(rand())/float(RAND_MAX)) );
 //	    objCenterY = int( float(24) * (float(rand())/float(RAND_MAX)) );
	//    getObjects = treeC.At( objCenterX, objCenterY );
	//	printf("Under %d %d found %d Objects:\n",objCenterX, objCenterY,getObjects.size());
	//	for( int j = 0; j < getObjects.size(); j++ )
	//	{
	//	   printf( "   Found Point2D %d %d %d\n", getObjects[j].x, getObjects[j].y, getObjects[j].hpxID );
	//	}
	//	getObjects.clear();
	//}

	// Do simple square region query
	rangeX = 8;
	rangeY = 8;
	cX = 7;
	cY = 7;
	printf("\n##### Simple Range Query on QuadTreeC#####\n");
	printf("  Centered at %d %d with dimensions %d x %d\n",cX,cY,rangeX,rangeY);
	for( int i = 0; i < rangeX; i++)
	{
		for( int j = 0; j < rangeY; j++)
		{
			qX = cX - int(0.5*rangeX) + i;
			qY = cY - int(0.5*rangeY) + j;
			fobjects = treeC.At( qX, qY );
			getObjects.insert(getObjects.end(),fobjects.begin(),fobjects.end());
		}
	}
	printf("Found %d Objects:\n",getObjects.size());
    for( int k = 0; k < getObjects.size(); k++ )
    {
      printf( "   %d %d %d\n", getObjects[k].ix, getObjects[k].iy, getObjects[k].mapidx );
    }

}

