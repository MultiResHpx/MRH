#include "MultiResHpx.h"
// STD includes
#include <time.h>
#include <vector>
#include <math.h>
#include <cstdlib>
#include "lsconstants.h"
#include "healpix_map.h"

#define MRHDEBUG true
#define EMPTY -99999


class HPXNode 
{
public:
	HPXNode(){longitude = -99999.0; latitude = -99999.0; data = EMPTY;};
	HPXNode(double _longitude, double _latitude, int _data)
	{
		longitude = _longitude;
		latitude = _latitude;
		data = _data;
	}
	~HPXNode(){};

	double longitude;
	double latitude;
	int data;
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


void TestBuildHpx() {
	Healpix_Custom hpxQ;
	Healpix_Map<HPXNode> mHPX(1,NEST);
	int64 nPix;
	//int64* bigarray = new int64[nPix];
	cout << "First check for 64bit compatability, should be able to represent Npix up to Order = 29!" << endl;
	for(int64 order = 0; order < 30; order++) {
		nPix = 12*((int64(1)<<order)<<order);
		cout << "Order: " << order << " , Npix: " << nPix << endl;
	}
	cout << "Second Show can set HPX Query object to full range of HPX Orders:" << endl;
	for(int64 order = 0; order < 30; order++) {
		cout << "Setting HPX Query to order: " << order << endl;
		hpxQ.Set(order,NEST);
		cout << "	SUCCESS!" << endl;
	}
	cout << "Next test the limit of HPX Map Array sizes:" << endl;
	for(int64 order = 0; order < 30; order++) {
		nPix = 12*((int64(1)<<order)<<order);
		cout << "Setting HPX Map to order: " << order << " , Npix: " << nPix << endl;
		mHPX.Set(order,NEST);
		//mHPX.SetNside(int64(1)<<order,NEST);
		//std::vector<int> bigVect(nPix);
		//cout << "Size of vector: " << bigVect.size() << endl;
		cout << "	SUCCESS!" << endl;

	}
			//mHPX.Set(order,NEST);

}


void TestMRHInsertNodeBase0(int max_depth)
{
    MultiResHpx myMLQT = MultiResHpx(max_depth,NEST);
     
     int index = 0;
	 int num_insert = 6;
	 pointing pt;
	 double nextLong[] = {50.63*D2R, 50.62*D2R, 50.63*D2R, 39.38*D2R, 45.00*D2R, 20.0*D2R};
	 double nextLat[] =  {53.44*D2R, 59.06*D2R, 33.75*D2R, 33.75*D2R, 78.75*D2R, 50.0*D2R};
	 int nextDat[] = {1,2,3,4,5,6};
	 if(MRHDEBUG)
	 {
		 cout << "Number Nodes to insert: " << num_insert << endl;
	 }
    for(int i = 0; i < num_insert; i++ )
	 {
     	if(MRHDEBUG)
		{
     	   cout << "Insert New Node: " << nextLong[i] << " " << nextLat[i] << " " << nextDat[i] << endl;
		}
		pt.theta = nextLat[i];
		pt.phi = nextLong[i];
		myMLQT.Insert(pt,nextDat[i]);
	 }
    myMLQT.PrintTreeAtIndex(0); 
}

void TestMRHInsertOneNodeSearchOneNode(int max_depth,double phi, double theta, int data, int64 hpxid, int order)
{

	vector<MortonNode> found;
	pointing pt;
	Healpix_Custom myHpxQ(1,NEST);

	MultiResHpx myMRH = MultiResHpx(max_depth,NEST);

    // Insert the GIS longitude, latitude pairs into MRH
    cout << "Insert Node: " << phi << " " << theta << data << endl << endl;
	pt.phi = phi;
	pt.theta = theta;
    myMRH.Insert(pt,data);

	// Now Search MRH for GIS longitude, latitude pairs
	cout << "\nSearching for: " << phi << " " << theta << endl;
	found = myMRH.Search(pt);
	for( unsigned int j = 0; j < found.size(); j++) {
		cout << "\nFound: " << found[j].phi << " " << found[j].theta << " ";
		for(unsigned int i = 0; i < found[j].data.size(); i++) cout << found[j].data[i] << " ";
		cout << endl;
	}


	// Now Search MRH via user provided HEALPix Index,Order pair.
	cout << "\nSearching for: " << hpxid << " " << order << endl;
	found.clear();
	found = myMRH.Search(hpxid,order,true);

	for( unsigned int i = 0; i < found.size(); i++ ) {
		cout << "\nFound: " << found[i].phi << " " << found[i].theta << " ";
		for(unsigned int j = 0; i < found[i].data.size(); i++) cout << found[i].data[j] << " ";
		cout << endl;
	}

	//Print all trees
	for(int i = 0; i < 12; i++ ) 
	{
		cout << "\n\n######################\n";
		cout << "#### BASE CELL " << i << " ####\n";
		cout << "######################" << endl;
        myMRH.PrintTreeAtIndex(i); 
	}
}


void TestMRHInsertNodeAllBase(int max_depth,int seed,int num_points)
{

	double LO_THETA = -90.0*D2R;
	double LO_PHI = -180.0*D2R;
	double HI_THETA = 90.0*D2R;
	double HI_PHI = 180.0*D2R;
	vector<pointing> points;
	pointing pt;
	vector<MortonNode> found;
	MultiResHpx myMRH = MultiResHpx(max_depth,NEST);
	srand(time(0));

	// Generate several random GIS longitude, latitude pairs
	for(int i = 0; i < num_points; i++)
	{
		pt.phi = LO_PHI + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(HI_PHI-LO_PHI)));
		pt.theta = LO_THETA + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(HI_THETA-LO_THETA)));
		if(MRHDEBUG)
		{
     	   cout << "Insert New Node: " << pt.phi << " " << pt.theta << i << endl;
		}	
		points.push_back(pt);
	}

    // Insert the GIS longitude, latitude pairs into MRH
	for(int i = 0; i < num_points; i++)
	{
		myMRH.Insert(points[i],i);
	}

	// Now Search MRH for GIS longitude, latitude pairs
	for(int i = 0; i < num_points; i++)
	{
		cout << "Searching for: " << points[i].phi << " " << points[i].theta << endl;
		found = myMRH.Search(points[i]);
		for( unsigned int j = 0; j < found.size(); j++) {
			cout << "Found: " << found[j].phi << " " << found[j].theta << endl;
		}
	}

	//Print all trees
	for(int i = 0; i < 12; i++ ) 
	{
		cout << "\n\n######################\n";
		cout << "#### BASE CELL " << i << " ####\n";
		cout << "######################" << endl;
        myMRH.PrintTreeAtIndex(i); 
	}
}

void TestMRHDeleteNode(int maxdepth,int num_insert, char* argv[])
{
	MultiResHpx myMRH = MultiResHpx(maxdepth,NEST);
	std::vector<double> temp;
	std::vector<pointing> points;
	std::vector<int> indices;
	pointing pt;
    int i,j,index = 0;
	int nextDat;

	if(MRHDEBUG)
	{
		cout << "Max Depth: " << maxdepth << " Num Points: " << num_insert << endl;
	}

	for( i = 0; i < num_insert; i++ )
	{
		temp.clear();
     	pt.phi = atof(argv[4+index]);
     	pt.theta  = atof(argv[4+index+1]);
     	nextDat  = atoi(argv[4+index+2]);
		points.push_back(pt);
		indices.push_back(nextDat);
     	index += 3;
     	if(MRHDEBUG)
		{
     	   cout << "Insert New Node: " 
			    << pt.phi << " "
				<< pt.theta  << " " 
				<< nextDat << endl;
		}
		myMRH.Insert(pt,nextDat);
	}     	    
	//Print all trees
	for( i = 0; i < 12; i++ ) 
	{
		cout << "\n\n######################\n";
		cout << "#### BASE CELL " << i << " ####\n";
		cout << "######################" << endl;
        myMRH.PrintTreeAtIndex(i); 
	}
    // Now delete nodes one by one and print tree each time
	cout << "\n\n### Now delete nodes one by one and print tree after each deletion.###\n\n";
    for( i = 0; i < points.size(); i++ )
	{
        if(MRHDEBUG)
		{
           cout << "\n\nDelete Node at: "
			    << points[i].phi << " "
			    << points[i].theta << " "
			    << indices[i];    
		}
		myMRH.Delete(points[i]);
		//Print all trees
		for( j = 0; j < 12; j++ ) 
		{
			cout << "\n\n######################\n";
			cout << "#### BASE CELL " << j << " ####\n";
			cout << "######################" << endl;
			myMRH.PrintTreeAtIndex(j); 
		}
	}
}

void TestMRHInsertNodeAllBaseThenDeleteAll(int max_depth,int seed,int num_points)
{

	double LO_THETA = -90.0*D2R;
	double LO_PHI = -180.0*D2R;
	double HI_THETA = 90.0*D2R;
	double HI_PHI = 180.0*D2R;
	pointing pt;
	int i,j = 0;
	vector<pointing> points;
	vector<MortonNode> found;
	MultiResHpx myMRH = MultiResHpx(max_depth,NEST);
	srand(time(0));

	// Generate several random GIS longitude, latitude pairs
	for(int i = 0; i < num_points; i++)
	{
		pt.phi = LO_PHI + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(HI_PHI-LO_PHI)));
		pt.theta = LO_THETA + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(HI_THETA-LO_THETA)));
		if(MRHDEBUG)
		{
     	   cout << "Insert New Node: " << pt.phi << " " << pt.theta << i << endl;
		}	
		points.push_back(pt);
	}

    // Insert the GIS longitude, latitude pairs into MRH
	for( i = 0; i < num_points; i++)
	{
		myMRH.Insert(points[i],i);
	}

	// Now Search MRH for GIS longitude, latitude pairs
	for( i = 0; i < num_points; i++)
	{
		cout << "Searching for: " << points[i].phi << " " << points[i].theta << endl;
		found = myMRH.Search(points[i]);
		for(  j = 0; j < found.size(); j++) {
			cout << "Found: " << found[j].phi << " " << found[j].theta << endl;
		}
	}

	//Print all trees
	for( i = 0; i < 12; i++ ) 
	{
		cout << "\n\n######################\n";
		cout << "#### BASE CELL " << i << " ####\n";
		cout << "######################" << endl;
        myMRH.PrintTreeAtIndex(i); 
	}

    // Now delete nodes one by one and print all trees each time
	cout << "\n\n### Now delete nodes one by one and print tree after each deletion.###\n\n";
    for( i = 0; i < points.size(); i++ )
	{
        if(MRHDEBUG)
		{
           cout << "\n\nDelete Node at: "
			    << points[i].phi << " "
			    << points[i].theta;    
		}
		myMRH.Delete(points[i]);
		//Print all trees
		for( j = 0; j < 12; j++ ) 
		{
			cout << "\n\n######################\n";
			cout << "#### BASE CELL " << j << " ####\n";
			cout << "######################" << endl;
			myMRH.PrintTreeAtIndex(j); 
		}
	}
}

void TestMRHTransformersAndInfo(int maxdepth,double longitude,double latitude,int order,int step)
{
	MultiResHpx myMRH = MultiResHpx(maxdepth,NEST);
    int ix,iy,fn;
	int64 ring;
	int64 hpxid,rhpxid,pid,spix,rpix;
	bool shift;
	double z,phi,ctheta,stheta,theta,maxrad;
	vec3 v;
	std::vector<vec3> out;
	pointing ang;

	//##### Test all the various HEALPix tranformations, back and forth.

	// First compute the appropriate HEALPix index given long,lat and order.
    ang.phi = longitude;
	ang.theta = latitude;

	cout << "#### Inputs: longitude = " << longitude << " latitude = " << latitude << " order = " << order
		 << endl << endl;

 	   // Ang2pix (const pointing &ang, int order, int64& hpxid);
	myMRH.Ang2pix(ang,order,hpxid);
	cout << "Convert GIS longitude,latitude to NESTED hpxid" << endl;
	cout << "  hpxid = " << hpxid << endl;

	   // Pix2ang (int64 hpxid, int order, pointing& pt);
	myMRH.Pix2ang(hpxid,order,ang);
	cout << "Now converting back to longitude,latitude:" << endl;
	cout << "  longitude = " << ang.phi << " latitude = " << ang.theta << endl << endl;

       // Boundaries (int pix, int order, tsize step, std::vector<vec3> &out)
	myMRH.Boundaries(hpxid,order,step,out);
	cout << "Boundaries of hpxid are: " << endl;
	for(unsigned int i =0; i < out.size(); i++) {
		cout << out[i].x << ", " << out[i].y << ", " << out[i].z << endl;
	}
	cout << endl;

	// Max_pixrad(order)
	maxrad = myMRH.Max_pixrad(order);
	cout << "Maximum angular distance (in radian) between any pixel center and its corners: \n" <<  maxrad
		 << endl << endl;

    // Pix2xyf(int64 hpxid,int order, int &ix, int &iy, int &fn);
    myMRH.Pix2xyf(hpxid,order,ix,iy,fn);
	cout << "Convert hpxid to ix,iy,fn:" << endl;
	cout << "  ix = " << ix << " iy = " << iy << " fn = " << fn << endl;
	
	// Xyf2pix(int ix, int iy, int fn, int order, int64& hpxid);
    myMRH.Xyf2pix(ix,iy,fn,order,hpxid);
	cout << "Now converting ix,iy,fn back to hpxid:" << endl;
	cout << "  hpxid = " << hpxid << endl << endl;


    // Pix2zphi (int64 hpxid, int order, double &z, double &phi);
	myMRH.Pix2zphi(hpxid,order,z,phi);
	cout << "Convert hpxid,order to z,phi:" << endl;
	cout << " z = " << z << " phi = " << phi << endl;

    // Zphi2pix (double z, double phi, int order, int64& hpxid);
	myMRH.Zphi2pix(z,phi,order,hpxid);
	cout << "Now converting z,phi back to hpxid:" << endl;
	cout << " hpxid = " << hpxid << endl << endl;
	

	// Pix2vec (int64 hpxid, int order, vec3& v);
	myMRH.Pix2vec(hpxid,order,v);
	cout << "Convert hpxid,order to vector:" << endl;
	cout << "  vector = " << v.x << "," << v.y << "," << v.z << endl;
	
    // Vec2pix (const vec3 &v,int order, int64& hpxid);
	myMRH.Vec2pix(v,order,hpxid);
	cout << "Converting vector back to hpxid:" << endl;
	cout << "  hpxid = " << hpxid << endl << endl;


	// Pix2ring (int64 hpxid, int order,int& ring);
    myMRH.Pix2ring(hpxid,order,ring);
	cout << "Compute what Ring NESTED hpxid maps to:" << endl;
	cout << "  ring = " << ring << endl;
	myMRH.Get_Ring_Info(ring,order,spix,rpix,ctheta,stheta,shift);
	cout << "Info about the ring:" << endl;
	cout << "  start hpxid = " << spix << " num pix = " << rpix 
		 << " cosine lat = " << ctheta << " sine lat = " << stheta
		 << " shifted? = " << shift << endl << endl;
	myMRH.Get_Ring_Info2(ring,order,spix,rpix,theta,shift);
	cout << "Less info about the ring:" << endl;
	cout << "  start hpxid = " << spix << " num pix = " << rpix 
		 << " ring lat = " << theta << " shifted? = " << shift  << endl << endl;
	myMRH.Get_Ring_Info_Small(ring,order,spix,rpix,shift);
	cout << "Even less info about the ring:" << endl;
	cout << "  start hpxid = " << spix << " num pix = " << rpix  << " shifted? = " << shift << endl << endl;

	// Max_pixrad(int ring,int order)
	maxrad = myMRH.Max_pixrad(ring,order);
	cout << "maximum angular distance (in radian) between any pixel center and its corners of ring: \n" <<  maxrad
		 << endl << endl;
   
	
    // Nest2ring (int64 hpxid, int order, int64& Rhpxid);
	myMRH.Nest2ring(hpxid,order,rhpxid);
	cout << "Convert NESTED hpxid to RING hpxid:" << endl;
	cout << "  ring hpxid = " << rhpxid << endl;
	
	// Ring2nest (int64 Rhpxid, int order, int64& Nhpxid);
    myMRH.Ring2nest(rhpxid,order,hpxid);
	cout << "Converting RING hpxid back to NESTED hpxid:" << endl;
	cout << " nested hpxid = " << hpxid << endl << endl;
   

    // Nest2peano (int64 hpxid, int order, int64& pid);
    myMRH.Nest2peano(hpxid,order,pid);
	cout << "Convert hpxid to peano index:" << endl;
	cout << " peano = " << pid << endl;

    // Peano2nest (int64 pid, int order, int64& Nhpxid);
    myMRH.Peano2nest(pid,order,hpxid);
	cout << "Converting peano index back to hpxid:" << endl;
	cout << " hpxid = " << hpxid << endl << endl;
}




void TestMRHQueryDiscRandom(int maxdepth,int seed,int num_points,double phi, double theta, double radius)
{
	double LO_THETA = 0.0;
	double LO_PHI = 0.0;
	double HI_THETA = pi;
	double HI_PHI = 2.0*pi;
	int i,j = 0;
	vector<pointing> points;
	pointing pt;
	vector<MortonNode> found;
	MultiResHpx myMRH = MultiResHpx(maxdepth,NEST);
	srand(time(0));

	// Generate several random HPX phi,theta pairs
	for(int i = 0; i < num_points; i++)
	{
		pt.phi = LO_PHI + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(HI_PHI-LO_PHI)));
		pt.theta = LO_THETA + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(HI_THETA-LO_THETA)));
		if(MRHDEBUG)
		{
     	   cout << "Insert New Node: " << pt.phi << " " << pt.theta << i << endl;
		}	
		points.push_back(pt);
	}

    // Insert the HPX phi, theta pairs into MRH
	for( i = 0; i < num_points; i++)
	{
		myMRH.Insert(points[i],i);
	}

	//Print all trees
	for( i = 0; i < 12; i++ ) 
	{
		cout << "\n\n######################\n";
		cout << "#### BASE CELL " << i << " ####\n";
		cout << "######################" << endl;
        myMRH.PrintTreeAtIndex(i); 
	}

	// Now perform disc query and output the results
	cout << endl << endl;
	cout << "################################\n";
	cout << "#### Now Perform Disc Query ####\n";
	cout << "################################\n";
	cout << "Query Center: phi: " << phi
		 << " theta: " << theta 
		 << " radius: " << radius << endl << endl;
	pt.phi = phi;
	pt.theta = theta;
    found = myMRH.QueryDisc(pt,radius);
	for(  j = 0; j < found.size(); j++) {
		cout << "Found: " << found[j].phi << " " << found[j].theta << " ";
		for( i = 0; i < found[j].data.size(); i++) cout << found[j].data[i] << " " << endl;
		cout << endl;
	}
	if(found.size() == 0) {
		cout << "No Results Found!\n";
	}
}

void TestMRHQueryDiscInclusiveRandom(int maxdepth,int seed,int num_points,double phi, double theta, double radius, int fact)
{
	double LO_THETA = -90.0*D2R;
	double LO_PHI = -180.0*D2R;
	double HI_THETA = 90.0*D2R;
	double HI_PHI = 180.0*D2R;
	int i,j = 0;
	vector<pointing> points;
	pointing pt;
	vector<MortonNode> found;
	MultiResHpx myMRH = MultiResHpx(maxdepth,NEST);
	srand(time(0));

	// Generate several random GIS longitude, latitude pairs
	for(int i = 0; i < num_points; i++)
	{
		pt.phi = LO_PHI + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(HI_PHI-LO_PHI)));
		pt.theta = LO_THETA + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(HI_THETA-LO_THETA)));
		if(MRHDEBUG)
		{
     	   cout << "Insert New Node: " << pt.phi << " " << pt.theta << i << endl;
		}	
		points.push_back(pt);
	}

    // Insert the GIS longitude, latitude pairs into MRH
	for( i = 0; i < num_points; i++)
	{
		myMRH.Insert(points[i],i);
	}

	//Print all trees
	for( i = 0; i < 12; i++ ) 
	{
		cout << "\n\n######################\n";
		cout << "#### BASE CELL " << i << " ####\n";
		cout << "######################" << endl;
        myMRH.PrintTreeAtIndex(i); 
	}
}

void TestMRHInsertOneNodeQueryDisc
(
 int max_depth,
 double phi, 
 double theta, 
 int data, 
 double qPhi, 
 double qTheta, 
 double qRad
)
{
	vector<MortonNode> found;
	pointing pt;
	Healpix_Custom myHpxQ(1,NEST);
	int i,j;

	MultiResHpx myMRH = MultiResHpx(max_depth,NEST);

    // Insert the GIS longitude, latitude pairs into MRH
	pt.phi = phi;
	pt.theta = theta;
    cout << "Insert Node: " << phi << " " << theta << data << endl << endl;
    myMRH.Insert(pt,data);

	//Print all trees
	for( i = 0; i < 12; i++ ) 
	{
		cout << "\n\n######################\n";
		cout << "#### BASE CELL " << i << " ####\n";
		cout << "######################" << endl;
        myMRH.PrintTreeAtIndex(i); 
	}

	// Now perform disc query and output the results
	cout << endl << endl;
	cout << "################################\n";
	cout << "#### Now Perform Disc Query ####\n";
	cout << "################################\n";
	cout << "Query Center: longitude: " << qPhi
		 << " latitude: " << qTheta 
		 << " radius: " << qRad << endl << endl;
	pt.phi = qPhi;
	pt.theta = qTheta;
    found = myMRH.QueryDisc(pt,qRad);
	for(  j = 0; j < found.size(); j++) {
		cout << "Found: " << found[j].phi << " " << found[j].theta << " ";
		for( i = 0; i < found[j].data.size(); i++) cout << found[j].data[i] << " " << endl;
		cout << endl;
	}
	if(found.size() == 0) {
		cout << "No Results Found!\n";
	}
}


void TestMRHInsertOneNodeQueryDiscInclusive
(
 int max_depth,
 double phi, 
 double theta, 
 int data, 
 double qPhi, 
 double qTheta, 
 double qRad,
 int fact
)
{
	vector<MortonNode> found;
	pointing pt;
	Healpix_Custom myHpxQ(1,NEST);

	MultiResHpx myMRH = MultiResHpx(max_depth,NEST);

    // Insert the GIS phi, latitude pairs into MRH
    cout << "Insert Node: " << phi << " " << theta << data << endl << endl;
	pt.phi = phi;
	pt.theta = theta;
    myMRH.Insert(pt,data);
    int i;

	//Print all trees
	for( i = 0; i < 12; i++ ) 
	{
		cout << "\n\n######################\n";
		cout << "#### BASE CELL " << i << " ####\n";
		cout << "######################" << endl;
        myMRH.PrintTreeAtIndex(i); 
	}


}


void TestMRHInsertOneNodeQueryPolygon
(
 int max_depth,
 double phi,
 double theta,
 int data,
 int num_points,
 char* argv[]
)
{
	vector<MortonNode> found;
	Healpix_Custom myHpxQ(1,NEST);
    vector<pointing> poly;
	pointing pt;
	MultiResHpx myMRH = MultiResHpx(max_depth,NEST);

    // Insert the GIS longitude, latitude pairs into MRH
    cout << "Insert Node: " << phi << " " << theta << data << endl << endl;
	pt.phi = phi;
	pt.theta = theta;
    myMRH.Insert(pt,data);
    int i,j,index;

	//Print all trees
	for( i = 0; i < 12; i++ ) 
	{
		cout << "\n\n######################\n";
		cout << "#### BASE CELL " << i << " ####\n";
		cout << "######################" << endl;
        myMRH.PrintTreeAtIndex(i); 
	}
   	

	// Next perform polygon query and output the results
	cout << endl << endl;
	cout << "####################################\n";
	cout << "#### Now Perform Polygon Query  ####\n";
	cout << "####################################\n";
	// Get query polygon definition 
	index = 0;
	poly.clear();
	for( i = 0; i < num_points; i++ )
	{
     	pt.phi = atof(argv[7+index]);
     	pt.theta  = atof(argv[7+index+1]);
		poly.push_back(pt);
     	index += 2;
     	if(MRHDEBUG)
		{
     	   cout << "Next Polygon Point: " 
			    << pt.phi << " "
				<< pt.theta  << " " << endl;
		}
	}
    found = myMRH.QueryPolygon(poly);
	for(  j = 0; j < found.size(); j++) {
		cout << "Found: " << found[j].phi << " " << found[j].theta << " ";
		for( i = 0; i < found[j].data.size(); i++) cout << found[j].data[i] << " " << endl;
		cout << endl;
	}
	if(found.size() == 0) {
		cout << "No Results Found!\n";
	}
}

void TestMRHInsertOneNodeQueryPolygonInclusive
(
 int max_depth,
 double phi,
 double theta,
 int data,
 int num_points,
 int fact,
 char* argv[]
)
{
	vector<MortonNode> found;
	Healpix_Custom myHpxQ(1,NEST);
    vector<pointing> poly;
	pointing pt;
	MultiResHpx myMRH = MultiResHpx(max_depth,NEST);

    // Insert the GIS longitude, latitude pairs into MRH
    cout << "Insert Node: " << phi << " " << theta << data << endl << endl;
	pt.phi = phi;
	pt.theta = theta;
    myMRH.Insert(pt,data);
    int i,index;

	//Print all trees
	for( i = 0; i < 12; i++ ) 
	{
		cout << "\n\n######################\n";
		cout << "#### BASE CELL " << i << " ####\n";
		cout << "######################" << endl;
        myMRH.PrintTreeAtIndex(i); 
	}

	// Now get query polygon definition 
	index = 0;
	poly.clear();
	for( i = 0; i < num_points; i++ )
	{
     	pt.phi = atof(argv[8+index]);
     	pt.theta  = atof(argv[8+index+1]);
		poly.push_back(pt);
     	index += 2;
     	if(MRHDEBUG)
		{
     	   cout << "Next Polygon Point: " 
			    << pt.phi << " "
				<< pt.theta  << " " << endl;
		}
	}     	


}

void TestMRHInsertOneNodeQueryStrip
(
 int max_depth,
 double phi,
 double theta,
 int data,
 double theta1,
 double theta2
)
{
	vector<MortonNode> found;
	Healpix_Custom myHpxQ(1,NEST);
	MultiResHpx myMRH = MultiResHpx(max_depth,NEST);
	pointing pt;

    // Insert the GIS longitude, latitude pairs into MRH
    cout << "Insert Node: " << phi << " " << phi << data << endl << endl;
	pt.phi = phi;
	pt.theta = theta;
    myMRH.Insert(pt,data);
    int i,j;

	//Print all trees
	for( i = 0; i < 12; i++ ) 
	{
		cout << "\n\n######################\n";
		cout << "#### BASE CELL " << i << " ####\n";
		cout << "######################" << endl;
        myMRH.PrintTreeAtIndex(i); 
	}

	// Next perform query strip and output the results
	cout << endl << endl;
	cout << "#################################\n";
	cout << "#### Now Perform Strip Query ####\n";
	cout << "#################################\n";
    found = myMRH.QueryStrip(theta1,theta2);
	for(  j = 0; j < found.size(); j++) {
		cout << "Found: " << found[j].phi << " " << found[j].theta << " ";
		for( i = 0; i < found[j].data.size(); i++) cout << found[j].data[i] << " " << endl;
		cout << endl;
	}
	if(found.size() == 0) {
		cout << "No Results Found!\n";
	}
}

void TestMRHNeighbors
(
 int max_depth,
 double phi,
 double theta,
 int data,
 double qPhi,
 double qTheta
)
{
	vector<MortonNode> found;
	int64 order;
	pointing pt;
	Healpix_Custom myHpxQ(1,NEST);
	MultiResHpx myMRH = MultiResHpx(max_depth,NEST);

   // Insert the GIS longitude, latitude pair into MRH
   cout << "Insert Node: " << phi << " " << theta << data << endl << endl;
   pt.phi = phi;
   pt.theta = theta;
   myMRH.Insert(pt,data);
   int i,j;

	//Print all trees
	for( i = 0; i < 12; i++ ) 
	{
		cout << "\n\n######################\n";
		cout << "#### BASE CELL " << i << " ####\n";
		cout << "######################" << endl;
      myMRH.PrintTreeAtIndex(i); 
	}

	// Next perform neighbor query and output the results
	cout << endl << endl;
	cout << "#####################################\n";
	cout << "#### Now Perform Neighbors Query ####\n";
	cout << "#####################################\n";
    found = myMRH.Neighbors(pt,order);
	for(  j = 0; j < found.size(); j++) {
		cout << "Found: " << found[j].phi << " " << found[j].theta << " ";
		for( i = 0; i < found[j].data.size(); i++) cout << found[j].data[i] << " " << endl;
		cout << endl;
	}
	if(found.size() == 0) {
		cout << "No Results Found!\n";
	}
}


void TestHPXIndexToBase0Normalization(int order,int numPoints)
{
	double LO_THETA = -89.9*D2R; double LO_PHI = -179.9*D2R;
	double HI_THETA = 89.9*D2R; double HI_PHI = 179.9*D2R;
	double Q_LO_THETA = -40*D2R; double Q_LO_PHI = -179.9*D2R;
	double Q_HI_THETA = 40*D2R; double Q_HI_PHI = 179.9*D2R;
	pointing ptGIS,ptHPX,ptHPXCent,ptHPXNorm,ptHPXDeNorm;
	int fn;
    pair<int64,int> hpxId,hpxIDNorm,hpxIDDeNorm;
	Morton m;
	int64 HpxFromMorton;
	int OrderFromMorton;

	// Create MRH data structure 
	MultiResHpx mMRH = MultiResHpx(order,NEST);

	// Create HPX Query object
	Healpix_Custom hpxQ(order,NEST);

	// Generate random GIS longitude, latitude pairs
 	srand(time(NULL));
	printf("Order: %d\n",order);
	printf("GIS_Long,GIS_Lat,HPX_P,HPX_Z,HPX_Cent_P,HPX_Cent_Z,HPX_Norm_P,HPX_Norm_Z,HPX_DNorm_P,HPX_DNorm_Z,HPXID,Face_Num,NPFACE,HPXID_Shift,HPXID_Morton\n");
	for( int i = 0; i < numPoints; i++)	{
		ptGIS.phi = LO_PHI + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(HI_PHI-LO_PHI)));
		ptGIS.theta = LO_THETA + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(HI_THETA-LO_THETA)));
		//Convert GIS to HPX coordinates
		ptHPX = GIStoHPX(ptGIS); 

		//Compute HPX coordinates and get coordinates of cell center
	    hpxId.first = hpxQ.ang2pix(ptHPX);
	    hpxId.second = order;
		ptHPXCent = hpxQ.pix2ang (hpxId.first);

		//Compute Morton Node
		m = HpxToMorton(hpxId.first,hpxId.second);

		//Re-Compute HPX index from Morton Node
		MortonToHpx(m,HpxFromMorton,OrderFromMorton);

		// Normalize the HEALPix index to Base 0 indexing through shift
		// and compute cell center angle
		fn = hpxQ.FaceNum(ptHPX);
		hpxIDNorm.first = hpxId.first - fn*order_to_npface[order];
		hpxIDNorm.second = order;
		ptHPXNorm = hpxQ.pix2ang (hpxIDNorm.first);

		// Now De-Normalize the point (shift back the other way)
		hpxIDDeNorm.first = hpxIDNorm.first + fn*order_to_npface[order];
		hpxIDDeNorm.second = order;
		ptHPXDeNorm = hpxQ.pix2ang (hpxIDDeNorm.first);

		printf("%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%d,%d,%d,%d,%d\n",
				ptGIS.phi,ptGIS.theta,
				ptHPX.phi/pi,cos(ptHPX.theta),
				ptHPXCent.phi/pi,cos(ptHPXCent.theta),
				ptHPXNorm.phi/pi,cos(ptHPXNorm.theta),
				ptHPXDeNorm.phi/pi,cos(ptHPXDeNorm.theta),
				hpxId.first,fn,order_to_npface[order],
				hpxIDNorm.first,
				HpxFromMorton);
	}
}

void OutputHPXCellBoundaries(int64 hpxid,int order,int step)
{
	pair<int64,int> hpxidx;
	hpxidx.first = hpxid;
	hpxidx.second = order;
	HpxPZCellBoundaries(hpxidx,step);
}

void TestRadialSeperation(double phi1, double theta1, double phi2, double theta2)
{
	pointing p1,p2;
	p1.phi = phi1; p1.theta = theta1;
	p2.phi = phi2; p2.theta = theta2;
	cout << RadialDist(p1,p2) << endl;
}


void TestHPXInstantiation()
{
	int64 nside_,npface_,npix_;
	cout << "#### Test Creating HPX Custom of various orders ####\n";
	for(int i = 0; i <= 29; i++) {
		// Create HPX Query object
		Healpix_Custom hpxQ(i,NEST);
		cout << "	Successfully created HPX Custom order = " << i << endl;
	}

	cout << "\n#### Now Test Creating array of various sizes ####\n";
	//for(int64 i = 0; i <= 29; i++) {
	//	nside_ = int64(1) << i;
	//	npface_ = nside_ << i;
	//	npix_   = 12*npface_;

	//	// Create HPX Map object
	//	cout << "	Attempting to create array of size = " << npix_ << endl;
	//	HPXNode* newarray = new int[npix_];
	//}


	cout << "\n#### Now Test Creating HPX Map of various orders ####\n";
	for(int64 i = 0; i <= 29; i++) {
		nside_ = int64(1) << i;
		npface_ = nside_ << i;
		npix_   = 12*npface_;

		// Create HPX Map object
		cout << "	Attempting to create HPX Map order = " << i << " npix = " << npix_ << endl;
		Healpix_Map<Measurement> mHPX(i,NEST);
	}

}


int main(int64 argc, char* argv[])

{
	int testNum, order,max_depth,data;
	double longitude,latitude;
	Morton morton;
	int64 hpxIdx;
	std::string filename;

	testNum = atoi(argv[1]);


//####
//#### MULTIRESHPX TESTS
//####

	if(testNum == -1) {
		TestHPXInstantiation();
	}

	if(testNum == 0) {
		double phi1 = atof(argv[2]);
		double theta1 = atof(argv[3]);
		double phi2 = atof(argv[4]);
		double theta2 = atof(argv[5]);
		TestRadialSeperation(phi1,theta1,phi2,theta2);
	}

	if(testNum == 1) {
		cout << "Set HPX Map to different orders and sizes:" << endl;
		TestBuildHpx();
	}

	// Insert Data Index into Base 0 Cell
	if(testNum == 2) {
  	   max_depth = atoi(argv[2]);
       TestMRHInsertNodeBase0(max_depth);
	}

	// Insert Data Index into any Base Cell and Search for it
	if(testNum == 3) {
  	   max_depth = atoi(argv[2]);
	   longitude = atof(argv[3]);
	   latitude = atof(argv[4]);
       data = atoi(argv[5]);
       hpxIdx = StringToHpx(argv[6]);
	   order = atoi(argv[7]);
       TestMRHInsertOneNodeSearchOneNode(max_depth,longitude,latitude,data,hpxIdx,order);
	}

	// Insert N Data Indices from any Base Cell
	if(testNum == 4) {
  	   max_depth = atoi(argv[2]);
	   int seed = atoi(argv[3]);
       int num_points = atoi(argv[4]);
	   TestMRHInsertNodeAllBase(max_depth,seed,num_points);
	}

	// Insert N Data Indices from any Base Cell then Delete them all, one by one.
	if(testNum == 5) {
  	   max_depth = atoi(argv[2]);
	   int seed = atoi(argv[3]);
       int num_points = atoi(argv[4]);
	   TestMRHInsertNodeAllBaseThenDeleteAll(max_depth,seed,num_points);
	}

	// Insert N Data Indices from any Base Cell then do Disc Query and output results
	if(testNum == 6) {
  	   max_depth = atoi(argv[2]);
	   int seed = atoi(argv[3]);
       int num_points = atoi(argv[4]);
	   double longitude = atof(argv[5]);
	   double latitude = atof(argv[6]);
	   double radius = atof(argv[7]);
  	   TestMRHQueryDiscRandom(max_depth,seed,num_points,longitude,latitude,radius);
	}

	// Insert N Data Indices from any Base Cell then do Disc Query Inclusive and output results
	if(testNum == 7) {
  	   max_depth = atoi(argv[2]);
	   int seed = atoi(argv[3]);
       int num_points = atoi(argv[4]);
	   double longitude = atof(argv[5]);
	   double latitude = atof(argv[6]);
	   double radius = atof(argv[7]);
	   int fact = atoi(argv[8]);
  	   TestMRHQueryDiscInclusiveRandom(max_depth,seed,num_points,longitude,latitude,radius,fact);
	}


	// Insert One Data Index then do Disc Query and output results
	if(testNum == 8) {
  	   max_depth = atoi(argv[2]);
	   double longitude = atof(argv[3]);
	   double latitude = atof(argv[4]);
	   int data = atoi(argv[5]);
	   double qlong = atof(argv[6]);
	   double qlat = atof(argv[7]);
	   double radius = atof(argv[8]);
  	   TestMRHInsertOneNodeQueryDisc(max_depth,longitude,latitude,data,qlong,qlat,radius);
	}

	// Insert One Data Index then do Disc Query Inclusive and output results
	if(testNum == 9) {
  	   max_depth = atoi(argv[2]);
	   double longitude = atof(argv[3]);
	   double latitude = atof(argv[4]);
	   int data = atoi(argv[5]);
	   double qlong = atof(argv[6]);
	   double qlat = atof(argv[7]);
	   double radius = atof(argv[8]);
	   int fact = atoi(argv[9]);
  	   TestMRHInsertOneNodeQueryDiscInclusive(max_depth,longitude,latitude,data,qlong,qlat,radius,fact);
	}

	// Insert One Data Index then do Polygon Query and output results
	if(testNum == 10) {
  	   max_depth = atoi(argv[2]);
	   double longitude = atof(argv[3]);
	   double latitude = atof(argv[4]);
	   int data = atoi(argv[5]);
	   int numpoints = atoi(argv[6]);
  	   TestMRHInsertOneNodeQueryPolygon(max_depth,longitude,latitude,data,numpoints,argv);
	}


	// Insert One Data Index then do Polygon Inclusive Query and output results
	if(testNum == 11) {
  	   max_depth = atoi(argv[2]);
	   double longitude = atof(argv[3]);
	   double latitude = atof(argv[4]);
	   int data = atoi(argv[5]);
	   int numpoints = atoi(argv[6]);
	   int fact = atoi(argv[7]);
  	   TestMRHInsertOneNodeQueryPolygonInclusive(max_depth,longitude,latitude,data,numpoints,fact,argv);
	}

	// Insert One Data Index then do Strip Query and output results
	if(testNum == 12) {
  	   max_depth = atoi(argv[2]);
	   double longitude = atof(argv[3]);
	   double latitude = atof(argv[4]);
	   int data = atoi(argv[5]);
	   double lat1 = atof(argv[6]);
	   double lat2 = atof(argv[7]);
  	   TestMRHInsertOneNodeQueryStrip(max_depth,longitude,latitude,data,lat1,lat2);
	}


	// Insert One Data Index then do Neighbors Query and output results
	if(testNum == 13) {
  	   max_depth = atoi(argv[2]);
	   double longitude = atof(argv[3]);
	   double latitude = atof(argv[4]);
	   int data = atoi(argv[5]);
	   double qLong = atof(argv[6]);
	   double qLat = atof(argv[7]);
  	   TestMRHNeighbors(max_depth,longitude,latitude,data,qLong,qLat);
	}

	// Test of various HEALPix coordinate transformers and informers
	if(testNum == 14) {
  	   max_depth = atoi(argv[2]);
	   double longitude = atof(argv[3]);
	   double latitude = atof(argv[4]);
	   int order = atoi(argv[5]);
	   int step = atoi(argv[6]);
       TestMRHTransformersAndInfo(max_depth,longitude,latitude,order,step);
	}

	// Test HPX index to Base 0 HPX Normalization
	if(testNum == 15 ) {
  	   order = atoi(argv[2]);
       int numPoints = atoi(argv[3]);
	   TestHPXIndexToBase0Normalization(order,numPoints);
	}

	// Output the defining shape of HPXID at particular resolution
	if(testNum == 16) {
		int64 hpxid = atoi(argv[2]);
		order = atoi(argv[3]);
		int step = atoi(argv[4]);
		OutputHPXCellBoundaries(hpxid,order,step);
	}

}

