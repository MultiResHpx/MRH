#include "MultiResHpx.h"
// STD includes
#include <time.h>
#include <vector>
#include <math.h>
#include <cstdlib>
#include "lsconstants.h"
#include "healpix_map.h"

#define MORTONLQTDEBUG true


void TestMask()
{
   int max_order = 64;
	uint64 mask = 0;
	cout << sizeof(int64) << " " << sizeof(int) << endl;
	uint64 to_shift = 0;
	for(uint64 pos = 1; pos <= max_order; pos++)
	{
	   //shift_amount = ((2*i)-2);
      //   mask = 3 << shift_amount;
	   //cout << "Mask @ Order: " << i << " Shift Amount: " << shift_amount << " : " << mask << endl;
		if (pos == 1){ to_shift = 1 << pos; }
		 else if (pos == 2){ to_shift = 1 << pos;  }
		 else if (pos == 3){ to_shift = 1 << pos;   }
		 else if (pos == 4){ to_shift = 1 << pos;   }
		 else if (pos == 5){ to_shift = 1 << pos;   }
		 else if (pos == 6){ to_shift = 1 << pos;   }
		 else if (pos == 7){ to_shift = 1 << pos;   }
		 else if (pos == 8){ to_shift = 1 << pos;   }
		 else if (pos == 9){ to_shift = 1 << pos;   }
		 else if (pos == 10){ to_shift = 1 << pos;   }
		 else if (pos == 11){ to_shift = 1 << pos;   }
		 else if (pos == 12){ to_shift = 1 << pos;   }
		 else if (pos == 13){ to_shift = 1 << pos;   }
		 else if (pos == 14){ to_shift = 1 << pos;   }
		 else if (pos == 15){ to_shift = 1 << pos;   }
		 else if (pos == 16){ to_shift = 1 << pos;   }
		 else if (pos == 17){ to_shift = 1 << pos;   }
		 else if (pos == 18){ to_shift = 1 << pos;   }
		 else if (pos == 19){ to_shift = 1 << pos;   }
		 else if (pos == 20){ to_shift = 1 << pos;    }
		 else if (pos == 21){ to_shift = 1 << pos;    }
		 else if (pos == 22){ to_shift = 1 << pos;    }
		 else if (pos == 23){ to_shift = 1 << pos;    }
		 else if (pos == 24){ to_shift = 1 << pos;    }
		 else if (pos == 25){ to_shift = 1 << pos;    }
		 else if (pos == 26){ to_shift = 1 << pos;   }
		 else if (pos == 27){ to_shift = 1 << pos;   }
		 else if (pos == 28){ to_shift = 1 << pos;   }
		 else if (pos == 29){ to_shift = 1 << pos;  }
		 else if (pos == 30){ to_shift = 1 << pos;   }
		 else if (pos == 31){ to_shift = 1 << 1; to_shift = to_shift << 30;   }
		 else if (pos == 32){ to_shift = 1 << 1; to_shift = to_shift << 31;   }
		 else if (pos == 33){ to_shift = 1 << 1; to_shift = to_shift << 32;   }
		 else if (pos == 34){ to_shift = 1 << 2; to_shift = to_shift << 32;   }
		 else if (pos == 35){ to_shift = 1 << 3; to_shift = to_shift << 32;   }
		 else if (pos == 36){ to_shift = 1 << 4; to_shift = to_shift << 32;   }
		 else if (pos == 37){ to_shift = 1 << 5; to_shift = to_shift << 32;   }
		 else if (pos == 38){ to_shift = 1 << 6; to_shift = to_shift << 32;   }
		 else if (pos == 39){ to_shift = 1 << 7; to_shift = to_shift << 32;   }
		 else if (pos == 40){ to_shift = 1 << 8; to_shift = to_shift << 32;   }
		 else if (pos == 41){ to_shift = 1 << 9; to_shift = to_shift << 32;   }
		 else if (pos == 42){ to_shift = 1 << 10; to_shift = to_shift << 32;   }
		 else if (pos == 43){ to_shift = 1 << 11; to_shift = to_shift << 32;   }
		 else if (pos == 44){ to_shift = 1 << 12; to_shift = to_shift << 32;   }
		 else if (pos == 45){ to_shift = 1 << 13; to_shift = to_shift << 32;   }
		 else if (pos == 46){ to_shift = 1 << 14; to_shift = to_shift << 32;   }
		 else if (pos == 47){ to_shift = 1 << 15; to_shift = to_shift << 32;   }
		 else if (pos == 48){ to_shift = 1 << 16; to_shift = to_shift << 32;   }
		 else if (pos == 49){ to_shift = 1 << 17; to_shift = to_shift << 32;   }
		 else if (pos == 50){ to_shift = 1 << 18; to_shift = to_shift << 32;   }
		 else if (pos == 51){ to_shift = 1 << 19; to_shift = to_shift << 32;   }
		 else if (pos == 52){ to_shift = 1 << 20; to_shift = to_shift << 32;   }
		 else if (pos == 53){ to_shift = 1 << 21; to_shift = to_shift << 32;   }
		 else if (pos == 54){ to_shift = 1 << 22; to_shift = to_shift << 32;   }
		 else if (pos == 55){ to_shift = 1 << 23; to_shift = to_shift << 32;   }
	  	 else if (pos == 56){ to_shift = 1 << 24; to_shift = to_shift << 32;   }
		 else if (pos == 57){ to_shift = 1 << 25; to_shift = to_shift << 32;   }
		 else if (pos == 58){ to_shift = 1 << 26; to_shift = to_shift << 32;   }
		 else if (pos == 59){ to_shift = 1 << 27; to_shift = to_shift << 32;   }
		 else if (pos == 60){ to_shift = 1 << 28; to_shift = to_shift << 32;   }
		 else if (pos == 61){ to_shift = 1 << 29; to_shift = to_shift << 32;   }
		 else if (pos == 62){ to_shift = to_shift << 30; to_shift = to_shift << 32;   }
		 else if (pos == 63){ to_shift = to_shift << 31; to_shift = to_shift << 32;   }
		 else if (pos == 64){ to_shift = to_shift << 32; to_shift = to_shift << 32;   }
       //mask = 3 << to_shift;
	    cout << "Mask @ Order: " << pos << " Mask: " << to_shift << endl;
	
	}
}



//
// TEST HpxToMorton, MortonToHpx
//
void TestHpxToMorton(int64 hpxID, int order)
{
  //Inputs
  // arg[1] = HPX index
  // arg[2] = HPX order 
  // 2-bit address to Morton table
  //
  // MSB  LSB  Morton
  //============
  // 0    0   1
  // 0    1   2
  // 1    0   3
  // 1    1   4
  cout << "#### TEST #1 FROM COMMAND LINE ####" << endl;
  cout << "Inputs:" << endl;
  cout << "  HPX id: " << hpxID << endl;
  cout << "  Order: " << order << endl;
  cout << "Number 2-bit groups: " << order << endl;
  cout << "*** Convert HPX -> Morton ***" << endl;
  Morton morton = HpxToMorton(hpxID,order);                   
  cout << "Morton Code: "; PrintMorton(morton); cout << endl;
}

void TestMortonToHpx(Morton morton)
{
  int64 hpxId;
  int order;
  cout << "#### TEST #2 FROM COMMAND LINE ####" << endl;
  cout << "Inputs:" << endl;
  cout << "  Morton: "; PrintMorton(morton); cout << endl;
  cout << "*** Convert Morton -> HPX ***" << endl;
  MortonToHpx(morton,hpxId,order);
       
  cout << "HPX id: " << hpxId << " order: " << order << endl;
}


void TestPhiThetaToMorton(double phi, double theta, int level)
{     
	pointing pt;
	pt.phi = phi;
	pt.theta = theta;
	cout << "Inputs: " << endl;
	cout << "   Phi: " << phi << endl;
	cout << "   Theta:  " << theta << endl;
	cout << "   Morton Level: " << level << endl;

	Morton morton = PhiThetaToMorton(pt.phi,pt.theta,level);
	cout << "==================" << endl;
	cout << "   Unpacked Morton Code: "; PrintMorton(morton); cout << endl;
	cout << "SizeOfMorton: " << sizeof(morton);
}

void TestMortonToPhiTheta(Morton morton)
{    
	pointing pt;
	cout << "Inputs:" << endl;
	cout << "   Morton: "; PrintMorton(morton); cout << endl;

	pt = MortonToPhiTheta(morton);
	cout << "==================" << endl;
	cout << "   Phi: " << pt.phi << endl;
	cout << "   Theta: " << pt.theta << endl;
}

//
// TEST HpxToMorton, MortonToHpx
//



void TestAllNSIDE8HpxToMorton(int order)
{
  Morton morton;
  int npixel;
  int64 hpxId;
  cout << "#### TEST #3 Convert all normalized Order= " << order << " HPXids to Morton####" << endl;
  npixel = (1<<order)*(1<<order);
  for(int i = 0; i < npixel; i++ )
  {
     morton = HpxToMorton(i,order);
     MortonToHpx(morton,hpxId,order);
     cout << "HPX id: " << i << " Morton: "; 
				PrintMorton(morton); cout << " back to HPX id: " 
				<< hpxId << " of order: " << order << endl;
     cout << endl;
  }
}


void TestParentOfMorton(Morton morton)
{
	Morton parent;
	int mBit;
	int level = GetMortonLevel(morton);
	parent = ParentOfMorton(morton);

	cout << "Morton code: "; PrintMorton(morton); cout << " level: " << level << endl;
	cout << "Bits of Morton code:" << endl;
	for(int i = 1; i < level+1; i++)
	{
		mBit = GetMortonBit(morton,i);
		cout << "Level " << i << " : " << mBit << endl;
	}
	cout << "  Parent: "; PrintMorton(parent); cout << endl;
}

void TestChildrenOfMorton(Morton morton)
{
	int mBit;
	std::vector<Morton> children;
	int level = GetMortonLevel(morton);
	cout << "Morton code: "; PrintMorton(morton); cout << " level: " << level << endl;
	cout << "Bits of Morton code:" << endl;
	for(int i = 1; i < level+1; i++)
	{
		mBit = GetMortonBit(morton,i);
		cout << "Level " << i << " : " << mBit << endl;
	}
	children = ChildrenOfMorton(morton);
	cout << "  Children: "; PrintMorton(children[0]); cout << endl;
	cout << "            "; PrintMorton(children[1]); cout << endl;
	cout << "            "; PrintMorton(children[2]); cout << endl;
	cout << "            "; PrintMorton(children[3]); cout << endl;
}

void TestSiblingsOfMorton(Morton morton)
{
	int mBit;
	std::vector<Morton > siblings;
	int level = GetMortonLevel(morton);
	cout << "Morton code: "; PrintMorton(morton); cout << " level: " << level << endl;
	cout << "Bits of Morton code:" << endl;
	for(int i = 1; i < level+1; i++)
	{
		mBit = GetMortonBit(morton,i);
		cout << "Level " << i << " : " << mBit << endl;
	}
   siblings = SiblingsOfMorton(morton);
	cout << "  Siblings: "; PrintMorton(siblings[0]); cout << endl;
								   PrintMorton(siblings[1]); cout << endl;
								   PrintMorton(siblings[2]); cout << endl;
}

void TestInsertNode(int maxdepth,int num_insert, char* argv[],int isRadians)
{
    MortonLQT myMLQT = MortonLQT(maxdepth);
	int index = 0;
	int64 nextDat;
	std::vector<int64> datList;
	pointing pt;
	if(MORTONLQTDEBUG)
	{
	 cout << "Number Nodes to insert: " << num_insert << endl;
	}
    for(int i = 0; i < num_insert; i++ )
	{
     	pt.phi = atof(argv[4+index]);
     	pt.theta  = atof(argv[4+index+1]);
     	nextDat  = atoi(argv[4+index+2]);

     	if(MORTONLQTDEBUG)
		{
     	   cout << "Insert New Node: " << pt.phi << " " << pt.theta << " " << nextDat << endl;
		}
		if( isRadians == 0 ) {
		   pt.phi *= D2R;
		   pt.theta *= D2R;
		}
     	index += 3;
		datList.push_back(nextDat);
		myMLQT.InsertMortonNodeAtPhiTheta(pt,datList);
		datList.clear();
	}
    myMLQT.PrintMortonTree(); 
}

void TestSearchNode(int maxdepth,int num_insert, char* argv[])
{
     MortonLQT myMLQT = MortonLQT(maxdepth);
	 std::vector<MortonNode> FoundNodes;
     
    unsigned int i,j,index = 0;
	 int64 nextDat;
	 std::vector<int64> datList;
	 Morton morton;
	 pointing pt;
     for( i = 0; i < num_insert; i++ )
	 {
     	pt.phi = atof(argv[4+index]);
     	pt.theta  = atof(argv[4+index+1]);
     	nextDat  = atoi(argv[4+index+2]);
     	index += 3;
     	if(MORTONLQTDEBUG)
		{
     	   cout << "Insert New Node:" << pt.phi << " " << pt.theta << " " << nextDat << endl;
		}
		datList.push_back(nextDat);
		myMLQT.InsertMortonNodeAtPhiTheta(pt,datList);
		datList.clear();
	 }     	    
     myMLQT.PrintMortonTree(); 
     
     // Now search for all data nodes that fall under query node
    SetMortonBit(morton,0,1);
	 cout << endl << "Now search Morton LQT..." << endl;
     FoundNodes = myMLQT.SearchMortonNode(morton,1,0)  ;  
     cout << "Found the following data nodes under Morton Node: "; PrintMorton(morton); cout << endl;
     for( i = 0; i < FoundNodes.size(); i ++ )
	 {
        PrintMorton(FoundNodes[i].m); cout << " " 
			    << FoundNodes[i].sub << " ";
				for( j = 0; j < FoundNodes[i].data.size(); j++ ) cout << FoundNodes[i].data[j] << " ";
				cout << endl;
	 }

    SetMortonBit(morton,1,1);
     FoundNodes = myMLQT.SearchMortonNode(morton,1,0)  ;  
     cout << "Found the following data nodes under Morton Node: "; PrintMorton(morton); cout << endl;
     for( i = 0; i < FoundNodes.size(); i ++ )
	 {
        PrintMorton(FoundNodes[i].m); cout << " " 
			    << FoundNodes[i].sub << " ";
				for( j = 0; j < FoundNodes[i].data.size(); j++ ) cout << FoundNodes[i].data[j] << " ";
				cout << endl;
	 }

    SetMortonBit(morton,2,1);
     FoundNodes = myMLQT.SearchMortonNode(morton,1,0)  ;  
     cout << "Found the following data nodes under Morton Node: "; PrintMorton(morton); cout << endl;
     for( i = 0; i < FoundNodes.size(); i ++ )
	 {
        PrintMorton(FoundNodes[i].m); cout << " " 
			    << FoundNodes[i].sub << " ";
				for( j = 0; j < FoundNodes[i].data.size(); j++ ) cout << FoundNodes[i].data[j] << " ";
				cout << endl;
	 }

    SetMortonBit(morton,3,1);
     FoundNodes = myMLQT.SearchMortonNode(morton,1,0)  ;  
     cout << "Found the following data nodes under Morton Node: "; PrintMorton(morton); cout << endl;
     for( i = 0; i < FoundNodes.size(); i ++ )
	 {
        PrintMorton(FoundNodes[i].m); cout << " " 
			    << FoundNodes[i].sub << " ";
				for( j = 0; j < FoundNodes[i].data.size(); j++ ) cout << FoundNodes[i].data[j] << " ";
				cout << endl;  
	 }

  //   morton = -1;
  //   FoundNodes = myMLQT.SearchMortonNode(morton,1)  ;  
  //   cout << "Found the following data nodes under root Morton Node: " << morton << endl;
  //   for( i = 0; i < FoundNodes.size(); i ++ )
	 //{
  //      PrintMorton(FoundNodes[i].morton); cout << " " 
		//	    << FoundNodes[i].sub << " " 
		//		 << FoundNodes[i].data << endl;
	 //}
}


void TestDeleteNode(int maxdepth,int num_insert, char* argv[])
{
     MortonLQT myMLQT = MortonLQT(maxdepth);
	 std::vector<pointing> points;
	 std::vector<int> indices;
	 pointing pt;
	 Morton morton;
     
     unsigned int i,index = 0;
	 int64 nextDat;
	 std::vector<int64> datList;
     for( i = 0; i < num_insert; i++ )
	 {
     	pt.phi = atof(argv[4+index]);
     	pt.theta  = atof(argv[4+index+1]);
     	nextDat  = atoi(argv[4+index+2]);
		indices.push_back(nextDat);
		points.push_back(pt);
     	index += 3;
     	if(MORTONLQTDEBUG)
		{
     	   cout << "Insert New Node: " << pt.phi << " "
												 << pt.theta << " "
												 << nextDat << endl;
		}
		datList.push_back(nextDat);
		myMLQT.InsertMortonNodeAtPhiTheta(pt,datList);
		datList.clear();
	 }     	    
     myMLQT.PrintMortonTree(); 

     // Now delete nodes one by one and print tree each time
	 cout << endl << endl << "### Now delete nodes one by one and print tree after each deletion.###" << endl << endl;
     for( i = 0; i < points.size(); i++ )
	 {
        morton = myMLQT.FindMortonAtPhiTheta(points[i]);
		cout << "Next Morton to delete: "; PrintMorton(morton); cout << endl;
        if(MORTONLQTDEBUG)
		{
           cout << endl << endl << "Delete Node at: phi: " << points[i].phi 
			                       << " theta " << points[i].theta
										  << " data " << indices[i]
										  << " "; PrintMorton(morton); cout << endl;
		}
		myMLQT.DeleteMortonNodeAtPhiTheta(points[i]);
        myMLQT.PrintMortonTree(); 
	 }
}

void TestWriteTree(int maxdepth,int num_insert, std::string filename,  char* argv[])
{
     MortonLQT myMLQT = MortonLQT(maxdepth);
	 std::vector<pointing> points;
	 pointing pt;
	 std::vector<int> indices;
	 std::vector<int> morton_sub;
     
     unsigned int i,index = 0;
	 int64 nextDat;
	 std::vector<int64> datList;
     for( i = 0; i < num_insert; i++ )
	 {
     	pt.phi = atof(argv[5+index]);
     	pt.theta  = atof(argv[5+index+1]);
     	nextDat  = atoi(argv[5+index+2]);
		points.push_back(pt);
		indices.push_back(nextDat);
     	index += 3;
     	if(MORTONLQTDEBUG)
		{
     	   cout << "Insert New Node: " << pt.phi << " " 
				                         << pt.theta << " "
												 << nextDat << endl;
		}
		datList.push_back(nextDat);
		myMLQT.InsertMortonNodeAtPhiTheta(pt,datList);
		datList.clear();
	 }     	    
     myMLQT.PrintMortonTree();      

     // Write out MLQT to file
	 cout << "Now output MortonLQT to file: " << filename.c_str() << endl;
	ofstream fp;
	fp.open(filename.c_str() );
	if(fp.fail())
	{
	  cout << "Unable to open MortonLQT output file: " << filename.c_str() << endl;
	  exit(1);
	}
    myMLQT.SaveTreeToFile(fp);
	fp.close();
}

void TestLoadTree(int maxdepth,std::string filename)
{
     Morton morton;
	 unsigned int i,j;
	 std::vector<MortonNode> FoundNodes;
	 MortonLQT myMLQT = MortonLQT(maxdepth);

     // Read a MLQT from file
	 cout << endl << "### Load Morton Linear Quadtree from " << filename.c_str() << " ###" << endl << endl;
	 ifstream fp;
	 fp.open(filename.c_str());

	 if(fp.fail())
	 {
	   cout << "Unable to open MortonLQT input file: " << filename.c_str() << endl;
	   exit(1);     
	 }

     myMLQT.LoadTreeFromFile(fp);    

	 fp.close();
     
     // Print the tree for file output verification
     myMLQT.PrintMortonTree();

     // Now search for all data nodes that fall under query node
	  SetMortonBit(morton,1,1);
     FoundNodes = myMLQT.SearchMortonNode(morton,-1,0)  ;  
     cout << "Found the following data nodes under Morton Node: "; PrintMorton(morton); cout << endl;
     for( i = 0; i < FoundNodes.size(); i ++ )
	 {
        PrintMorton(FoundNodes[i].m); cout << " "
			    << FoundNodes[i].sub << " ";
				for( j = 0; j < FoundNodes[i].data.size(); j++ ) cout << FoundNodes[i].data[j] << " ";
				cout << endl;
	 }

	  SetMortonBit(morton,2,1);
     FoundNodes = myMLQT.SearchMortonNode(morton,-1,0)  ;  
     cout << "Found the following data nodes under Morton Node:"; PrintMorton(morton); cout << endl;
     for( i = 0; i < FoundNodes.size(); i ++ )
	 {
        PrintMorton(FoundNodes[i].m); cout << " "
			    << FoundNodes[i].sub << " ";
				for( j = 0; j < FoundNodes[i].data.size(); j++ ) cout << FoundNodes[i].data[j] << " ";
				cout << endl;   
	 }

	  SetMortonBit(morton,3,1);
     FoundNodes = myMLQT.SearchMortonNode(morton,-1,0)  ;  
     cout << "Found the following data nodes under Morton Node: "; PrintMorton(morton); cout << endl;
     for( i = 0; i < FoundNodes.size(); i ++ )
	 {
        PrintMorton(FoundNodes[i].m); cout << " "
			    << FoundNodes[i].sub << " ";
				for( j = 0; j < FoundNodes[i].data.size(); j++ ) cout << FoundNodes[i].data[j] << " ";
				cout << endl; 
	 }

	  SetMortonBit(morton,4,1);
     FoundNodes = myMLQT.SearchMortonNode(morton,-1,0)  ;  
     cout << "Found the following data nodes under Morton Node: "; PrintMorton(morton); cout << endl;
     for( i = 0; i < FoundNodes.size(); i ++ )
	 {
        PrintMorton(FoundNodes[i].m); cout << " "
			    << FoundNodes[i].sub << " ";
				for( j = 0; j < FoundNodes[i].data.size(); j++ ) cout << FoundNodes[i].data[j] << " ";
				cout << endl; 
	 }


}

    
int main(int64 argc, char* argv[])

{
	int testNum, order, level,max_depth,num_insert;
	double longitude,latitude;
	Morton morton;
	int64 hpxIdx;
	std::string filename;

	testNum = atoi(argv[1]);

	if(MORTONLQTDEBUG) {
		cout << "Test #" << testNum << endl;
	}

//####
//#### MORTONLQT TESTS
//####


	if(testNum == 1) {
		morton = StringToMorton(argv[2]);
		cout << "Morton: "; PrintMorton(morton); cout << endl;
		cout << "Morton Str: " << MortonToString(morton) << endl;
		cout << "Size of Morton: " << sizeof(morton) << endl;

		uint64 mWord = GetMortonWord(morton);

		cout << "Morton Word: " << mWord << endl;
	}

	if(testNum == 2) {
		//order = atoi(argv[2]);
		//int face_num = atoi(argv[3]);
		//int step = atoi(argv[4]);
  //      OutputBaseCellBoundary(order,face_num,step);
		TestMask();
	}

    // HEALPix to Morton
	if(testNum == 3 ) {
		hpxIdx = StringToHpx(argv[2]);
		order = atoi(argv[3]);
      TestHpxToMorton(hpxIdx,order);
	}
	
	// Morton to HEALPix
	if(testNum == 4 ) {
		morton = StringToMorton(argv[2]);
        TestMortonToHpx(morton);
	}

	// All NSIDE HEALPix to Morton
	if(testNum == 5 ) {
		order = atoi(argv[2]);
        TestAllNSIDE8HpxToMorton(order);
	}

	// Longitude, Latitude to Morton
	if(testNum == 6 ) {
	   longitude = atof(argv[2]);
	   latitude = atof(argv[3]);
	   level = atoi(argv[4]);
       TestPhiThetaToMorton(longitude,latitude,level);
	}

	// Morton to Longitude, Latitude
	if(testNum == 7 ) {
	   morton = StringToMorton(argv[2]);
       TestMortonToPhiTheta(morton);
	}

	// Parent of Morton
	if(testNum == 8) {
	   morton = StringToMorton(argv[2]);
       TestParentOfMorton(morton);
	}

	// Children of Morton
	if(testNum == 9) {
	   morton = StringToMorton(argv[2]);
       TestChildrenOfMorton(morton);
	}

	// Siblings of Morton
	if(testNum == 10) {
	   morton = StringToMorton(argv[2]);
       TestSiblingsOfMorton(morton);
	}

	// Insert Morton Node
	if(testNum == 11) {
	   max_depth = atoi(argv[2]);
	   num_insert = atoi(argv[3]);
      TestInsertNode(max_depth,num_insert,argv,1);
	}

	// Insert Morton Node DEGREES
	if(testNum == 12) {
	   max_depth = atoi(argv[2]);
	   num_insert = atoi(argv[3]);
       TestInsertNode(max_depth,num_insert,argv,0);
	}

	// Search Morton Node
	if(testNum == 13) {
	   max_depth = atoi(argv[2]);
	   num_insert = atoi(argv[3]);
       TestSearchNode(max_depth,num_insert,argv);
	}

	// Delete Morton Node
	if(testNum == 14) {
	   max_depth = atoi(argv[2]);
	   num_insert = atoi(argv[3]);
      TestDeleteNode(max_depth,num_insert,argv);
	}

	// Write Morton LQT
	if(testNum == 15) {
	   max_depth = atoi(argv[2]);
	   num_insert = atoi(argv[3]);
	   filename = argv[4];
      TestWriteTree(max_depth,num_insert,filename,argv);
	}

	// Load Morton LQT
	if(testNum == 16) {
	   max_depth = atoi(argv[2]);
	   filename = argv[3];
      TestLoadTree(max_depth,filename);
	}
}