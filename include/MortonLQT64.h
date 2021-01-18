#ifndef __MORTONLQT64_H__
#define __MORTONLQT64_H__

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <math.h>
#include <limits>
#include <cmath>
#include <iomanip>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

//HEALPix Library
#include "healpix_custom.h"
#include "MultiResHpx_tables.h"
#include "geom_utils.h"

//using std::vector;
#include "lsconstants.h"
using namespace::std;

//GLOBALS
//#define MAXDEPTH 13 //32BIT machines
#define MAXDEPTH 29 //64BIT machines
//#define SIBLINGNODES 
//#define MLQT_VERBOSE 
#define R2D 57.29577951308
#define D2R  0.01745329252

#define WORDMETHOD
//#define CRITCOUNT
//#define BENCHMARKING

const double cos45 = cos(degr2rad*45.0);

// HPX Resolution Table
const double order_to_cellres[] = {
	1.023326707946480, 0.511663353973244, 0.255831676986622, 0.127915838493311,
	0.063957919246656, 0.031978959623328, 0.015989479811664, 0.007994739905832,
	0.003997369952916, 0.001998684976458, 0.000999342488229, 0.000499671244114,
	0.000249835622057, 0.000124917811029, 0.000062458905514, 0.000031229452757,
	0.000015614726379, 0.000007807363189, 0.000003903681595, 0.000001951840797,
	0.000000975920399, 0.000000487960199, 0.000000243980100, 0.000000121990050,
	0.000000060995025, 0.000000030497512, 0.000000015248756, 0.000000007624378,
	0.000000003812189, 0.000000001906095
};

const int64 order_to_npface[] = { 1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576,
 4194304, 16777216, 67108864, 268435456, 1073741824, 4294967296, 17179869184, 68719476736, 274877906944,
 1099511627776, 4398046511104, 17592186044416, 70368744177664, 281474976710656, 1125899906842624,
 4503599627370496, 18014398509481984, 72057594037927936, 288230376151711744 };

const double basex_to_base0[13][2] = { 
	{0.0,	 0.0}, //Base 0 to Base 0
	{-0.5*pi,0.0}, //Base 1 to Base 0
	{-pi,	 0.0}, //Base 2 to Base 0
	{-1.5*pi,0.0}, //Base 3 to Base 0
	
	{-1.75*pi,-0.7297723785}, //Base 4b to Base 0 (Long > 315 and Long < 360)
	{0.25*pi, -0.7297723785},   //Base 4a to Base 0 (Long > 0 and Long < 45)
	{-0.25*pi,-0.7297723785}, //Base 5 to Base 0
	{-0.75*pi,-0.7297723785}, //Base 6 to Base 0
	{-1.25*pi,-0.7297723785}, //Base 7 to Base 0
	
	{0.0,	 -1.4595447570}, //Base 8 to Base 0
	{-0.5*pi,-1.4595447570}, //Base 9 to Base 0
	{-pi,	 -1.4595447570},  //Base 10 to Base 0
	{-1.5*pi,-1.4595447570}}; //Base 11 to Base 0
   
struct Morton {
public:
	Morton() :	LEVEL(0),
				L1(0),L2(0),L3(0),L4(0),L5(0),L6(0),
				L7(0),L8(0),L9(0),L10(0),L11(0),L12(0),
				L13(0),L14(0),L15(0),L16(0),L17(0),L18(0),
				L19(0),L20(0),L21(0),L22(0),L23(0),L24(0),
				L25(0),L26(0),L27(0),L28(0),L29(0)
	{	L1 = 0;	L2 = 0;	L3 = 0;	L4 = 0;	L5 = 0;	L6 = 0;
		L7 = 0;	L8 = 0;	L9 = 0;	L10 = 0;	L11 = 0;	L12 = 0;	
		L13 = 0;	L14 = 0;	L15 = 0;	L16 = 0;	L17 = 0;	L18 = 0;	
		L19 = 0;	L20 = 0;	L21 = 0;	L22 = 0;	L23 = 0;	L24 = 0;	
		L25 = 0;	L26 = 0;	L27 = 0;	L28 = 0;	L29 = 0;	
		
	}


			// two-bit unsigned short, allows values 0,1,2,3
			unsigned short L29 : 2; //Morton Level 1 Bit
			unsigned short L28 : 2; //Morton Level 2 Bit
			unsigned short L27 : 2; //Morton Level 3 Bit
			unsigned short L26 : 2; //Morton Level 4 Bit
			unsigned short L25 : 2; //Morton Level 5 Bit
			unsigned short L24 : 2; //Morton Level 6 Bit
			unsigned short L23 : 2; //Morton Level 7 Bit
			unsigned short L22 : 2; //Morton Level 8 Bit
			unsigned short L21 : 2; //Morton Level 9 Bit
			unsigned short L20 : 2; //Morton Level 10 Bit
			unsigned short L19 : 2; //Morton Level 11 Bit
			unsigned short L18 : 2; //Morton Level 12 Bit
			unsigned short L17 : 2; //Morton Level 13 Bit
			unsigned short L16 : 2; //Morton Level 14 Bit
			unsigned short L15 : 2; //Morton Level 15 Bit
			unsigned short L14 : 2; //Morton Level 16 Bit
			unsigned short L13 : 2; //Morton Level 17 Bit
			unsigned short L12 : 2; //Morton Level 18 Bit
			unsigned short L11 : 2; //Morton Level 19 Bit
			unsigned short L10 : 2; //Morton Level 20 Bit
			unsigned short L9 : 2; //Morton Level 21 Bit
			unsigned short L8 : 2; //Morton Level 22 Bit
			unsigned short L7 : 2; //Morton Level 23 Bit
			unsigned short L6 : 2; //Morton Level 24 Bit
			unsigned short L5 : 2; //Morton Level 25 Bit
			unsigned short L4 : 2; //Morton Level 26 Bit
			unsigned short L3 : 2; //Morton Level 27 Bit
			unsigned short L2 : 2; //Morton Level 28 Bit
			unsigned short L1 : 2; //Morton Level 29 Bit
			unsigned short LEVEL : 6; //Level of Morton Code

} ;

//Morton Data Type Utility Functions
Morton		StringToMorton(std::string mStr);
int64		StringToHpx(std::string mStr);
std::string MortonToString(Morton m);
void		SetMortonBit(Morton& m,int64 mBit,int bitLevel);
int64		GetMortonBit(Morton m,int bitLevel);
int			GetMortonLevel(Morton m);   
void		PrintMorton(Morton m);
uint64		GetMortonWord(Morton m);
int64		GetBitMask(int order);
Morton		AppendMortonBit(Morton m,int bit);	
int64		AppendHPXBit(int64 hpxid,int bit); 
std::vector<Morton> SiblingsOfMorton(Morton m);	
Morton		ParentOfMorton(Morton m);   	
std::vector<Morton> ChildrenOfMorton(Morton m);		
Morton		HpxToMorton(int64 hpxid,int order);	   
void		MortonToHpx(Morton m,int64& hpxid, int& order);	 	
Morton		PhiThetaToMorton(double phi,double theta,int level);	
pointing	MortonToPhiTheta(Morton m);
bool		Equals(Morton a, Morton b);
bool		GreaterThan(Morton a, Morton b);
bool		LessThan(Morton a, Morton b);
int         LessThanV2(Morton a, Morton b);
static		Healpix_Custom hpxQ;


class MortonNode 
{
public:
	MortonNode();
	MortonNode(Morton _m, int _sub, int _childrenYN, double _phi, double _theta, int _data);
	~MortonNode(){};

//private:
	Morton m;
	int sub;
	int childrenYN;
	double phi;
	double theta;
	std::vector<int64> data;
	int sentinel; 
};

class MortonLQT 
{
public:
	 int CRIT_COUNT;

// Methods
public:
	MortonLQT(){ userMaxDepth = 29; faceNum = 0;}
	MortonLQT(int depth){userMaxDepth = depth; faceNum = 0; CRIT_COUNT = 0; };
	//MortonLQT(int depth){userMaxDepth = depth; faceNum = 0; };
	~MortonLQT(){};

	int GetCritCount();
	void ResetCritCount();

	std::vector<int> GetMortonNodeSizeHistogram();

	int GetNumMortonNodes();
	
	int GetMortonTreeDepth();

	int GetMortonTreeMinDepth();

	double GetAvgMortonTreeDepth();
	 
	void SetFaceNum(int face_num);

	void ClearTree();

	//// FIND / SEARCH NODE
	
	int FindIndexAtMorton(Morton m);
	
	int FindIndexAtMortonSub(Morton m,int sub);
	
	std::vector<MortonNode> SearchMortonNodeAtPhiTheta(pointing pt,int sentinel);

	std::vector<MortonNode> SearchMortonNodeHpx(int64 hpxid,int order,int sentinel);
	
	std::vector<MortonNode> SearchMortonNode(Morton m,int sub,int sentinel);
		
	Morton FindMortonAtPhiTheta(pointing pt);

	Morton ClosestMortonToPhiTheta(std::vector<MortonNode> NodeList,pointing pt);

	//// GET, APPEND, ADD INSERT NODE 
	MortonNode GetNodeAtIndex(unsigned int index);

	//std::vector<MortonNode> GetClosestMortonNodeAtPhiTheta(pointing pt);

	std::vector<MortonNode> GetNodeAtMorton(Morton m,int sub);
	    
	void UpdateMortonNode(MortonNode node);

	void UpdateMortonNodeSentinel(MortonNode node);
	    
	void AppendMortonNode(MortonNode node);
		
	void AddMortonNode(MortonNode node);
	void AddMortonNode2(MortonNode node);
	
	Morton InsertMortonNodeAtMorton(Morton m,std::vector<int64> data);
	
	Morton InsertMortonNodeAtHpx(int64 hpxid,int order,std::vector<int64> data);
	
	Morton InsertMortonNodeAtPhiTheta(pointing pt,std::vector<int64> data);

	//// DELETE NODE

	void DeleteMortonNodeAtPhiTheta(pointing pt);
	
	void DeleteMortonNodeAtHpx(int64 hpxid,int order);
	
	void DeleteMortonNodeAtMorton(Morton m,int sub);

	//// WRITE TREE TO FILE
	
	void SaveTreeToFile(ofstream& fp);
	
	void LoadTreeFromFile(ifstream& fp);

	//// PRINT TREE
	void PrintMortonTree();
	void PrintMortonList();

					
private:
	void SearchMortonNode_internal(std::vector<MortonNode>& found_nodes,Morton m,int sub,int sentinel);
	Morton InsertMortonNode_internal_1(pointing pt,std::vector<int64> data);
	int FindIndexAtMortonSub_internal(Morton sMorton, int sSubs);
	int FindIndexAtMortonSub_BinarySearch_internal(Morton sMorton, int sSubs);
    void WriteHeader(ofstream& fp);
	void WriteData(ofstream& fp);


// Attributes
private:
	std::vector<MortonNode> mTree;
	int userMaxDepth;
	int curTreeDepth;
	int faceNum;
};


inline void skip_line(ifstream &IN)
{
  char next_char;	
  next_char = '\0'; //reset
  while( next_char != '\n' )
  {
    IN.get(next_char);
  }
}

inline bool SortFunctionMorton2( MortonNode a, MortonNode b );

inline bool SortFunctionMorton( MortonNode a, MortonNode b ) { return( LessThan(a.m,b.m) ); }

inline bool SortFunctionMapIdx( MortonNode a, MortonNode b ) { return( a.data[0] < b.data[0] ); }

inline int QSortFunctionMorton( const void* a, const void* b ); 

inline int MortonLQT::GetCritCount( ){ return( CRIT_COUNT); }

inline void MortonLQT::ResetCritCount(){ CRIT_COUNT = 0; }


#endif 