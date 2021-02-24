#ifndef __MORTONLQT_H__
#define __MORTONLQT_H__

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

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

//HEALPix Library
#include "healpix_custom.h"
#include "MultiResHpx_tables.h"

//using std::vector;
#include "lsconstants.h"
using namespace::std;

//GLOBALS
#define R2D 57.29577951308
#define D2R  0.01745329252
   
// Compute constants
const double NegSin45 = sin(-45.0*D2R);
const double NegCos45 = cos(-45.0*D2R);
const double initLeft = 0.0;
const double initRight = (sqrt(2.0)/4.0)*pi;
const double initTop = 0.0;
const double initBottom = (sqrt(2.0)/4.0)*pi;

class MortonNode 
{
public:
	MortonNode();
	MortonNode(int64 _morton, int _sub, int _childrenYN, float _longitude, float _latitude, int _data, int _facenum);
	~MortonNode(){};

//private:
	int64 morton;
	int sub;
	int childrenYN;
	float longitude;
	float latitude;
	int data;
	int facenum;
};


class MortonLQT 
{
// Methods
public:
	MortonLQT();
	MortonLQT(int depth);
	~MortonLQT();

	unsigned int GetNumMortonNodes();
	
	int GetMortonTreeDepth();
	
	void PrintMortonTree();
	
	int FindIndexAtMorton(int64 morton);
	
	int FindIndexAtMortonSub(int64 morton,int sub);
	
	std::vector<MortonNode> SearchMortonNodeAtLongLat(float longitude,float latitude);

	std::vector<MortonNode> SearchMortonNodeHpx(int64 hpxid,int order);
	
	std::vector<MortonNode> SearchMortonNode(int morton,int sub);
	
	std::vector<int64> FindMortonAtLongLat(float longitude,float latitude);
	
	std::vector<MortonNode> GetNodeAtMorton(int morton,int sub);
    
	void UpdateNode(MortonNode node);
    
	void AppendNode(MortonNode node);
	
	void AddNode(MortonNode node);

	int64 InsertMortonNodeAtMorton(int64 morton,int data);
	
	int64 InsertMortonNodeAtHpx(int64 hpxid,int order,int data);
	
	int64 InsertMortonNodeAtLongLat(float longitude,float latitude,int data);

	void DeleteMortonNodeAtLongLat(float longitude,float latitude);
	
	void DeleteMortonNodeAtHpx(int64 hpxid,int order);
	
	void DeleteMortonNodeAtMorton(int64 morton,int sub);
	
	int SaveTreeToFile(std::string filename);
	
	int LoadTreeFromFile(std::string filename);
	
	int GetMortonLevel(int64 morton);   
	
	int64 GetMortonBit(int64 morton,int bitLevel);
    	
	int64 GetBitMask(int order);
 
	int64 AppendMortonBit(int64 morton,int bit);    
	
	int64 AppendHPXBit(int64 hpxid,int bit);
    
	std::vector<int64> SiblingsOfMorton(int64 morton);   
	
	int64 ParentOfMorton(int64 morton);   
	
	std::vector<int64> ChildrenOfMorton(int64 morton);
	
	int64 HpxToMorton(int64 hpxid,int order);	  
	
	void MortonToHpx(int64 morton,int64& hpxid, int& order);	
		
	int64 LongLatToMorton(float longHpx,float latHpx,int level);
	
	int64 LongLatDegToMorton(float longHpx,float coLatHpx,int order);
		
	void MortonToLongLat(int64 morton,float& longHpx, float& coLatHpx);
		
	std::vector<int64> ClosestMortonToLongLat(std::vector<MortonNode> NodeList,float qLong,float qLat);

private:

	void SearchMortonNode_internal(std::vector<MortonNode>& found_nodes,int64 morton, int sub);
	int64 InsertNode_internal_1(float longitude,float latitude,int data);
	int FindIndexAtMortonSub_internal(int64 sMorton, int sSubs);
    void WriteHeader(FILE* fp);
	void WriteData(FILE* fp);

// Attributes
private:
	std::vector<MortonNode> mTree;
	int userMaxDepth;
	int curTreeDepth;
	Healpix_Custom hpxQ;
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

inline bool SortFunctionMorton( MortonNode a, MortonNode b ) { return( a.morton < b.morton ); }


#endif 