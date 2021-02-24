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

#include "MortonLQT64.h"


//#######
//####### MORTON DATA TYPE UTILITY FUNCTIONS 
//#######

int GetMortonLevel(Morton m)
{
	return m.LEVEL;
}

uint64 GetMortonWord(Morton m)
{
	uint64 MORTON = 0;

	// Get Morton Bits & Level
	memcpy(&MORTON,&m,sizeof(m));
	
	return MORTON;
}


void PrintMorton(Morton m)
{
   for(int i = 0; i < m.LEVEL; i++ )
   {
		cout << GetMortonBit(m,i+1)+1;
   }
}

bool GreaterThan(Morton a, Morton b)
{

	uint64 mA = GetMortonWord(a);
	uint64 mB = GetMortonWord(b);
	return mA > mB;
}


bool LessThan(Morton a, Morton b)
{
	uint64 mA = GetMortonWord(a);
	uint64 mB = GetMortonWord(b);
	return mA < mB;
}

bool Equals(Morton a, Morton b)
{
	uint64 mA = GetMortonWord(a);
	uint64 mB = GetMortonWord(b);
	return mA == mB;
}


bool AlmostEqual(double a,double b,int ulp) {
	return std::fabs(a-b) < std::numeric_limits<double>::epsilon()*std::abs(a+b)*ulp  || std::abs(a-b) < std::numeric_limits<double>::min();
}


int QSortFunctionMorton( const void* va, const void* vb )
{
	const MortonNode* a = static_cast<const MortonNode*>(va);
	const MortonNode* b = static_cast<const MortonNode*>(vb);
	//Morton a = na->m;
	//Morton b = nb->m;
	
	if( a->m.LEVEL < b->m.LEVEL ) { 
		return -1; 
	} else if( a->m.LEVEL > b->m.LEVEL ) {
		return 1;
	} else { // a and b are same level so now must compare each Morton Bit.
		// Compare bit level by bit level. Return once an equality is discovered.
		if ( a->m.L1 < b->m.L1 ) { return -1; } else if ( a->m.L1 > b->m.L1 ) { return 1; }
		if ( a->m.L2 < b->m.L2 ) { return -1; } else if ( a->m.L2 > b->m.L2 ) { return 1; }
		if ( a->m.L3 < b->m.L3 ) { return -1; } else if ( a->m.L3 > b->m.L3 ) { return 1; }
		if ( a->m.L4 < b->m.L4 ) { return -1; } else if ( a->m.L4 > b->m.L4 ) { return 1; }
		if ( a->m.L5 < b->m.L5 ) { return -1; } else if ( a->m.L5 > b->m.L5 ) { return 1; }
		if ( a->m.L6 < b->m.L6 ) { return -1; } else if ( a->m.L6 > b->m.L6 ) { return 1; }
		if ( a->m.L7 < b->m.L7 ) { return -1; } else if ( a->m.L7 > b->m.L7 ) { return 1; }
		if ( a->m.L8 < b->m.L8 ) { return -1; } else if ( a->m.L8 > b->m.L8 ) { return 1; }
		if ( a->m.L9 < b->m.L9 ) { return -1; } else if ( a->m.L9 > b->m.L9 ) { return 1; }
		if ( a->m.L10 < b->m.L10 ) { return -1; } else if ( a->m.L10 > b->m.L10 ) { return 1; }
		if ( a->m.L11 < b->m.L11 ) { return -1; } else if ( a->m.L11 > b->m.L11 ) { return 1; }
		if ( a->m.L12 < b->m.L12 ) { return -1; } else if ( a->m.L12 > b->m.L12 ) { return 1; }
		if ( a->m.L13 < b->m.L13 ) { return -1; } else if ( a->m.L13 > b->m.L13 ) { return 1; }
		if ( a->m.L14 < b->m.L14 ) { return -1; } else if ( a->m.L14 > b->m.L14 ) { return 1; }
		if ( a->m.L15 < b->m.L15 ) { return -1; } else if ( a->m.L15 > b->m.L15 ) { return 1; }
		if ( a->m.L16 < b->m.L16 ) { return -1; } else if ( a->m.L16 > b->m.L16 ) { return 1; }
		if ( a->m.L17 < b->m.L17 ) { return -1; } else if ( a->m.L17 > b->m.L17 ) { return 1; }
		if ( a->m.L18 < b->m.L18 ) { return -1; } else if ( a->m.L18 > b->m.L18 ) { return 1; }
		if ( a->m.L19 < b->m.L19 ) { return -1; } else if ( a->m.L19 > b->m.L19 ) { return 1; }
		if ( a->m.L20 < b->m.L20 ) { return -1; } else if ( a->m.L20 > b->m.L20 ) { return 1; }
		if ( a->m.L21 < b->m.L21 ) { return -1; } else if ( a->m.L21 > b->m.L21 ) { return 1; }
		if ( a->m.L22 < b->m.L22 ) { return -1; } else if ( a->m.L22 > b->m.L22 ) { return 1; }
		if ( a->m.L23 < b->m.L23 ) { return -1; } else if ( a->m.L23 > b->m.L23 ) { return 1; }
		if ( a->m.L24 < b->m.L24 ) { return -1; } else if ( a->m.L24 > b->m.L24 ) { return 1; }
		if ( a->m.L25 < b->m.L25 ) { return -1; } else if ( a->m.L25 > b->m.L25 ) { return 1; }
		if ( a->m.L26 < b->m.L26 ) { return -1; } else if ( a->m.L26 > b->m.L26 ) { return 1; }
		if ( a->m.L27 < b->m.L27 ) { return -1; } else if ( a->m.L27 > b->m.L27 ) { return 1; }
		if ( a->m.L28 < b->m.L28 ) { return -1; } else if ( a->m.L28 > b->m.L28 ) { return 1; }
		if ( a->m.L29 < b->m.L29 ) { return -1; } else if ( a->m.L29 > b->m.L29 ) { return 1; }
		return 0; //If made it here, a and b are EQUAL!
	}
}

void InitializeMorton(Morton& m)
{
	memcpy(&m,0,sizeof(m));
}


void SetMortonBit(Morton& m,int64 mBit,int bitLevel) 
{

	if( bitLevel > m.LEVEL ) {
		m.LEVEL = bitLevel;
	}
   switch( bitLevel ) {
		case 1:
		   m.L1 = mBit; return;
		case 2:
		   m.L2 = mBit; return;
		case 3:
		   m.L3 = mBit; return;
		case 4:
		   m.L4 = mBit; return;
		case 5:
		   m.L5 = mBit; return;
		case 6:
		   m.L6 = mBit; return;
		case 7:
		   m.L7 = mBit; return;
		case 8:
		   m.L8 = mBit; return;
		case 9:
		   m.L9 = mBit; return;
		case 10:
		   m.L10 = mBit; return;
		case 11:
		   m.L11 = mBit; return;
		case 12:
		   m.L12 = mBit; return;
		case 13:
		   m.L13 = mBit; return;
		case 14:
		   m.L14 = mBit; return;
		case 15:
		   m.L15 = mBit; return;
		case 16:
		   m.L16 = mBit; return;
		case 17:
		   m.L17 = mBit; return;
		case 18:
		   m.L18 = mBit; return;
		case 19:
		   m.L19 = mBit; return;
		case 20:
		   m.L20 = mBit; return;
		case 21:
		   m.L21 = mBit; return;
		case 22:
		   m.L22 = mBit; return;
		case 23:
		   m.L23 = mBit; return;
		case 24:
		   m.L24 = mBit; return;
		case 25:
		   m.L25 = mBit; return;
		case 26:
		   m.L26 = mBit; return;
		case 27:
		   m.L27 = mBit; return;
		case 28:
		   m.L28 = mBit; return;
		case 29:
		   m.L29 = mBit; return;
   }
}

int64 GetMortonBit(Morton m,int bitLevel)
{
   switch( bitLevel ) {
		case 1:
			return m.L1; 
		case 2:
			return m.L2; 
		case 3:
			return m.L3; 
		case 4:
			return m.L4; 
		case 5:
			return m.L5; 
		case 6:
			return m.L6; 
		case 7:
			return m.L7; 
		case 8:
			return m.L8; 
		case 9:
			return m.L9; 
		case 10:
			return m.L10; 
		case 11:
			return m.L11; 
		case 12:
			return m.L12; 
		case 13:
			return m.L13; 
		case 14:
			return m.L14; 
		case 15:
			return m.L15; 
		case 16:
			return m.L16; 
		case 17:
			return m.L17; 
		case 18:
			return m.L18; 
		case 19:
			return m.L19; 
		case 20:
			return m.L20; 
		case 21:
			return m.L21; 
		case 22:
			return m.L22; 
		case 23:
			return m.L23; 
		case 24:
			return m.L24; 
		case 25:
			return m.L25;
		case 26:
			return m.L26; 
		case 27:
			return m.L27; 
		case 28:
			return m.L28; 
		case 29:
			return m.L29; 
   }
}

Morton StringToMorton(std::string mStr)
{
   Morton m;
	int64 bit;
	stringstream ss;
	ss.clear();
	for(unsigned int i = 0; i < mStr.length(); i++ )
	{
		ss << mStr[i];
		ss >> bit;
		SetMortonBit(m,bit-1,i+1);
    	ss.clear();
	}
	return m;
}

int64 StringToHpx(std::string mStr)
{
   int64 hpx = 0;
	int64 bit;
	stringstream ss;
	ss.clear();
	for(unsigned int i = 0; i < mStr.length(); i++ )
	{
		ss << mStr[i];
		ss >> bit;
		hpx = (hpx*10)+bit;
    	ss.clear();
	}
	return hpx;
}

std::string MortonToString(Morton m)
{
	std::string mStr;
	mStr.clear();
	std::string tStr;
	std::ostringstream ss;
	for( int i = 1; i <= GetMortonLevel(m); i++)
	{
		ss << GetMortonBit(m,i)+1;
		tStr = ss.str();
	}
	mStr.append( tStr );
	return mStr;
}


Morton AppendMortonBit(Morton m,int bit)
{
    SetMortonBit(m,bit,m.LEVEL+1);
	return m;
} 

int64 AppendHPXBit(int64 hpxid,int bit)
{
    return (hpxid*10)+bit;
}

std::vector<Morton> SiblingsOfMorton(Morton m)
{
   std::vector<Morton> siblings;
   std::vector<Morton> c;
   Morton s1,s2,s3,s4;
   Morton mParent = ParentOfMorton(m);
   if( mParent.LEVEL == 0 )
   {
	   SetMortonBit(s1,1,1);
	   SetMortonBit(s2,2,1);
	   SetMortonBit(s3,3,1);
	   SetMortonBit(s4,4,1);
	   c.push_back(s1);
	   c.push_back(s2);
	   c.push_back(s3);
	   c.push_back(s4);
   }
   else
   {
      c = ChildrenOfMorton(mParent);
   }

   if( Equals(c[0],m) ) // Child #1
   {
		siblings.push_back(c[1]); // Child #2
		siblings.push_back(c[2]); // Child #3
		siblings.push_back(c[3]); // Child #4
   }
   if( Equals(c[1],m) ) // Child #2
   {
		siblings.push_back(c[0]); // Child #1
		siblings.push_back(c[2]); // Child #3
		siblings.push_back(c[3]); // Child #4
   }
   if( Equals(c[2],m) ) // Child #3
   {
		siblings.push_back(c[0]); // Child #1
		siblings.push_back(c[1]); // Child #2
		siblings.push_back(c[3]); // Child #4
   }
   if( Equals(c[3],m) ) // Child #4
   {
		siblings.push_back(c[0]); // Child #1
		siblings.push_back(c[1]); // Child #2
		siblings.push_back(c[2]); // Child #3
   }
   return siblings;
}

Morton ParentOfMorton(Morton m)
{
   Morton parent = m;
   SetMortonBit(parent,0,m.LEVEL);
   parent.LEVEL = m.LEVEL-1;
   if( parent.LEVEL < 1 ) {
	   parent.LEVEL = 0;
   }
   return parent;
}

std::vector<Morton> ChildrenOfMorton(Morton m)
{
	std::vector<Morton> children;
	children.push_back(AppendMortonBit(m,1));
	children.push_back(AppendMortonBit(m,2));
	children.push_back(AppendMortonBit(m,3));
	children.push_back(AppendMortonBit(m,4));
	return children;
}


int64 GetBitMask(int order)
{
	int64 mask=3;
	int64 shift = (2*order)-2;
	switch( shift )
	{
	case 1: mask = 3 << shift; break;
	case 2: mask = 3 << shift; break;
	case 3: mask = 3 << shift; break;
	case 4: mask = 3 << shift; break;
	case 5: mask = 3 << shift; break;
	case 6: mask = 3 << shift; break;
	case 7: mask = 3 << shift; break;
	case 8: mask = 3 << shift; break;
	case 9: mask = 3 << shift; break;
	case 10: mask = 3 << shift; break;
	case 11: mask = 3 << shift; break;
	case 12: mask = 3 << shift; break;
	case 13: mask = 3 << shift; break;
	case 14: mask = 3 << shift; break;
	case 15: mask = 3 << shift; break;
	case 16: mask = 3 << shift; break;
	case 17: mask = 3 << shift; break;
	case 18: mask = 3 << shift; break;
	case 19: mask = 3 << shift; break;
	case 20: mask = 3 << shift; break;
	case 21: mask = 3 << shift; break;
	case 22: mask = 3 << shift; break;
	case 23: mask = 3 << shift; break;
	case 24: mask = 3 << shift; break;
	case 25: mask = 3 << shift; break;
	case 26: mask = 3 << shift; break;
	case 27: mask = 3 << shift; break;
	case 28: mask = 3 << shift; break;
	case 29: mask = 3 << shift; break;
	case 30: mask = 3 << shift; break;
	case 31: mask = 3 << 1; mask = mask << 30; break;
	case 32: mask = 3 << 1; mask = mask << 31; break;
	case 33: mask = 3 << 1; mask = mask << 32; break;
	case 34: mask = 3 << 2; mask = mask << 32; break;
	case 35: mask = 3 << 3; mask = mask << 32; break;
	case 36: mask = 3 << 4; mask = mask << 32; break;
	case 37: mask = 3 << 5; mask = mask << 32; break;
	case 38: mask = 3 << 6; mask = mask << 32; break;
	case 39: mask = 3 << 7; mask = mask << 32; break;
	case 40: mask = 3 << 8; mask = mask << 32; break;
	case 41: mask = 3 << 9; mask = mask << 32; break;
	case 42: mask = 3 << 10; mask = mask << 32; break;
	case 43: mask = 3 << 11; mask = mask << 32; break;
	case 44: mask = 3 << 12; mask = mask << 32; break;
	case 45: mask = 3 << 13; mask = mask << 32; break;
	case 46: mask = 3 << 14; mask = mask << 32; break;
	case 47: mask = 3 << 15; mask = mask << 32; break;
	case 48: mask = 3 << 16; mask = mask << 32; break;
	case 49: mask = 3 << 17; mask = mask << 32; break;
	case 50: mask = 3 << 18; mask = mask << 32; break;
	case 51: mask = 3 << 19; mask = mask << 32; break;
	case 52: mask = 3 << 20; mask = mask << 32; break;
	case 53: mask = 3 << 21; mask = mask << 32; break;
	case 54: mask = 3 << 22; mask = mask << 32; break;
	case 55: mask = 3 << 23; mask = mask << 32; break;
	case 56: mask = 3 << 24; mask = mask << 32; break;
	case 57: mask = 3 << 25; mask = mask << 32; break;
	case 58: mask = 3 << 26; mask = mask << 32; break;
	case 59: mask = 3 << 27; mask = mask << 32; break;
	case 60: mask = 3 << 28; mask = mask << 32; break;
	case 61: mask = 3 << 29; mask = mask << 32; break;
	case 62: mask = mask << 30; mask = mask << 32; break;
	case 63: mask = mask << 31; mask = mask << 32; break;
	case 64: mask = mask << 32; mask = mask << 32; break;
	}
	return mask;
}

Morton HpxToMorton(int64 hpxid,int order)
{
	Morton m;
	int64 mask = 0;
	int64 nextBit = 0;
	// Process the 2-bit groups
	mask = GetBitMask(order);
	for(int i = order; i > 0; i--)  
	{
	  nextBit = ((hpxid & mask) >> ((2*i)-2) ); //Morton bits still in range 0,1,2,3 
	  SetMortonBit(m,nextBit,(order-i)+1);
	  mask = mask >> 2; //operator.rshift(mask,2) 
	}
	return m;	
}

void MortonToHpx(Morton m,int64& hpxid, int& order)
{
	int64 nextBit = 0;
	order = GetMortonLevel(m);
	hpxid = 0;
	int temp = 0;
	int64 shift = 0;
	for(int i = order; i >= 1; i--)
	{
	  nextBit = GetMortonBit(m,i);
	  nextBit = nextBit << shift; //operator.lshift(nextBit,shift)
	  hpxid += nextBit;  
	  shift += 2;
	}
}
 
Morton PhiThetaToMorton(double longHpx,double coLatHpx,int order)
{
   pointing pt;
   pair<int64,int> hpxId;
   pt.theta = coLatHpx;
   pt.phi = longHpx;
   hpxQ.Set(order,NEST);
   hpxId.first = hpxQ.ang2pix(pt);
   hpxId.second = order;

   // Normalize the HEALPix index to Base 0 indexing through shift
   hpxId.first = hpxId.first - hpxQ.FaceNum(hpxId.first,hpxId.second)*order_to_npface[order];

   Morton m = HpxToMorton(hpxId.first,hpxId.second);
   return m;
}

pointing MortonToPhiTheta(Morton m)
{
   pointing pt;
   int64 hpxId;
   int order;
   MortonToHpx(m,hpxId,order);
   hpxQ.Set(order,NEST);
   
   // De-Normalize the HEALPix index from Base 0 indexing to global indexing 
   hpxId = hpxId + hpxQ.FaceNum(hpxId,order)*order_to_npface[order];

   pt = hpxQ.pix2ang(hpxId);

   return pt;
}


bool SortFunctionMorton2( MortonNode a, MortonNode b )
{
	// If a.m < b.m return true
	if( LessThan( a.m, b.m ) ) {
		return true;
	}

	// If a.m == b.m need to check the sub index 
	if( Equals( a.m, b.m ) ) {
		if( a.sub < b.sub ) {
			return true;
		}
	}
	return false;
}


//#######
//####### MORTON NODE METHODS 
//#######

MortonNode::MortonNode() :
	sub(1),
	childrenYN(0),
	phi(0.0),
	theta(0.0),
	sentinel(-1)
{
}

MortonNode::MortonNode( Morton _m, int _sub, int _childrenYN, double _phi, double _theta, int _data, int _facenum) :
	m(_m),
    sub(_sub),
	childrenYN(_childrenYN),
	phi(_phi),
	theta(_theta),
	facenum(_facenum)
{
	data.push_back(_data);
	sentinel = -1;
}

template<class MortonNode>
class HasPackedId {
    Morton _m;
	int _sub;
public:
    HasPackedId(Morton m,int sub) : _m(m),_sub(sub) {}
    bool operator()(const MortonNode & o) const {
        return ( o.m == _m && o.sub == _sub );
    }
};

void MortonLQT::ClearTree()
{
	mTree.clear();
}

int MortonLQT::GetNumMortonNodes()
{
  return( mTree.size() );
}

std::vector<int> MortonLQT::GetMortonNodeSizeHistogram()
{
	std::vector<int> node_size_histo;
	node_size_histo.clear();
    node_size_histo.push_back(0);
    node_size_histo.push_back(0);

	unsigned int i;
	bool done = false;

	Morton current,next;
	int proximate_count = 1;

	// Find the first leaf node, else just
	// find last node in tree.
	unsigned int n = 0;
	for( i = 0; i < mTree.size(); i++  )
	{
		// Look for FIRST leaf node
		if( mTree[i].childrenYN == false )
		{
			current = mTree[i].m;
			n = i+1; //Location of node next to FIRST leaf node
			done = true;
		}
		
		// Check for out of bounds conditions, if so, just return empty histogram
		if( n >= mTree.size() )
		{
			return node_size_histo;
		}

		if( done == true) 
		{
			i = mTree.size(); //break out of loop
		}
	}

	for( i = n; i < mTree.size(); i++)
	{
		// Check if current Morton Node is a leaf node
		if( mTree[i].childrenYN == false )
		{
			// If leaf node, count number of matching, proximate, Morton 
			// addresses.  This count is then used to update the node size histogram.  
			// If current proximate count is larger than the histogram, expand histogram 
			// until histogram is at least the same size as proximate count.
			// Otherwise increment histogram element n where n is the current proximate count 
			// i.e. element 0 represents Morton addresses with zero proximate count 
			// (Morton address is unique), element 8 represents Morton addresses with a 
			// proximate count of 8, etc...
			next = mTree[i].m;
			if( Equals(next,current) == true )
			{
				proximate_count++;
			}
			// New Morton address so update histogram and set new "current"
			else 
			{				
				// Now update histogram, expanding it if need be...
				if( proximate_count >= node_size_histo.size() )
				{
					while( proximate_count >= node_size_histo.size() )
					{
						node_size_histo.push_back(0);
					}
				}
				node_size_histo[proximate_count]++;
				proximate_count = 1;
				current = next;
			}	
		}
	}
	return node_size_histo;
}

double MortonLQT::GetAvgMortonTreeDepth()
{
	int64 sum = 0;
	if( GetNumMortonNodes() == 0 ) {
		return 0;
	}
	for( int i = 0; i < GetNumMortonNodes(); i++ ) {
		sum = sum + GetMortonLevel(mTree[i].m);
	}
	return (double)sum / (double)GetNumMortonNodes();
}

int MortonLQT::GetMortonTreeDepth()
{
  int maxdepth = 0;
  int level = 0;
  // Check for zero size tree
  if( GetNumMortonNodes() == 0)
	  return 0;
  for( int i = 0; i < GetNumMortonNodes(); i++ ) {     
     level = GetMortonLevel(mTree[i].m);
     if( level > maxdepth)
	 {
        maxdepth = level;
	 }
  }
  return( maxdepth ); 
}

int MortonLQT::GetMortonTreeMinDepth()
{
	int mindepth = 0;
	// Check for zero size tree
	if( GetNumMortonNodes() == 0)
		return 0;
	for( int i = 0; i < GetNumMortonNodes(); i++ ) {     
		if( mTree[i].childrenYN == false) {
			return GetMortonLevel(mTree[i].m);
		}
	}
}


void MortonLQT::SetFaceNum(int face_num)
{
	faceNum = face_num;
}

int MortonLQT::FindIndexAtMorton(Morton m)
{
  return( FindIndexAtMortonSub_BinarySearch_internal(m,1) );
}

int MortonLQT::FindIndexAtMortonSub(Morton m,int sub)
{
	// Do Linear search if list size is relatively "small", quicker
	// than binary search for "small" lists. Otherwise, do binary search
	// for larger lists.
	if( mTree.size() < 100 ) {
		return( FindIndexAtMortonSub_internal(m,sub) );
	}
	return( FindIndexAtMortonSub_BinarySearch_internal(m,sub) );
}

// Simple Linear Search of sTree vector until matching Morton,Sub combination is found. Returns index of sTree where
// match occured. If no match is found return -1.
int MortonLQT::FindIndexAtMortonSub_internal(Morton sMorton, int sSub)
{
   for(unsigned int i = 0; i < mTree.size(); i++ )
   {
     if( (Equals(mTree[i].m,sMorton) && mTree[i].sub == sSub) )
	  {
         return i;
	  }

   }
   return -1;
}

int MortonLQT::FindIndexAtMortonSub_BinarySearch_internal(Morton sMorton, int sSub)
{
	int min = 0, max = mTree.size()-1, mid;

	// Check for sSub = -1, set to 1 to find first in block of appended MortonNodes
	if( sSub == -1 ) sSub = 1;

	while ( min <= max ) {
		mid = (int) ((min+max)/2);
		if( LessThan(sMorton,mTree[mid].m) ) {
			max = mid - 1;
		}
		else if ( GreaterThan(sMorton,mTree[mid].m) ) {
			min = mid + 1;
		}
		else {
			// Now must compare the Morton Subaddress
			if( sSub < mTree[mid].sub ) {
				max = mid - 1;
			} else if ( sSub > mTree[mid].sub ) {
				min = mid + 1;
			}
			else {
				return mid;
			}
		}
	}
	return -1;
}

//Inputs: pointing angle = HPX phi(phi)(radians), HPX theta (theta) (radians)
std::vector<MortonNode> MortonLQT::SearchMortonNodeAtPhiTheta(pointing pt,int sentinel)
{

  Morton m = FindMortonAtPhiTheta(pt);
  return SearchMortonNode(m,-1,sentinel);
}

std::vector<MortonNode> MortonLQT::SearchMortonNodeHpx(int64 hpxid,int order,int sentinel)
{
	Morton m;
	std::vector<MortonNode> temp,found;
	//Check for possible base node search: hpxid = 0 and order = 0
	if( hpxid == 0 && order == 0 ) 	   
	{
	   temp.clear();
	   m = HpxToMorton(0,1);
	   temp = SearchMortonNode(m,-1,sentinel);
	   found.insert(found.end(),temp.begin(),temp.end());
	   
	   temp.clear();
   	   m = HpxToMorton(1,1);
	   temp = SearchMortonNode(m,-1,sentinel);
	   found.insert(found.end(),temp.begin(),temp.end());

	   temp.clear();
   	   m = HpxToMorton(2,1);
	   temp = SearchMortonNode(m,-1,sentinel);
	   found.insert(found.end(),temp.begin(),temp.end());

	   temp.clear();
   	   m = HpxToMorton(3,1);
	   temp = SearchMortonNode(m,-1,sentinel);
	   found.insert(found.end(),temp.begin(),temp.end());
	   temp.clear();
	   return found;
	}
	else
	{
	   m = HpxToMorton(hpxid,order);
	   return SearchMortonNode(m,-1,sentinel);
	}
}

std::vector<MortonNode> MortonLQT::SearchMortonNode(Morton m,int sub,int sentinel)
{
	std::vector<MortonNode> found_nodes;
	// If specified root node then search entire tree
	if( m.LEVEL == 0 )
	{
		m.L1 = 1;
		SearchMortonNode_internal(found_nodes,m,sub,sentinel);
		m.L1 = 2;
		SearchMortonNode_internal(found_nodes,m,sub,sentinel);
		m.L1 = 3;
		SearchMortonNode_internal(found_nodes,m,sub,sentinel);
		m.L1 = 4;
		SearchMortonNode_internal(found_nodes,m,sub,sentinel);
	}
	else
	{
	   SearchMortonNode_internal(found_nodes,m,sub,sentinel);
	}

	return( found_nodes );
}

void MortonLQT::SearchMortonNode_internal(std::vector<MortonNode>& found_nodes,Morton m,int sub,int sentinel)
{

	// If Node has children, descend to them while keeping list of found data nodes
	// along the way. Return when no children are found.
	std::vector<MortonNode> nodes = GetNodeAtMorton(m,sub);
	  
	for(unsigned int i = 0; i < nodes.size(); i++)
	{
	   if( nodes[i].data.size() > 0 ) 
	   {
		   // Check for current sentinel value, if found don't add to found
		   // list as this node has already been discovered.
		   if( nodes[i].sentinel != sentinel ) {
				found_nodes.push_back(nodes[i]);

				// Set sentinel value on data node marking as FOUND
				MortonNode update_node;
				update_node.m = nodes[i].m;
				update_node.sentinel = sentinel;
				UpdateMortonNodeSentinel(update_node);
		   }
	   }
	}
	if( nodes.size() > 0 && nodes[0].childrenYN == 1 )
	{
		// Compute children of node
		std::vector<Morton> children = ChildrenOfMorton(m);	   
		SearchMortonNode_internal(found_nodes,children[0],sub,sentinel);
		SearchMortonNode_internal(found_nodes,children[1],sub,sentinel);
		SearchMortonNode_internal(found_nodes,children[2],sub,sentinel);
		SearchMortonNode_internal(found_nodes,children[3],sub,sentinel);
	}
}

Morton MortonLQT::ClosestMortonToPhiTheta(std::vector<MortonNode> NodeList,pointing pt)
{
   double deltaNet;
   double minDist = 99999999.0;
   Morton closest;
   int sub = -1;
   for( unsigned int i = 0; i < NodeList.size(); i++ )
   {
   	  deltaNet = acos(fabs(cosdist_zphi(cos(pt.theta),pt.phi,cos(NodeList[i].theta),NodeList[i].phi)));
      if( deltaNet < minDist )
	  {
         minDist = deltaNet;
         closest = NodeList[i].m;
   	  }
   }
   return closest;
}

Morton MortonLQT::FindMortonAtPhiTheta(pointing pt)
{
	int sub = 1;	
	bool done = false;
	int level = 1;
	Morton m;
	std::vector<MortonNode> node;

	while( done == false )
	{
		m = PhiThetaToMorton(pt.phi,pt.theta,level);
		node = GetNodeAtMorton(m,-1);   
		// Return empty list if nothing found
		if( node.size() == 0 )
			return m;
		if( node[0].childrenYN == 1 )
		{
		  level += 1;
		}
		else
		{
		  // Possibility of multiple data nodes with duplicate
		  // Morton Code, must examine each duplicate to find
		  // data node with closest matching phi, theta
		  if( node.size() > 1 )
		  {
			 m = ClosestMortonToPhiTheta(node,pt);
		  }
		  done = true;
		}
	}
	return( m );
}


MortonNode MortonLQT::GetNodeAtIndex(unsigned int index)
{
	return mTree[index];
}


std::vector<MortonNode> MortonLQT::GetNodeAtMorton(Morton m,int sub)
{
	std::vector<MortonNode> found_node_list;
	MortonNode found_node;
	int index;
	found_node_list.clear();

	// Do we want ALL indices under morton returned OR
	// just the actual morton,sub pair?
	// Return all matching morton if "sub" equals -1
	if( sub == -1 )
	{
		index = FindIndexAtMorton(m);
	}
	else
	{
		index = FindIndexAtMortonSub(m,sub);
	}

	// Node NOT Found, return empty list
	if( index == -1 )
	{
		return found_node_list;
	}
	
	// Otherwise, found the Node so populate the returnable
	found_node.m = mTree[index].m;
	found_node.sub = mTree[index].sub;
	found_node.phi = mTree[index].phi;
	found_node.theta = mTree[index].theta;
	found_node.childrenYN = mTree[index].childrenYN;
	found_node.data = mTree[index].data;
	found_node.sentinel = mTree[index].sentinel;
	found_node.facenum = mTree[index].facenum;
	found_node_list.push_back(found_node);

    // Because of possibility of concatenated, duplicate Morton Nodes we must
	// check to see if there are any other Morton Nodes mapped to same Morton Code
	if( sub == -1)
	{
		bool done = false;
		int numNodes = GetNumMortonNodes();
		while( done == false )
		{
		   index += 1;
		   if( index < numNodes)
		   {
			  if( Equals(mTree[index].m,m) )
			  {
				 found_node.m = mTree[index].m;
				 found_node.sub = mTree[index].sub;
				 found_node.phi = mTree[index].phi;
				 found_node.theta = mTree[index].theta;
				 found_node.childrenYN = mTree[index].childrenYN;
				 found_node.data = mTree[index].data;
				 found_node.sentinel = mTree[index].sentinel;
				 found_node.facenum = mTree[index].facenum;
				 found_node_list.push_back(found_node);
			  }
			  else
			  {
				 done = true;
			  }
		   }
		   else
		   {
			  done = true;
		   }
		}
	}

	return( found_node_list );
}

void MortonLQT::UpdateMortonNodeSentinel(MortonNode node)
{
    // Find index of list MortonNode objects that contain
	// Morton code = m, could include one or more subaddresses with same m!
	int index = FindIndexAtMorton(node.m);

	if( index != -1 )
	{
		// Update sentinel of first MortonNode
		mTree[index].sentinel = node.sentinel;

		// Now check for possible duplicate Morton addresses and update 
		// their sentinel values as well.
		bool done = false;
		int numNodes = GetNumMortonNodes();
		while( done == false )
		{
		   index += 1;
		   if( index < numNodes)
		   {
			  if( Equals(mTree[index].m,node.m) )
			  {
				 mTree[index].sentinel = node.sentinel;
			  }
			  else
			  {
				 done = true;
			  }
		   }
		   else
		   {
			  done = true;
		   }
		}
	}
}

void MortonLQT::UpdateMortonNode(MortonNode node)
{
    // Find index of list MortonNode object that contains 
	// Morton code = m
	int index = FindIndexAtMorton(node.m);
	if( index != -1 )
	{
		mTree[index].m = node.m;
		mTree[index].sub = node.sub;
		mTree[index].childrenYN = node.childrenYN;
		mTree[index].phi = node.phi;
		mTree[index].theta = node.theta;
		mTree[index].data = node.data;
		mTree[index].sentinel = node.sentinel;
		mTree[index].facenum = node.facenum;
	}
}



void MortonLQT::AppendMortonNode(MortonNode node)
{
    // Must compute the sub of this duplicate m code
    // Subaddress is based on first come first served
	std::vector<MortonNode> found_nodes;
	found_nodes = GetNodeAtMorton(node.m,-1);

	// Want to know the highest sub-address so take the last MortonNode in the found_nodes
	// list and return it and use its "sub" address.
	int sub = found_nodes.back().sub;
      
    MortonNode new_node;
	new_node.m = node.m;
	new_node.sub = sub+1;
	new_node.childrenYN = node.childrenYN;
	new_node.phi = node.phi;
	new_node.theta = node.theta;
	new_node.data = node.data;
	new_node.sentinel = node.sentinel;
	new_node.facenum = node.facenum;
	mTree.push_back(new_node);
	// Sort the tree
	std::sort( mTree.begin(), mTree.end(), SortFunctionMorton2);
}

void MortonLQT::AddMortonNode2(MortonNode node)
{
	MortonNode new_node;
	new_node.m = node.m;
	new_node.sub = node.sub;
	new_node.childrenYN = node.childrenYN;
	new_node.phi = node.phi;
	new_node.theta = node.theta;
	new_node.data = node.data;
	new_node.sentinel = node.sentinel;
	new_node.facenum = node.facenum;
	mTree.push_back(new_node);
	// Sort the tree
	std::sort( mTree.begin(), mTree.end(), SortFunctionMorton2);
}

void MortonLQT::AddMortonNode(MortonNode node)
{
	MortonNode new_node;
	new_node.m = node.m;
	new_node.sub = node.sub;
	new_node.childrenYN = node.childrenYN;
	new_node.phi = node.phi;
	new_node.theta = node.theta;
	new_node.data = node.data;
	new_node.sentinel = node.sentinel;
	new_node.facenum = node.facenum;
	mTree.push_back(new_node);
	// Now must add the other sibling nodes at same level of quad tree
	// Sort the tree
	std::sort( mTree.begin(), mTree.end(), SortFunctionMorton2);

	// Lastly, make sure parent (if exists) of new nodes is marked as having children
	if( ParentOfMorton(node.m).LEVEL != 0 )
	{
		MortonNode update_node;
		update_node.m = ParentOfMorton(node.m);
		update_node.phi = 0.0;
		update_node.theta = 0.0;
		update_node.childrenYN = 1;
		update_node.data.clear();
		update_node.facenum = faceNum;
		UpdateMortonNode(update_node);
	}
}

Morton MortonLQT::InsertMortonNodeAtMorton(Morton m,std::vector<int64> data)
{
	// Convert Morton code to phi & theta
	// then insert
	pointing pt;
	pt = MortonToPhiTheta(m);
	return( InsertMortonNodeAtPhiTheta(pt,data) );
}

Morton MortonLQT::InsertMortonNodeAtHpx(int64 hpxid,int order,std::vector<int64> data)
{
	// Convert HEAPix index to Morton code then insert via InsertNodeMorton
	Morton m = HpxToMorton(hpxid,order);
	return( InsertMortonNodeAtMorton(m,data) );
}

Morton MortonLQT::InsertMortonNodeAtPhiTheta(pointing pt,std::vector<int64> data)
{
	return InsertMortonNode_internal_1(pt,data);
}

Morton MortonLQT::InsertMortonNode_internal_1(pointing pt,std::vector<int64> data)
{
	int treeDepth = 1;
	Morton new_m;
	bool done = false;
	pointing pt_save;
	int childrenYN_save;
	std::vector<int64> data_save;


	// Create new Morton Node to be inserted based on its phi & theta
	MortonNode insert_node;
	insert_node.m = PhiThetaToMorton(pt.phi,pt.theta,treeDepth);
	insert_node.phi = pt.phi;
	insert_node.theta = pt.theta;
	insert_node.data = data;
	insert_node.childrenYN = 0;
	insert_node.facenum = faceNum;

	// Insert the new Morton node checking for data collision (Morton
	// node who's "data" element is NOT equal to -1). When collision
	// occurs create new MortonNode's who's Morton codes are the children
	// of the Morton code to be inserted. We recompute the Morton code to 
	// next level deeper in tree so will match one of the computed child
	// Morton codes. If there is no data collision for matched Morton node
	// then we simply update the matched Morton node with data to be inserted.
	// Last possibility is that the Morton node is not found in the list. If
	// so we create new Morton node and add it to the list. In all cases of 
	// inserted Morton nodes we finish by sorted the list by ascending Morton
	// code
	while( done == false ){
       // Search for Morton node in tree, node returned if found, empty vector otherwise
	   std::vector<MortonNode> node = GetNodeAtMorton(insert_node.m,1);
 
       // If Morton node is found we then need to check for data collision, else update
       // Morton node with data
       if( node.size() > 0 ) {
         
         // If Morton node has children then need to re-calculate Morton
         // code of node to be inserted to next level deeper.
         if( node[0].childrenYN == 1 ){
            
            treeDepth += 1;
            new_m = PhiThetaToMorton(pt.phi,pt.theta,treeDepth);
            insert_node.m = new_m;  
		 }
         // Otherwise check to see if there is a data node insertion collision or not
         else {
            
			// If "data" is populated in MortonNode (not empty vector) then
            // we have a node insertion collision!
            if( node[0].data.size() > 0 ){
              
              // Store collided Morton Node's information for re-insertion
              pt_save.phi = node[0].phi;
              pt_save.theta = node[0].theta;
              childrenYN_save = node[0].childrenYN;
              data_save = node[0].data;
                  
              // Re-calculate next level deeper Morton code for Morton Node to be inserted. 
              treeDepth += 1;

			  if( abs_approx(node[0].phi,insert_node.phi) && abs_approx(node[0].theta,insert_node.theta) ) {
					node[0].data.push_back(insert_node.data[0]);
					UpdateMortonNode(node[0]);
					done = true;
			  } else {
				  // Check if have reached user specified maximum tree depth.
				  // If so will overwrite the OLD Node with NEW Node's phi,theta and data index.
				  if( treeDepth > userMaxDepth ){
	                 
					 AppendMortonNode(insert_node);
					 done = true;
				  }
				  // Re-insert the new MortonNode and collided MortonNode at next level in Morton Tree (higher resolution)
				  else {
					   new_m = PhiThetaToMorton(insert_node.phi,insert_node.theta,treeDepth);
					   insert_node.m = new_m;
	                     
					   // Add children Morton nodes to Morton list, populate the matching
					   // child Morton == re-calculated Morton code node with "data" to be
					   // inserted, other child Morton nodes have "data" == -1
					   AddMortonNode(insert_node);	

					   // Now must prepare the collided node(s) to move to new node at next 
					   // lower level of the tree.
					   new_m = PhiThetaToMorton(pt_save.phi,pt_save.theta,treeDepth);
					   insert_node.m = new_m;
					   insert_node.phi = pt_save.phi;
					   insert_node.theta = pt_save.theta;
					   insert_node.childrenYN = childrenYN_save;
					   insert_node.data = data_save;
				  }
			  }
			}
            // Otherwise no insertion collision so update the data,phi,theta,etc
            // attributes of MortonNode.
            else{
              UpdateMortonNode(insert_node);
              if( treeDepth > curTreeDepth ){
                curTreeDepth = treeDepth;
			  }
              done = true;
			}
		 }
	   }            
       else{
         // If Morton code not found in list then append it to the end
         // of the Morton list then sort the list ascending by Morton code
         AddMortonNode(insert_node);
         if( treeDepth > curTreeDepth ){
           curTreeDepth = treeDepth;
		 }
         done = true; 
	   }
	 }      
     return( insert_node.m );
}


void MortonLQT::DeleteMortonNodeAtPhiTheta(pointing pt)
{
	Morton m = FindMortonAtPhiTheta(pt);
	DeleteMortonNodeAtMorton(m,1);
}

void MortonLQT::DeleteMortonNodeAtHpx(int64 hpxid,int order)
{
	Morton m = HpxToMorton(hpxid,order);			      
	DeleteMortonNodeAtMorton(m,1);      
}

void MortonLQT::DeleteMortonNodeAtMorton(Morton m,int sub)
{
	pointing pt;
	// Find the tree index of the node to be deleted
	int index = FindIndexAtMortonSub(m,sub);

	if( index == -1 )
	{
		return;
	}
	// Remove the node from the tree
	mTree.erase(mTree.begin()+index);

	// Run through entire tree and find leaf nodes
	// If leaf node is found, insert it into new tree
	std::vector<MortonNode> oldTree;
	oldTree = mTree;
	mTree.clear();
	for(unsigned int i = 0; i < oldTree.size(); i++)
	{
		if( oldTree[i].data.size() > 0 )
		{
			pt.phi = oldTree[i].phi;
			pt.theta = oldTree[i].theta;
			InsertMortonNodeAtPhiTheta(pt,oldTree[i].data);
		}
	}
}

void MortonLQT::SaveTreeToFile(ofstream& fp)
{
	// Write out header
	WriteHeader(fp);

	// Write out data
	WriteData(fp);	 

}

void MortonLQT::LoadTreeFromFile(ifstream& fp)
{
	int numNodes, nodeNumber, sub, num_data_indices, facenum;
	int64 next_data;
	std::string m;
	int childrenYN;
	double phi,theta;
	std::string temp;
          
    // Reset current tree
    mTree.clear();    

	// First line is Base Number, skip
    skip_line(fp);

    // Next line is number of records
	fp >> temp >> numNodes;  skip_line(fp);

    // Skip header line
	skip_line(fp);

    // Next lines are data records (the nodes of the tree)
    for(int i = 0; i < numNodes; i++)
	{
      MortonNode next_node;

	  fp >> nodeNumber >> facenum >> m >> sub >> childrenYN >> phi >> theta;
		 
	  next_node.m = StringToMorton(m);
      next_node.sub = sub;
	  next_node.childrenYN = childrenYN;
      next_node.phi = phi;
      next_node.theta = theta;
	  next_node.facenum = facenum;

	  // Now read in count of data indices mapped to this record
	  fp >> num_data_indices;
		 
	  for(int j = 0; j < num_data_indices; j++) {
		fp >> next_data;
		next_node.data.push_back(next_data);
      }

	  // Skip the rest of the line
	  skip_line(fp);


      mTree.push_back(next_node);
	}
}

void MortonLQT::WriteHeader(ofstream& fp)
{     
     // Write out count of total number of nodes (records in file)
	fp << "BASE\t" << faceNum << '\n';
	fp << "NODES\t" << GetNumMortonNodes() << '\n';
	fp << "INDEX\t";
	fp << "MORTON\t";
	fp << "SUB\t";
	fp << "CHILDRENYN\t";
	fp << "PHI\t";
	fp << "THETA\t";
	fp << "#DATA\t";
	fp << "DATA\t\n";
}

void MortonLQT::WriteData(ofstream& fp)
{
	unsigned int j;
	for( int i = 0; i < GetNumMortonNodes(); i++)
	{
		fp << i << "\t"; // Index
		fp << MortonToString(mTree[i].m).c_str() << "\t"; // Morton
		fp << mTree[i].sub << "\t"; // Sub Index
		fp << mTree[i].childrenYN << "\t"; // Children YN
		fp << setprecision(16) << mTree[i].phi << "\t"; // Phi
		fp << setprecision(16) << mTree[i].theta << "\t"; // Theta
		fp << mTree[i].data.size() << "\t"; // # of Data Indices
		for( j = 0; j < mTree[i].data.size(); j++) {
			fp << mTree[i].data[j] << "\t"; // Data Indices
		}
		fp << "\n";
	}
}

void MortonLQT::PrintMortonTree()
{
	 pointing p1;
	 pair<int64,int> hpxId;
	 unsigned int j;
     cout << endl;
	 cout << "#################################################################################" << endl;
     cout << "#######                MORTON    LINEAR    QUAD    TREE                   #######" << endl;
	 cout << "#################################################################################" << endl;
     cout << "Index, Facenum, Morton, SubAddr, ChildrenYN, Phi, Theta, Data, HPXid, HPXorder, Sentinel" << endl;
     cout << "=================================================================================" << endl;
     for( int i = 0; i < GetNumMortonNodes(); i++)
	 {
		// Get corresponding HPX index
		p1.phi = mTree[i].phi;
		p1.theta = mTree[i].theta;
        hpxQ.Set(GetMortonLevel(mTree[i].m),NEST);
        hpxId.first = hpxQ.ang2pix(p1);
        hpxId.second = GetMortonLevel(mTree[i].m);
		   cout << i << ", ";
		   cout << mTree[i].facenum << ", ";
		   PrintMorton(mTree[i].m); cout << ", ";
		   cout	<< mTree[i].sub << ", "
			    << mTree[i].childrenYN << ", "
			    << setprecision(16) << mTree[i].phi << ", "
			    << setprecision(16) << mTree[i].theta << ", ";
        		for(j = 0; j < mTree[i].data.size(); j++)cout << mTree[i].data[j] << ", "; // Data Indices
		   cout << hpxId.first << ", "
				<< hpxId.second << ", ";
		   cout << mTree[i].sentinel << endl;

	 }
     cout << "======================================================" << endl;
     cout << "Number of nodes: " << GetNumMortonNodes() << endl;
     cout << "Maximum depth of tree: " << GetMortonTreeDepth() << endl;
}


void MortonLQT::PrintMortonList()
{
	for( int i = 0; i < GetNumMortonNodes(); i++)
	{
		PrintMorton(mTree[i].m); cout << endl;
	}
	cout << endl;
}