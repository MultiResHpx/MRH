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

#include "MortonLQT.h"


//#######
//####### MORTON NODE METHODS 
//#######

MortonNode::MortonNode() :
	morton(0),
	sub(1),
	childrenYN(0),
	longitude(0.0),
	latitude(0.0),
	data(-1),
	facenum(0)
{
}

MortonNode::MortonNode( int64 _morton, int _sub, int _childrenYN, float _longitude, float _latitude, int _data, int _facenum) :
	morton(_morton),
	sub(_sub),
	childrenYN(_childrenYN),
	longitude(_longitude),
	latitude(_latitude),
	data(_data),
	facenum(_facenum)
{
}

template<class MortonNode>
class HasId {
    int64 _morton;
	int _sub;
public:
    HasId(int64 morton,int sub) : _morton(morton),_sub(sub) {}
    bool operator()(const MortonNode & o) const {
        return ( o.morton == _morton && o.sub == _sub );
    }
};



//#######
//####### MORTON LINEAR QUAD TREE METHODS 
//#######

MortonLQT::MortonLQT() :
   userMaxDepth(1)
{
	hpxQ.Set(1,NEST);
}

MortonLQT::MortonLQT( int _depth ) :
   userMaxDepth(_depth)
{
	hpxQ.Set(1,NEST);
}


MortonLQT::~MortonLQT()
{
  //Clean up memory
  mTree.clear();
}

int MortonLQT::GetMortonLevel(int64 morton)
{
   int length = 0;
   while( morton >= 1)
	{
      morton /= 10;
      length += 1;
	}
   return length;
}

int64 MortonLQT::GetMortonBit(int64 morton,int bitLevel)
{
    int mBit = 0;
	int curlvel = 0;
	int64 digit = 0;
	int mortonLevel = 0;
    mortonLevel = GetMortonLevel(morton);
    if(mortonLevel == 1)
       return morton;
       
    int curLevel = mortonLevel;
    while( morton )
	 {
       digit = morton % 10;
       morton /= 10;
       if( curLevel == bitLevel )
		 {
          return digit;
		 }
       curLevel -= 1;
	 }
    return mBit;
}

int64 MortonLQT::AppendMortonBit(int64 morton,int bit)
{
    return (morton*10)+bit;
}   

int64 MortonLQT::AppendHPXBit(int64 hpxid,int bit)
{
    return (hpxid*10)+bit;
}

std::vector<int64> MortonLQT::SiblingsOfMorton(int64 morton)
{
   std::vector<int64> siblings;
   std::vector<int64> c;
   int64 mParent = ParentOfMorton(morton);
   if( mParent == -1 )
   {
	   c.push_back(1);
	   c.push_back(2);
	   c.push_back(3);
	   c.push_back(4);
   }
   else
   {
      c = ChildrenOfMorton(mParent);
   }
   if( c[0] == morton ) // Child #1
   {
		siblings.push_back(c[1]); // Child #2
		siblings.push_back(c[2]); // Child #3
		siblings.push_back(c[3]); // Child #4
   }
   if( c[1] == morton ) // Child #2
   {
		siblings.push_back(c[0]); // Child #1
		siblings.push_back(c[2]); // Child #3
		siblings.push_back(c[3]); // Child #4
   }
   if( c[2] == morton ) // Child #3
   {
		siblings.push_back(c[0]); // Child #1
		siblings.push_back(c[1]); // Child #2
		siblings.push_back(c[3]); // Child #4
   }
   if( c[3] == morton ) // Child #4
   {
		siblings.push_back(c[0]); // Child #1
		siblings.push_back(c[1]); // Child #2
		siblings.push_back(c[2]); // Child #3
   }
   return siblings;
}


int64 MortonLQT::ParentOfMorton(int64 morton)
{
   int64 mDiv10 = int64(morton/10);
   
   if( mDiv10 == 0 )
   {
      return -1;
   }
   return mDiv10;
}

std::vector<int64> MortonLQT::ChildrenOfMorton(int64 morton)
{
	std::vector<int64> children;
	children.push_back(AppendMortonBit(morton,1));
	children.push_back(AppendMortonBit(morton,2));
	children.push_back(AppendMortonBit(morton,3));
	children.push_back(AppendMortonBit(morton,4));
	return children;
}

int64 MortonLQT::GetBitMask(int order)
{
	int64 mask=3;
	int64 shift = (2*order)-2;
	if (shift == 1){ mask = 3 << shift; }
	else if (shift == 2){ mask = 3 << shift; }
	else if (shift == 3){ mask = 3 << shift; }
	else if (shift == 4){ mask = 3 << shift; }
	else if (shift == 5){ mask = 3 << shift; }
	else if (shift == 6){ mask = 3 << shift; }
	else if (shift == 7){ mask = 3 << shift; }
	else if (shift == 8){ mask = 3 << shift; }
	else if (shift == 9){ mask = 3 << shift; }
	else if (shift == 10){ mask = 3 << shift; }
	else if (shift == 11){ mask = 3 << shift; }
	else if (shift == 12){ mask = 3 << shift; }
	else if (shift == 13){ mask = 3 << shift; }
	else if (shift == 14){ mask = 3 << shift; }
	else if (shift == 15){ mask = 3 << shift; }
	else if (shift == 16){ mask = 3 << shift; }
	else if (shift == 17){ mask = 3 << shift; }
	else if (shift == 18){ mask = 3 << shift; }
	else if (shift == 19){ mask = 3 << shift; }
	else if (shift == 20){ mask = 3 << shift; }
	else if (shift == 21){ mask = 3 << shift; }
	else if (shift == 22){ mask = 3 << shift; }
	else if (shift == 23){ mask = 3 << shift; }
	else if (shift == 24){ mask = 3 << shift; }
	else if (shift == 25){ mask = 3 << shift; }
	else if (shift == 26){ mask = 3 << shift; }
	else if (shift == 27){ mask = 3 << shift; }
	else if (shift == 28){ mask = 3 << shift; }
	else if (shift == 29){ mask = 3 << shift; }
	else if (shift == 30){ mask = 3 << shift; }
	else if (shift == 31){ mask = 3 << 1; mask = mask << 30; }
	else if (shift == 32){ mask = 3 << 1; mask = mask << 31; }
	else if (shift == 33){ mask = 3 << 1; mask = mask << 32; }
	else if (shift == 34){ mask = 3 << 2; mask = mask << 32; }
	else if (shift == 35){ mask = 3 << 3; mask = mask << 32; }
	else if (shift == 36){ mask = 3 << 4; mask = mask << 32; }
	else if (shift == 37){ mask = 3 << 5; mask = mask << 32; }
	else if (shift == 38){ mask = 3 << 6; mask = mask << 32; }
	else if (shift == 39){ mask = 3 << 7; mask = mask << 32; }
	else if (shift == 40){ mask = 3 << 8; mask = mask << 32; }
	else if (shift == 41){ mask = 3 << 9; mask = mask << 32; }
	else if (shift == 42){ mask = 3 << 10; mask = mask << 32; }
	else if (shift == 43){ mask = 3 << 11; mask = mask << 32; }
	else if (shift == 44){ mask = 3 << 12; mask = mask << 32; }
	else if (shift == 45){ mask = 3 << 13; mask = mask << 32; }
	else if (shift == 46){ mask = 3 << 14; mask = mask << 32; }
	else if (shift == 47){ mask = 3 << 15; mask = mask << 32; }
	else if (shift == 48){ mask = 3 << 16; mask = mask << 32; }
	else if (shift == 49){ mask = 3 << 17; mask = mask << 32; }
	else if (shift == 50){ mask = 3 << 18; mask = mask << 32; }
	else if (shift == 51){ mask = 3 << 19; mask = mask << 32; }
	else if (shift == 52){ mask = 3 << 20; mask = mask << 32; }
	else if (shift == 53){ mask = 3 << 21; mask = mask << 32; }
	else if (shift == 54){ mask = 3 << 22; mask = mask << 32; }
	else if (shift == 55){ mask = 3 << 23; mask = mask << 32; }
	else if (shift == 56){ mask = 3 << 24; mask = mask << 32; }
	else if (shift == 57){ mask = 3 << 25; mask = mask << 32; }
	else if (shift == 58){ mask = 3 << 26; mask = mask << 32; }
	else if (shift == 59){ mask = 3 << 27; mask = mask << 32; }
	else if (shift == 60){ mask = 3 << 28; mask = mask << 32; }
	else if (shift == 61){ mask = 3 << 29; mask = mask << 32; }
	else if (shift == 62){ mask = mask << 30; mask = mask << 32; }
	else if (shift == 63){ mask = mask << 31; mask = mask << 32; }
	else if (shift == 64){ mask = mask << 32; mask = mask << 32; }
	return mask;
}

int64 MortonLQT::HpxToMorton(int64 hpxid,int order)
{
	if(VERBOSE)
	{
	  cout << "Input: hpxid = " << hpxid << " of order " << order << endl;
	}
	int64 morton = 0;
	int64 mask = 0;
	int64 nextBit = 0;
	// Process the 2-bit groups
	mask = 3 << ((2*order)-2); 
	for(int i = order; i > 0; i--)  
	{
	  if(VERBOSE) cout << "  mask #" << i <<": " << mask << endl;
	  nextBit = ((hpxid & mask) >> ((2*i)-2) ) + 1; 
	  if(VERBOSE)  cout << "     value: " << nextBit << endl;
	  morton = AppendMortonBit(morton,nextBit);
	  mask = mask >> 2; 
	}
	return morton;	
}


void MortonLQT::MortonToHpx(int64 morton,int64& hpxid, int& order)
{
	int nextBit = 0;
	if(VERBOSE) cout << "Input: Morton = " << morton << endl;
	order = GetMortonLevel(morton);
	if(VERBOSE) printf("  Order: %d\n",order);
	hpxid = 0;
	int temp = 0;
	int shift = 0;
	for(int i = 0; i < order; i++)
	{
	  nextBit = GetMortonBit(morton,order-i)-1;
	  if(VERBOSE) cout << "  Morton digit #" << i << ": " << nextBit << endl;
	  nextBit = nextBit << shift; 
	  if(VERBOSE) cout << "  after shift: " << nextBit << endl;
	  hpxid += nextBit;  
	  if(VERBOSE) cout << "hpxid " << hpxid << endl;
	  shift += 2;
	}
}

int64 MortonLQT::LongLatDegToMorton(float longHpx,float coLatHpx,int order)
{
	return LongLatToMorton(longHpx*D2R,coLatHpx*D2R,order);
}

int64 MortonLQT::LongLatToMorton(float longHpx,float coLatHpx,int order)
{
   pointing pt;
   pt.theta = coLatHpx;
   pt.phi = longHpx;
   hpxQ.Set(order,NEST);
   int64 hpxId = hpxQ.ang2pix(pt);
   int64 morton = HpxToMorton(hpxId,order);
	if (VERBOSE) cout << "HpxId: " << hpxId << " Order: " << order << " Morton: " << morton << endl;
   return morton;
}

void MortonLQT::MortonToLongLat(int64 morton,float& longHpx, float& coLatHpx)
{
   pointing pt;
   int64 hpxId;
   int order;
   MortonToHpx(morton,hpxId,order);
   hpxQ.Set(order,NEST);
   pt = hpxQ.pix2ang(hpxId);
   longHpx = pt.phi;
   coLatHpx = pt.theta;
   if(VERBOSE) cout << "LongHpx: " << pt.phi << " CoLatHpx: " << pt.theta << endl;
}

// Given list of candidate data nodes (all have Morton Code that maps to longitude,latitude)
// Determine which data node closest matches the query Longitude,Latitude. 
std::vector<int64> MortonLQT::ClosestMortonToLongLat(std::vector<MortonNode> NodeList,float qLong,float qLat)
{
   float deltaX,deltaY,deltaNet;
   float minDist = 99999999.0;
   int64 morton = -1;
   int sub = -1;
   std::vector<int64> closest;
   for( unsigned int i = 0; i < NodeList.size(); i++ )
   {
      deltaX = (qLong-NodeList[i].longitude);
      deltaY = (qLat-NodeList[i].latitude);
      deltaNet = deltaX*deltaX + deltaY*deltaY;
      if( deltaNet < minDist )
	  {
         minDist = deltaNet;
         morton = NodeList[i].morton;
         sub = NodeList[i].sub;
   	  }
   }
   closest.push_back(morton);
   closest.push_back(sub);
   return closest;
}



unsigned int MortonLQT::GetNumMortonNodes()
{
  return( mTree.size() );
}

int MortonLQT::GetMortonTreeDepth()
{
  int maxdepth = 0;
  int level = 0;
  // Check for zero size tree
  if( GetNumMortonNodes() == 0)
	  return 0;
  for( unsigned int i = 0; i < GetNumMortonNodes(); i++ )     
     level = GetMortonLevel(mTree[i].morton);
     if( level > maxdepth)
	 {
        maxdepth = level;
	 }
  return( maxdepth ); 
}

void MortonLQT::PrintMortonTree()
{
     cout << endl;
     cout << "######################################################" << endl;
     cout << "#######    MORTON    LINEAR    QUAD    TREE    #######" << endl;
     cout << "######################################################" << endl;
     cout << "Index Morton SubAddr ChildrenYN Longitude Latitude Data" << endl;
     cout << "======================================================" << endl;
     for(unsigned int i = 0; i < GetNumMortonNodes(); i++)
	 {
		   cout << i << " "
			      << mTree[i].facenum << " "
				  << mTree[i].morton << " "
				  << mTree[i].sub << " "
			     << mTree[i].childrenYN << " "
			     << mTree[i].longitude << " "
			     << mTree[i].latitude << " "
			     << mTree[i].data << endl;

	 }
     cout << "======================================================" << endl;
     cout << "Number of nodes: " << GetNumMortonNodes() << endl;
     cout << "Maximum depth of tree: " << GetMortonTreeDepth() << endl;
}


int MortonLQT::FindIndexAtMorton(int64 morton)
{
  return( FindIndexAtMortonSub_internal(morton,1) );
}

int MortonLQT::FindIndexAtMortonSub(int64 morton,int sub)
{
  return( FindIndexAtMortonSub_internal(morton,sub) );
}

// Simple Linear Search of sTree vector until matching Morton,Sub combination is found. Returns index of sTree where
// match occured. If no match is found return -1.
int MortonLQT::FindIndexAtMortonSub_internal( int64 sMorton, int sSub)
{
   for(unsigned int i = 0; i < mTree.size(); i++ )
   {
      if( (mTree[i].morton == sMorton && mTree[i].sub == sSub) )
	  {
         return i;
	  }
   }
   return -1;
}


std::vector<MortonNode> MortonLQT::SearchMortonNodeAtLongLat(float longitude,float latitude)
{
  std::vector<int64> morton_sub = FindMortonAtLongLat(longitude,latitude);
  return SearchMortonNode(morton_sub[0],morton_sub[1]);
}

std::vector<MortonNode> MortonLQT::SearchMortonNodeHpx(int64 hpxid,int order)
{
	int64 morton = HpxToMorton(hpxid,order);
	return SearchMortonNode(morton,-1);
}


std::vector<MortonNode> MortonLQT::SearchMortonNode(int morton,int sub)
{
	std::vector<MortonNode> found_nodes;
	// If specified root node then search entire tree
	if( sub == -1 )
	{
		SearchMortonNode_internal(found_nodes,1,sub);
		SearchMortonNode_internal(found_nodes,2,sub);
		SearchMortonNode_internal(found_nodes,3,sub);
		SearchMortonNode_internal(found_nodes,4,sub);
	}
	else
	{
	   SearchMortonNode_internal(found_nodes,morton,sub);
	}
	return( found_nodes );
}

void MortonLQT::SearchMortonNode_internal(std::vector<MortonNode>& found_nodes,int64 morton, int sub)
{
	// If Node has children, descend to them while keeping list of found data nodes
	// along the way. Return when no children are found.
	std::vector<MortonNode> nodes = GetNodeAtMorton(morton,sub);
	  
	for(unsigned int i = 0; i < nodes.size(); i++)
	{
	   if( nodes[i].data != -1 ) 
	   {
		  found_nodes.push_back(nodes[i]);
	   }
	}
	if( nodes[0].childrenYN == 1 )
	{
		// Compute children of node
		std::vector<int64> children = ChildrenOfMorton(morton);	   
		SearchMortonNode_internal(found_nodes,children[0],sub);
		SearchMortonNode_internal(found_nodes,children[1],sub);
		SearchMortonNode_internal(found_nodes,children[2],sub);
		SearchMortonNode_internal(found_nodes,children[3],sub);
	}
}


std::vector<int64> MortonLQT::FindMortonAtLongLat(float longitude,float latitude)
{
	int sub = 1;	
	bool done = false;
	int level = 1;
	int64 morton;
	std::vector<MortonNode> node;
	std::vector<int64> morton_sub(2);
	while( done == false )
	{
		morton = LongLatToMorton(longitude,latitude,level);
		morton_sub[0] = morton;
		morton_sub[1] = sub;
		node = GetNodeAtMorton(morton,-1);   
		// Return empty list if nothing found
		if( node.size() == 0 )
			return morton_sub;
		if( node[0].childrenYN == 1 )
		{
		  level += 1;
		}
		else
		{
		  // Possibility of multiple data nodes with duplicate
		  // Morton Code, must examine each duplicate to find
		  // data node with closest matching longitude, latitude
		  if( node.size() > 1 )
		  {
			 morton_sub = ClosestMortonToLongLat(node,longitude,latitude);
		  }
		  done = true;
		}
	}
	return( morton_sub );
}

std::vector<MortonNode> MortonLQT::GetNodeAtMorton(int morton,int sub)
{
	std::vector<MortonNode> found_node_list;
	MortonNode found_node;
	int index;

	// Do we want ALL indices under morton returned OR
	// just the actual morton,sub pair?
	// Return all matching morton if "sub" equals -1
	if( sub == -1 )
	{
		index = FindIndexAtMorton(morton);
	}
	else
	{
		index = FindIndexAtMortonSub(morton,sub);
	}

	// Node NOT Found, return empty list
	if( index == -1 )
	{
		return found_node_list;
	}
	
	// Otherwise, found the Node so populate the returnable
	found_node.morton = mTree[index].morton;
	found_node.sub = mTree[index].sub;
	found_node.longitude = mTree[index].longitude;
	found_node.latitude = mTree[index].latitude;
	found_node.childrenYN = mTree[index].childrenYN;
	found_node.data = mTree[index].data;
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
			  if( mTree[index].morton == morton )
			  {
				 found_node.morton = mTree[index].morton;
				 found_node.sub = mTree[index].sub;
				 found_node.longitude = mTree[index].longitude;
				 found_node.latitude = mTree[index].latitude;
				 found_node.childrenYN = mTree[index].childrenYN;
				 found_node.data = mTree[index].data;
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


void MortonLQT::UpdateNode(MortonNode node)
{
    // Find index of list MortonNode object that contains 
	// Morton code = morton
	int index = FindIndexAtMorton(node.morton);
	if( index != -1 )
	{
		mTree[index].morton = node.morton;
		mTree[index].sub = node.sub;
		mTree[index].childrenYN = node.childrenYN;
		mTree[index].longitude = node.longitude;
		mTree[index].latitude = node.latitude;
		mTree[index].data = node.data;
		mTree[index].facenum = node.facenum;
		if(VERBOSE)
		{
			cout << "Updated node: " << mTree[index].facenum << " "
									  << mTree[index].morton << " "
				                      << mTree[index].sub << " "
											 << mTree[index].childrenYN << " "
											 << mTree[index].longitude << " "
											 << mTree[index].latitude << " "
											 << mTree[index].data << endl;
		}
	}
}

void MortonLQT::AppendNode(MortonNode node)
{
    // Must compute the sub of this duplicate morton code
    // Subaddress is based on first come first served
	std::vector<MortonNode> found_nodes;
	found_nodes = GetNodeAtMorton(node.morton,-1);

	// Want to know the highest sub-address so take the last MortonNode in the found_nodes
	// list and return it and use its "sub" address.
	int sub = found_nodes.back().sub;
      
    MortonNode new_node;
	new_node.morton = node.morton;
	new_node.sub = sub+1;
	new_node.childrenYN = node.childrenYN;
	new_node.longitude = node.longitude;
	new_node.latitude = node.latitude;
	new_node.data = node.data;
	new_node.facenum = node.facenum;
	mTree.push_back(new_node);
	if(VERBOSE)
	{
		cout << "Append duplicate new node: " << new_node.facenum << " "
														  << new_node.morton << " "
														  << new_node.sub << " "
														  << new_node.childrenYN << " "
														  << new_node.longitude << " "
													     << new_node.latitude << " "
												 	     << new_node.data << endl;   
	}
	// Sort the tree
	std::sort( mTree.begin(), mTree.end(), SortFunctionMorton);
}

void MortonLQT::AddNode(MortonNode node)
{
	MortonNode new_node;
	new_node.morton = node.morton;
	new_node.sub = node.sub;
	new_node.childrenYN = node.childrenYN;
	new_node.longitude = node.longitude;
	new_node.latitude = node.latitude;
	new_node.data = node.data;
	new_node.facenum = node.facenum;
	mTree.push_back(new_node);
	if(VERBOSE)
	{
		cout << "Added new node: " << new_node.facenum << " "
									<< new_node.morton << " "
											<< new_node.sub << " "
											<<	new_node.childrenYN << " "
											<< new_node.longitude << " "
											<< new_node.latitude << " "
											<< new_node.data << endl;
	}
	// Now must add the other sibling nodes at same level of quad tree
	std::vector<int64> siblings = SiblingsOfMorton(node.morton);     
	if(VERBOSE)
	{
		cout << "Create siblings of " << node.morton << " "
												<< siblings[0] << " "
												<< siblings[1] << " "
												<< siblings[2] << endl;
	}
	MortonNode empty_node;
	empty_node.morton = siblings[0];
	mTree.push_back(empty_node);
	if(VERBOSE)
	{
		cout << "Added sibling node: " << empty_node.facenum << " "
										<< empty_node.morton << " "
												 << empty_node.sub << " "
												 << empty_node.childrenYN << " "
												 << empty_node.longitude << " "
												 << empty_node.latitude << " "
												 << empty_node.data << endl;
	}
	empty_node.morton = siblings[1];
	mTree.push_back(empty_node);
	if(VERBOSE)
	{
		cout << "Added sibling node: " << empty_node.facenum << " "
												<< empty_node.morton << " "
												 << empty_node.sub << " "
												 << empty_node.childrenYN << " "
												 << empty_node.longitude << " "
												 << empty_node.latitude << " "
												 << empty_node.data << endl;
	}
	empty_node.morton = siblings[2];
	mTree.push_back(empty_node);
	if(VERBOSE)
	{
		cout << "Added sibling node: " << empty_node.facenum << " "
												<< empty_node.morton << " "
												 << empty_node.sub << " "
												 << empty_node.childrenYN << " "
												 << empty_node.longitude << " "
												 << empty_node.latitude << " "
												 << empty_node.data << endl;
	}
	// Sort the tree
	std::sort( mTree.begin(), mTree.end(), SortFunctionMorton);

	// Lastly, make sure parent (if exists) of new nodes is marked as having children
	if( ParentOfMorton(node.morton) != -1 )
	{
		MortonNode update_node;
		update_node.morton = ParentOfMorton(node.morton);
		update_node.longitude = 0.0;
		update_node.longitude = 0.0;
		update_node.childrenYN = 1;
		update_node.data = -1;
		UpdateNode(update_node);
	}
}

int64 MortonLQT::InsertMortonNodeAtMorton(int64 morton,int data)
{
	// Convert Morton code to longitude & latitude
	// then insert
	float longHpx,coLatHpx;
	MortonToLongLat(morton,longHpx,coLatHpx);
	return( InsertMortonNodeAtLongLat(longHpx,coLatHpx,data) );
}

int64 MortonLQT::InsertMortonNodeAtHpx(int64 hpxid,int order,int data)
{
	// Convert HEAPix index to Morton code then insert via InsertNodeMorton
	int64 morton = HpxToMorton(hpxid,order);
	return( InsertMortonNodeAtMorton(morton,data) );
}

int64 MortonLQT::InsertMortonNodeAtLongLat(float longitude,float latitude,int data)
{
	return InsertNode_internal_1(longitude,latitude,data);
}

int64 MortonLQT::InsertNode_internal_1(float longitude,float latitude,int data)
{
	int treeDepth = 1;
	int64 new_morton;
	bool done = false;
	float longitude_save,latitude_save;
	int childrenYN_save;
	int data_save;
	// Create new Morton Node to be inserted based on its longitude & latitude
	MortonNode insert_node;
	int64 morton = LongLatToMorton(longitude,latitude,treeDepth);
	insert_node.morton = morton;
	insert_node.longitude = longitude;
	insert_node.latitude = latitude;
	insert_node.data = data;
	insert_node.childrenYN = 0;

	// Insert the new Morton node checking for data collision (Morton
	// node who's "data" element is NOT equal to -1. When collision
	// occurs create new MortonNode's who's Morton codes are the children
	// of the Morton code to be inserted. We recompute the Morton code to 
	// next level deeper in tree so will match one of the computed child
	// Morton codes. If there is no data collision for matched Morton node
	// then we simply update the matched Morton node with data to be inserted.
	// Last possibility is that the Morton node is not found in the list. If
	// so we create new Morton node and add it to the list. In all cases of 
	// inserted Morton nodes we finish by sorted the list by ascending Morton
	// code
	while( done == false )
	{
       // Search for Morton node in tree, node returned if found, empty vector otherwise
		 std::vector<MortonNode> node = GetNodeAtMorton(insert_node.morton,1);
 
       // If Morton node is found we then need to check for data collision, else update
       // Morton node with data
       if( node.size() > 0 )
		 {
         if(VERBOSE)
			{
           cout << "Morton Node " << node[0].morton << " found!" << endl;
			}
         // If Morton node has children then need to re-calculate Morton
         // code of node to be inserted to next level deeper.
         if( node[0].childrenYN == 1 )
			{
            if(VERBOSE)
			   {
               cout << "Found Morton Node " << node[0].morton << " has children." << endl;
			   }
            treeDepth += 1;
            new_morton = LongLatToMorton(longitude,latitude,treeDepth);
            insert_node.morton = new_morton;  
			}
         // Otherwise check to see if there is a data node insertion collision or not
         else
			{
            // If "data" is populated in MortonNode (not "None") then
            // we have a node insertion collision!
            if( node[0].data != -1 )
			   {
              if(VERBOSE)
				  {
                cout << "Collision! Morton Node " << node[0].morton << " already has data: " << node[0].data << endl;
				  }
              // Store collided Morton Node's information for re-insertion
              longitude_save = node[0].longitude;
              latitude_save = node[0].latitude;
              childrenYN_save = node[0].childrenYN;
              data_save = node[0].data;
                  
              // Re-calculate next level deeper Morton code for Morton Node to be inserted. 
              treeDepth += 1;
                  
              // Check if have reached user specified maximum tree depth.
              // If so will create a NEW Node with SAME Morton Code and add to the tree.
              if( treeDepth > userMaxDepth )
				  {
                 if(VERBOSE)
					  {
                    cout << "Reached Max Tree Depth: Append Node with identical Morton Code: " << insert_node.morton << endl;
                    cout << "Specify deeper Max Tree Depth!" << endl;
					  }
                 AppendNode(insert_node);
                 done = true;
			     }
              else
				  {
                   new_morton = LongLatToMorton(insert_node.longitude,insert_node.latitude,treeDepth);
                   insert_node.morton = new_morton;
                     
                   // Add children Morton nodes to Morton list, populate the matching
                   // child Morton == re-calculated Morton code node with "data" to be
                   // inserted, other child Morton nodes have "data" == "None"
                   AddNode(insert_node);	
                   // Now must prepare the collided node to move to new node at next 
                   // lower level of the tree.
                   new_morton = LongLatToMorton(longitude_save,latitude_save,treeDepth);
                   insert_node.morton = new_morton;
                   insert_node.longitude = longitude_save;
                   insert_node.latitude = latitude_save;
                   insert_node.childrenYN = childrenYN_save;
                   insert_node.data = data_save;
                   if(VERBOSE)
					    {
                      cout << "Now re-insert the collided node with data " << data_save << " to " << new_morton << endl;
					    }
				  }
			   }
            // Otherwise no insertion collision so update the data,longitude,latitude,etc
            // attributes of MortonNode.
            else
			   {
              if(VERBOSE)
				  {
                cout << "No insertion collision found, updating Morton Node: " << node[0].morton << " with " << insert_node.data << endl;
				  }
              UpdateNode(insert_node);
              if( treeDepth > curTreeDepth )
				  {
                curTreeDepth = treeDepth;
				  }
              done = true;
				}
			}
		 }            
       else
		 {
         if(VERBOSE)
		   {
           cout << "Morton Node " << insert_node.morton << " not found, Adding new Morton Node!" << endl;
			}
         // If Morton code not found in list then append it to the end
         // of the Morton list then sort the list ascending by Morton code
         AddNode(insert_node);
         if( treeDepth > curTreeDepth )
			{
           curTreeDepth = treeDepth;
			}
         done = true; 
		 }
	  }      
     return( insert_node.morton );
}


void MortonLQT::DeleteMortonNodeAtLongLat(float longitude,float latitude)
{
	std::vector<int64> morton_sub = FindMortonAtLongLat(longitude,latitude);
	DeleteMortonNodeAtMorton(morton_sub[0],morton_sub[1]);
}

void MortonLQT::DeleteMortonNodeAtHpx(int64 hpxid,int order)
{
	int64 morton = HpxToMorton(hpxid,order);			      
	DeleteMortonNodeAtMorton(morton,1);      
}

void MortonLQT::DeleteMortonNodeAtMorton(int64 morton,int sub)
{
	// Find the tree index of the node to be deleted
	int index = FindIndexAtMortonSub(morton,sub);

	if( index == -1 )
	{
		return;
	}
	if(VERBOSE)
	{
	   cout << "DeleteNodeAtMorton: " << morton << " "
			                            << sub << " "
												 << index << endl;
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
		if( oldTree[i].data != -1 )
		{
			InsertMortonNodeAtLongLat(oldTree[i].longitude,oldTree[i].latitude,oldTree[i].data);
		}
	}
}

int MortonLQT::SaveTreeToFile(std::string filename)
{
	FILE* fp = fopen(filename.c_str(), "w");
	
	if(fp == NULL)
	{
	  cout << "Unable to open MLQT output file: " << filename.c_str() << endl;
	  return false;        
	}
	// Write out header
	WriteHeader(fp);

	// Write out data
	WriteData(fp);	 

	return true;
}

int MortonLQT::LoadTreeFromFile(std::string filename)
{
	  char buffer[1024],nextLine[1024];
	  long int res;
	  int numNodes,nodeNumber,facenum, morton,sub,data;
	  int childrenYN;
	  float longitude,latitude;
	  ifstream fp;
	  fp.open(filename.c_str());
      
      if(fp.fail())
	  {
          cout << "Unable to open MLQT output file: " << filename.c_str() << endl;
          return false;     
	  }
          
      // Reset current tree
      mTree.clear();    
          
      // Skip comment line
	  skip_line(fp);
            
      // Next line is number of records
	  fp >> numNodes;  skip_line(fp);


      // Skip header line
	  skip_line(fp);

      // Next lines are data records (the nodes of the tree)
      for(int i = 0; i < numNodes; i++)
	  {
	     fp >> nodeNumber >> facenum >> morton >> sub >> childrenYN >> longitude >> latitude >> data;
         
		 // Skip the rest of the line
		 skip_line(fp);

      	 MortonNode next_node;
		 next_node.facenum = facenum;
         next_node.morton = morton;
         next_node.sub = sub;
		 next_node.childrenYN = childrenYN;
         next_node.longitude = longitude;
         next_node.latitude = latitude;
         next_node.data = data;
         mTree.push_back(next_node);
	  }
      fp.close();
	  return true;
}


void MortonLQT::WriteHeader(FILE* fp)
{
     // Write out creation origin, date, timestamp
	 time_t rawtime;
	 struct tm* timestamp;
	 time(&rawtime);
	 timestamp = localtime(&rawtime);
     fprintf(fp,"*** PyMortonLQT File Created on %s ***\n",timestamp);      
      
     // Write out count of total number of nodes (records in file)
     fprintf(fp,"%d\n",GetNumMortonNodes());
  
     // Write out data record tags
     fprintf(fp,"INDEX\t");
	 fprintf(fp, "FACENUM\t");
     fprintf(fp,"MORTON\t");
     fprintf(fp,"SUB\t");
     fprintf(fp,"CHILDRENYN\t");
     fprintf(fp,"LONGITUDE\t");
     fprintf(fp,"LATITUDE\t");
     fprintf(fp,"DAT\t");
     fprintf(fp,"\n");
}

void MortonLQT::WriteData(FILE* fp)
{
      for(int i = 0; i < GetNumMortonNodes(); i++)
	  {
      	 fprintf(fp,"%d\t",i); // Index
		 fprintf(fp, "%d\t", mTree[i].facenum); // Facenum
      	 fprintf(fp,"%ld\t",mTree[i].morton); // Morton
      	 fprintf(fp,"%d\t",mTree[i].sub); // Sub Index
      	 fprintf(fp,"%d\t",mTree[i].childrenYN); // Children YN
      	 fprintf(fp,"%f\t",mTree[i].longitude); // Longitude
      	 fprintf(fp,"%f\t",mTree[i].latitude); // Latitude
      	 fprintf(fp,"%d\t",mTree[i].data); // Data Index
         fprintf(fp,"\n");
	  }
      fclose(fp);   	
}