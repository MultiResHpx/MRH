#ifndef MULTIRESHPX_MAP_H
#define MULTIRESHPX_MAP_H


//#define BENCHMARKING

#include "MultiResHpx.h"

/*! A MultiResHpx map of a given datatype */
template<typename T> class MultiResHpx_Map: public MultiResHpx
{
 
  private:
	  std::vector<T> map;
      int64 numrec;

  public:
    /*! Constructs an unallocated map. */
	  MultiResHpx_Map () 
		  : MultiResHpx(MAX_HPX_ORDER64,scheme) {numrec = 0; Healpix_Ordering_Scheme scheme = NEST; map.clear(); }

    /*! Constructs a map with a given \a order ordering
        scheme \a scheme and \a overwrite overwrite node permission */
    MultiResHpx_Map (int max_depth, Healpix_Ordering_Scheme scheme)
      : MultiResHpx ( max_depth, scheme) 
	{ 
		numrec = 0;
		map.clear();
	}
       
    /*! Adds \a T template type record to map and associates it with pointing vector
	    location which will produce unique HPX ID as key for record in map */
	void AddRecord (T data, pointing pt);

	void AddRecord2 (T data, pointing pt);

	int SaveMapToArchive(std::string filename);

	int LoadMapFromArchive(std::string filename);

	int64 NumRec();

	int64 SizeInMemKb();

    /*! Returns a constant reference to the map data. */
	const std::vector<T> &Map() const;

	/*! Clear map. */
	inline void ClearMap() { numrec = 0; map.clear(); }

	/*! Adds record to map. */
	void AddMapRecord(T data);

    /*! Returns a constant reference to the map element with \a index. */
    const T &operator[] (int index) const { return map[index]; }
    /*! Returns a reference to the map element with \a index. */
    T &operator[] (int index) { return map[index]; }

	/*! Search */
	std::vector<T> Search (pointing pt);
	std::vector<T> Search (int64 hpxid, int order, bool upsearch);

	/*! Range Queries */
    std::vector<T> QueryDisc ( pointing pt, double radius );

	std::vector<T> QueryPolygon ( std::vector<pointing> poly );

    std::vector<T> QueryStrip ( double theta1, double theta2 );

	std::vector<T> Neighbors( pointing pt, int64 order );

    std::vector<pair<T,T>> TwoPointCorrBin( double radius );

	// Insertion operator
	//friend std::ostream& operator<<(std::ostream& os, const MultiResHpx_Map& s)
	//{
	//	os << s.numrec << '\n';
	//	for(int i = 0; i < s.numrec; i++) {
	//		os << s.map[i];
	//	}
	//	return os;
	//}
	// Extraction operator
	//friend std::istream& operator>>(std::istream& is, MultiResHpx_Map& s)
	//{
	//	// read in individual members of s
	//	is >> s.numrec >> s.map;
	//	return is;
	//}

};


template<typename T> inline void MultiResHpx_Map<T>::AddRecord2 
(
 T record, 
 pointing pt
)
{
    // !!!!BENCHMARKING MODIFICATION!!! 
	// Since HEALPix_Map allows for overwriting map records when multiple data points map to same 
	// HEALPix index, we need to be able to mimic that behavior with MRH_Map. So BEFORE we add new
	// map record we FIRST determine if the new map record's spatial coordinates are already found
	// in a previously inserted map record. If so, we'll replace the previously inserted map record
	// with this new one. Otherewise, just append the new map record as always.
	bool overwrite = false;
	std::vector<MortonNode> found = MultiResHpx::Search(pt);

	// If found pt in MultiResHpx data structure:
	// 1. Are spatial coordinates IDENTICAL to within machine epsilon? OR
	// 2. Is Morton Code of "found" at maximum resolution? 
	//
	// If either #1 or #2 are true then new pt (and data) will overwrite existing map
	// element and reference in MultiResHpx data structure.
	//
	// Otherwise, don't overwrite the map record and proceed as normal, append new map record
	// and append spatial coordinate reference into MultiResHpx data structure
	if(found.size() > 0 ) {

		// Are spatial coordinates of "found" and "pt" IDENTICAL to within machine epsilon?
		if( abs_approx(pt.phi,found[0].phi) && abs_approx(pt.theta,found[0].theta) ) {
			overwrite = true;
		}

		// Is Morton Code of "found" at maximum resolution 
		if( found[0].m.LEVEL == MultiResHpx::GetTreeDepthLimit() ) {
			overwrite = true;
		}
	}

	if( overwrite == false )
	{
		// Add element to map
		map.push_back(record);

		// Add map record index to MRH data structure as leaf node
		MultiResHpx::Insert2( pt, numrec );

		// Increment number of map records
		numrec += 1;
	} else {
		// Overwrite existing record with new record
		map[found[0].data[0]] = record;
		MultiResHpx::Insert2( pt, found[0].data[0] );
	}

}


/*! Adds \a T template type record to map and associates it with pointing vector
    location which will produce unique Morton code as key for record in map */
template<typename T> inline void MultiResHpx_Map<T>::AddRecord 
(
 T record, 
 pointing pt
)
{

    // !!!!BENCHMARKING MODIFICATION!!! 
	// Since HEALPix_Map allows for overwriting map records when multiple data points map to same 
	// HEALPix index, we need to be able to mimic that behavior with MRH_Map. So BEFORE we add new
	// map record we FIRST determine if the new map record's spatial coordinates are already found
	// in a previously inserted map record. If so, we'll replace the previously inserted map record
	// with this new one. Otherewise, just append the new map record as always.
#ifdef BENCHMARKING
	bool overwrite = false;
	std::vector<MortonNode> found = MultiResHpx::Search(pt);

	// If found pt in MultiResHpx data structure:
	// 1. Are spatial coordinates IDENTICAL to within machine epsilon? OR
	// 2. Is Morton Code of "found" at maximum resolution? 
	//
	// If either #1 or #2 are true then new pt (and data) will overwrite existing map
	// element and reference in MultiResHpx data structure.
	//
	// Otherwise, don't overwrite the map record and proceed as normal, append new map record
	// and insert spatial coordinate reference into MultiResHpx data structure
	if(found.size() > 0 ) {

		// Are spatial coordinates of "found" and "pt" IDENTICAL to within machine epsilon?
		if( abs_approx(pt.phi,found[0].phi) && abs_approx(pt.theta,found[0].theta) ) {
			overwrite = true;
		}

		// Is Morton Code of "found" at maximum resolution 
		if( found[0].m.LEVEL == MultiResHpx::GetTreeDepthLimit() ) {
			overwrite = true;
		}

	} 
	
	if( overwrite == false )
	{
		// Add element to map
		map.push_back(record);

		// Add map record index to MRH data structure as leaf node
		MultiResHpx::Insert( pt, numrec );

		// Increment number of map records
		numrec += 1;
	} else {

		// Overwrite existing record with new record
		map[found[0].data[0]] = record;
		MultiResHpx::Insert( pt, found[0].data[0] );

	}
#else
	// Add data record to map
	map.push_back(record);

	// Add map record index to MRH data structure as leaf node
	MultiResHpx::Insert( pt, numrec );

	// Increment number of map records
	numrec += 1;
#endif
}

template<typename T> 
inline void MultiResHpx_Map<T>::AddMapRecord
(
 T record
)
{
	map.push_back(record);
	numrec += 1;
}

template<typename T> 
inline int MultiResHpx_Map<T>::SaveMapToArchive
(
 std::string filename
)
{
	std::ofstream fs;

	fs.open(filename.c_str(),ios::binary);

	// First write out number of records
	fs << map.size();

	// Now write out each template record
	for(unsigned int i = 0; i < map.size(); i++)
	{
		fs.write((char *)&map[i],sizeof(T));
	}
	fs.close();
	return 0;
}

template<typename T>
inline int MultiResHpx_Map<T>::LoadMapFromArchive
(
 std::string filename
)
{
	std::ifstream fs;
	fs.open(filename.c_str(),ios::binary);


	// First record is the number of records
	fs >> numrec;

	// Initialize MAP to hold num_rec records
	map.resize(numrec);

	// Now read in template records
	for(unsigned int i = 0; i < numrec; i++)
	{
		fs.read((char *)&map[i],sizeof(T));
	}
	fs.close();
	return 0;
}



/*! Returns a constant reference to the map data. */
template<typename T> inline const std::vector<T>& MultiResHpx_Map<T>::Map() const 
{ 
	return map; 
}


template<typename T> inline int64 MultiResHpx_Map<T>::NumRec()
{
	return numrec;
}

template<typename T> inline int64 MultiResHpx_Map<T>::SizeInMemKb()
{
	//Compute the total memory footprint of MRH_Map + MRH data structure
	int64 SizeOfMap = numrec*sizeof(T)/1024;  //Convert to Kb
	int64 SizeOfMortonLQT = MemSizeKb();

	return SizeOfMap+SizeOfMortonLQT;

}

/*! Returns all T data that maps to spatial coordinates specified */
template<typename T> std::vector<T> MultiResHpx_Map<T>::Search (pointing pt)
{
	std::vector<MortonNode> foundIDX;
	std::vector<T> foundT;
	foundIDX = MultiResHpx::Search (pt);
	for( unsigned int i = 0; i < foundIDX.size(); i++ ) {
		for( unsigned int j =0; j < foundIDX[i].data.size(); j++ ) {
			foundT.push_back( map[ foundIDX[i].data[j] ] );
		}
	}
	return foundT;
}

/*! Returns all T data that maps to spatial coordinates specified */
template<typename T> std::vector<T> MultiResHpx_Map<T>::Search (int64 hpxid, int order, bool upsearch)
{
	std::vector<MortonNode> foundIDX;
	std::vector<T> foundT;
	foundIDX = MultiResHpx::Search (hpxid,order,upsearch);
	for( unsigned int i = 0; i < foundIDX.size(); i++ ) {
		for( unsigned int j =0; j < foundIDX[i].data.size(); j++ ) {
			foundT.push_back( map[ foundIDX[i].data[j] ] );
		}
	}
	return foundT;
}



/*! Returns all T data whose spatial coordinates fall within specified query disc */
template<typename T> inline std::vector<T> MultiResHpx_Map<T>::QueryDisc (pointing pt, double radius)
{
	unsigned int i,j;
	std::vector<MortonNode> foundIDX;
	std::vector<T> foundT;
	foundIDX.clear(); foundT.clear();
	foundIDX = MultiResHpx::QueryDisc (pt, radius);
	for( i = 0; i < foundIDX.size(); i++ ) {
		for( j = 0; j < foundIDX[i].data.size(); j++ ) {
			foundT.push_back( map[ foundIDX[i].data[j] ] );
		}
	}
	return foundT;
}

/*! Returns all T data whose spatial coordinates fall within specified query polygon */
template<typename T> std::vector<T> MultiResHpx_Map<T>::QueryPolygon ( std::vector<pointing> poly)
{
	std::vector<MortonNode> foundIDX;
	unsigned int i,j;
	std::vector<T> foundT;
	foundIDX.clear(); foundT.clear();
	foundIDX = MultiResHpx::QueryPolygon (poly);
	for( i = 0; i < foundIDX.size(); i++ ) {
		for(  j = 0; j < foundIDX[i].data.size(); j++ ) {
			foundT.push_back( map[ foundIDX[i].data[j] ] );
		}
	}
	return foundT;
}

/*! Returns all T data whose spatial coordinates fall within specified query latitude strip */
template<typename T> std::vector<T> MultiResHpx_Map<T>::QueryStrip ( double theta1, double theta2)
{
	std::vector<MortonNode> foundIDX;
	std::vector<T> foundT;
	foundIDX = MultiResHpx::QueryStrip (theta1, theta2);
	for( unsigned int i = 0; i < foundIDX.size(); i++ ) {
		for( unsigned int j =0; j < foundIDX[i].data.size(); j++ ) {
			foundT.push_back( map[ foundIDX[i].data[j] ] );
		}
	}
	return foundT;
}

/*! Returns all T data whose spatial coordinates map to neighbor cells of specified spatial coordinates */
template<typename T> std::vector<T> MultiResHpx_Map<T>::Neighbors( pointing pt, int64 order )
{
	std::vector<MortonNode> foundIDX;
	std::vector<T> foundT;
	foundIDX = MultiResHpx::Neighbors (pt,order);
	for( unsigned int i = 0; i < foundIDX.size(); i++ ) {
		for( unsigned int j =0; j < foundIDX[i].data.size(); j++ ) {
			foundT.push_back( map[ foundIDX[i].data[j] ] );
		}
	}
	return foundT;
}

template<typename T> std::vector<pair<T,T>> MultiResHpx_Map<T>::TwoPointCorrBin( double radius )
{
	std::vector<std::pair<MortonNode,MortonNode>> foundIDX;
	std::vector<pair<T,T>> foundTT;
	pair<T,T> nPair;
	foundIDX = MultiResHpx::TwoPointCorrBin (radius);
	for( unsigned int i = 0; i < foundIDX.size(); i++ ) {
		for( unsigned int j =0; j < foundIDX[i].first.data.size(); j++ ) {
			nPair.first = map[ foundIDX[i].first.data[j] ];
			nPair.second = map[ foundIDX[i].second.data[j] ];
			foundTT.push_back(nPair);
			//foundTT.push_back( map[ foundIDX[i].data[j] ] );
		}
	}
	return foundTT;
}
#endif
