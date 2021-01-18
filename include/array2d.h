/*
 *  Copyright (C) 2007  Simon Perreault
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ARRAY2D_H
#define ARRAY2D_H

#include "array.h"

template< class T >
class Array2D : public Array<T,2>
{
public:
    Array2D();
    Array2D( int m, int n );
    Array2D( const Array<T,2>& array );

    int M() const;
    int N() const;

    const T& at( int i, int j ) const;
    const T& operator() ( int i, int j ) const;
    T& operator() ( int i, int j );

    Array2D<T> subarray( int iBegin, int jBegin, int iEnd, int jEnd );
};

template< class T >
Array2D<T>::Array2D()
    : Array<T,2>()
{
}

template< class T >
Array2D<T>::Array2D( int m, int n )
    : Array<T,2>()
{
    TinyVector<int,2> sizes;
    sizes(1) = m;
    sizes(0) = n;

    Array<T,2>::operator= ( Array<T,2>(sizes) );
}

template< class T >
Array2D<T>::Array2D( const Array<T,2>& array )
    : Array<T,2>(array)
{
}

template< class T >
int Array2D<T>::M() const
{
    return this->sizes()(1);
}

template< class T >
int Array2D<T>::N() const
{
    return this->sizes()(0);
}

template< class T >
const T& Array2D<T>::at( int i, int j ) const
{
    TinyVector<int,2> indices;
    indices(1) = i;
    indices(0) = j;
    return Array<T,2>::at(indices);
}

template< class T >
const T& Array2D<T>::operator() ( int i, int j ) const
{
    return at(i,j);
}

template< class T >
T& Array2D<T>::operator() ( int i, int j )
{
    TinyVector<int,2> indices;
    indices(1) = i;
    indices(0) = j;
    return Array<T,2>::operator()(indices);
}

template< class T >
Array2D<T> Array2D<T>::subarray( int iBegin, int jBegin, int iEnd, int jEnd )
{
    TinyVector<int,2> begin, end;
    begin(1) = iBegin;
    begin(0) = jBegin;
    end(1) = iEnd;
    end(0) = jEnd;
    return Array<T,2>::subarray( begin, end );
}

#endif
