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

#include "MultiResHpx_tables.h"
#include <math.h>

const double cos45 = cos(degr2rad*45.0);

const long int order_to_nside[] = { 1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576,
 4194304, 16777216, 67108864, 268435456, 1073741824, 4294967296, 17179869184, 68719476736, 274877906944,
 1099511627776, 4398046511104, 17592186044416, 70368744177664, 281474976710656, 1125899906842624,
 4503599627370496, 18014398509481984, 72057594037927936, 288230376151711744 };

const float basex_to_base0[13][2] = { 
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


