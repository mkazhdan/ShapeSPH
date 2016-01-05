/*
Copyright (c) 2013, Michael Kazhdan
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/
#ifndef CUBE_GRID_INCLUDED
#define CUBE_GRID_INCLUDED

#include "SignalProcessing/SphericalGrid.h"

#ifndef PI
#define PI 3.1415926535897932384
#endif 

// This class represents a set of regular samples of a real-valued, periodic, function in 2D
template<class Real=float>
class CubeGrid{
protected:
	Real* values;
	int res;
public:
	CubeGrid(void);
	CubeGrid( int res );
	~CubeGrid(void);

	// Returns the dimension of the array
	int resolution(void) const;
	// Allocates memory for the array
	int resize(const int& resolution);

	// Clears the values of the array to 0
	void clear(void);

	// Returns a reference to the indexed  array element
	Real& operator() (const int& x,const int& y,const int& z);
	Real operator() (const int& x,const int& y,const int& z) const;
	Real* operator[] ( int x );
	const Real* operator[] ( int x ) const;
	// Returns the linear interpolation of the value at the spedified index
	Real operator() (const double& x,const double& y,const double& z) const;

	// Returns the square of the L2-norm of the array elements
	Real squareNorm(void) const;

	// Reads in an array from the specified file
	int read(const char* fileName);
	int read(FILE* fp);

	// Writes out the array to the specified file
	int write(const char* fileName) const;
	int write(FILE* fp) const;

	// Returns the square of the L2-difference between two arrays
	static Real SquareDifference(const CubeGrid& g1,const CubeGrid& g2);

	// Returns the dot-product of two arrays
	static Real Dot(const CubeGrid& g1,const CubeGrid& g2);

	// Samples the grid over the specified sphere
	void SphereSample( const Real* center , Real radius , SphericalGrid< Real >& sGrid , int threads=1 ) const;
	// Samples the grid over the specified sphere, performing uniform sampling to approximate the monte-carlo integral
	// over the cell dual to the sphere vertices.
	void SphereSample( const Real* center , Real radius , SphericalGrid< Real >& sGrid , int subRes , Real thickness , int threads=1 ) const;
};
#include "CubeGrid.inl"
#endif // CUBE_GRID_INCLUDED
