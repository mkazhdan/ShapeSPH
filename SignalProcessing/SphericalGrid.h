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
#ifndef SPHERICAL_GRID_INCLUDED
#define SPHERICAL_GRID_INCLUDED

#ifndef PI
#define PI 3.1415926535897932384
#endif 

// This templated class represents a spherical function, where index (i,j) corresponds to the
// point on the sphere with spherical coordinates:
//			theta =2PI*i/r
//			phi   = PI*(2j+1)/(2r)
//			(theta,phi) -> [cos(theta)*sin(phi),cos(phi),sin(theta)*sin(phi)]
// where r is the resolution of the sampling.
template<class Real=float>
class SphericalGrid
{
protected:
	Real* values;
	int res;
	static void Transform(const Real matrix[3][3],const Real* in,Real* out);
public:
	SphericalGrid(void);
	SphericalGrid( int );
	~SphericalGrid(void);

	// Returns the dimension of the array
	int resolution(void) const;
	// Allocates memory for the array
	int resize(const int& resolution);

	// Clears the values of the array to 0
	void clear(void);

	// Returns a reference to the indexed array element
	Real& operator() ( const int& x , const int& y );
	Real  operator() ( const int& x , const int& y ) const;
	Real* operator[] ( const int& x );
	// Returns the linear interpolation of the value at the spedified index
	Real operator() (const double& x,const double& y);
	Real operator() (const double& x,const double& y) const ;

	// Returns the square of the L2-norm of the array elements
	Real squareNorm(void) const;

	// Reads in an array from the specified file
	int read(const char* fileName);
	int read(FILE* fp);

	// Writes out the array to the specified file
	int write(const char* fileName) const;
	int write(FILE* fp) const;

	// Sets the (x,y,z) coordinates of the spherical point indexed by (theta,phi)
	static void SetCoordinates(const Real& theta,const Real& phi,Real coords[3]);
	// Sets the (theta,phi) coordinates of the spherical point (x,y,z)
	static void SetCoordinates(const Real coords[3],Real& theta,Real& phi);

	// Sets the (x,y,z) coordinates of the spherical point indexed by (i,j)
	void setCoordinates(const int& i,const int& j,Real coords[3]) const;
	void setCoordinates(const Real& i,const Real& j,Real coords[3]) const;
	// Sets the (i,j) coordinates of the spherical point (x,y,z)
	void setCoordinates(const Real coords[3],Real& i,Real& j) const;

	// Returns the square of the L2-difference between two spherical grids
	static Real SquareDifference(const SphericalGrid& g1,const SphericalGrid& g2);

	// Returns the dot-product of two spherical grids
	static Real Dot(const SphericalGrid& g1,const SphericalGrid& g2);

	// Rotates a spherical grid
	static void Rotate(const SphericalGrid& in,const Real rotation[3][3],SphericalGrid& out);
};
#include "SphericalGrid.inl"
#endif // SPHERICAL_GRID_INCLUDED

