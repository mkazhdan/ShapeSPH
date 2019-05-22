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
#ifndef ROTATION_GRID_INCLUDED
#define ROTATION_GRID_INCLUDED


#include "CubeGrid.h"
#ifndef PI
#define PI 3.1415926535897932384
#endif 


// This templated class represents a function on the group of rotations, where index (i,j,k) corresponds
// to the rotation with Euler coefficients:
//			theta =2PI*j/r
//			phi   =PI*(2i+1)/(2r)
//			psi   =2PI*k/r
// where r is the resolution of the sampling.

template< class Real=float >
class RotationGrid : public CubeGrid< Real >{
public:
	using CubeGrid< Real >::res;
	using CubeGrid< Real >::values;
	using CubeGrid< Real >::resize;
	RotationGrid( void );
	RotationGrid( int res );
	~RotationGrid( void );

	// Returns a reference to the indexed  array element
	Real& operator() (const int& x,const int& y,const int& z);
	Real operator() (const int& x,const int& y,const int& z) const;
	// Returns the linear interpolation of the value at the spedified index
	Real sample( Real x , Real y , Real z ) const;

	// Sets the matrix coordinates of the Euler angle (theta,phi,psi)
	static void SetCoordinates(const Real& theta,const Real& phi,const Real& psi,Real matrix[3][3]);
	// Sets the coordinates (theta,phi,psi) associated to the specified rotation matrix
	static void SetCoordinates(const Real matrix[3][3],Real& theta,Real& phi,Real& psi);
	// Sets the matrix coordinates of the the rotation about the specified axis with the
	// specified angle
	static void SetCoordinates(const Real axis[3],const Real& angle,Real matrix[3][3]);

	// Sets the matrix coordinates of the Euler angle indexed by (i,j,k)
	void setCoordinates(const int& i,const int& j,const int& k,Real matrix[3][3]) const;
	void setCoordinates(const Real& i,const Real& j,const Real& k,Real matrix[3][3]) const;
	// Sets the coordinates (i,j,k) associated to the specified rotation matrix
	void setCoordinates(const Real matrix[3][3],Real& i,Real& j,Real& k) const;

	// Returns the square of the L2-norm of the array elements
	Real squareNorm(void) const;

	// Returns the square of the L2-difference between two arrays
	static Real SquareDifference(const RotationGrid& g1,const RotationGrid& g2);

	// Returns the dot-product of two arrays
	static Real Dot(const RotationGrid& g1,const RotationGrid& g2);
};
#include "RotationGrid.inl"
#endif // ROTATION_GRID_INCLUDED

