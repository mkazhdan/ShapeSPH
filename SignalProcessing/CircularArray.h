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
#ifndef CIRCULAR_ARRAY_INCLUDED
#define CIRCULAR_ARRAY_INCLUDED

#ifndef PI
#define PI 3.1415926535897932384
#endif 

// This templated class represents a set of regular samples of a real-valued function on the unit-circle
template<class Real=float>
class CircularArray{
protected:
	int res;
	Real* values;
public:
	CircularArray( void );
	CircularArray( int res );
	~CircularArray( void );

	// Returns the dimension of the array
	int resolution( void ) const;
	// Allocates memory for the array
	int resize( int resolution );

	// Clears the values of the array to 0
	void clear(void);

	// Returns a reference to the indexed array element
	Real& operator() (const int& x);
	Real operator() (const int& x) const;
	// Returns the linear interpolation of the value at the spedified index
	Real operator() (const double& x);
	Real operator() (const double& x) const;

	// Returns the square of the L2-norm of the array elements
	Real squareNorm(void) const;

	// Reads in an array from the specified file
	int read(const char* fileName);
	int read(FILE* fp);

	// Writes out the array to the specified file
	int write(const char* fileName) const;
	int write(FILE* fp) const;

	// Returns the square of the L2-difference between two arrays
	static Real SquareDifference(const CircularArray& g1,const CircularArray& g2);

	// Returns the dot-product of two arrays
	static Real Dot(const CircularArray& g1,const CircularArray& g2);
};
#include "CircularArray.inl"
#endif // CIRCULAR_ARRAY_INCLUDED
