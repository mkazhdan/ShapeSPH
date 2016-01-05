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
#ifndef FOURIER_INCLUDED
#define FOURIER_INCLUDED
#include <fftw3.h>

#include <Util/Algebra.h>

#ifndef PI
#define PI     3.1415926535897932384
#endif
#ifndef SQRT_2
#define SQRT_2 1.4142135623709504880
#endif

#include "RotationGrid.h"
#include "SphericalGrid.h"
#include "SquareGrid.h"
#include "CircularArray.h"
#include "Complex.h"

// This templated class represents the fourier coefficients of a real valued, 1D signal
// Because the input signal generating the key is assumed to be real, we know that the
// negative Fourier coefficients are just conjugates of their positive counterparts, and
// we only store the non-negative coefficients.
template< class Real=float >
class FourierKey1D : public InnerProductSpace< Real , FourierKey1D< Real > >
{
	int res , bw;
	Complex< Real >* values;
public:
    /////////////////////////////////
    // Inner product space methods //
    void Add            ( const FourierKey1D& key );
    void Scale          ( Real s );
    Real InnerProduct   ( const FourierKey1D& key ) const;
    /////////////////////////////////
	
	FourierKey1D( void );
	FourierKey1D( const FourierKey1D& key );
	FourierKey1D( int res );
	~FourierKey1D( void );
	FourierKey1D& operator = ( const FourierKey1D& key );

	// Returns the complex dimension of the array
	int size( void ) const;
	// Returns the resolution of the signal
	int resolution( void ) const;
	// Allocates memory for the array
	int resize( int resolution , bool clr=true );

	// Clears the values of the array to 0
	void clear( void );

	// Returns a reference to the indexed array element
	// 0 <= i <= end()
	Complex< Real >& operator() ( int i );
	// start() <= i <= end()
	Complex< Real >  operator() ( int i ) const;

	// Reads in an array from the specified file
	int read( const char* fileName );
	int read( FILE* fp );

	// Writes out the array to the specified file
	int write( const char* fileName ) const;
	int write( FILE* fp ) const;

	static int Entries( int bw );

	int start( void ) const;
	int end( void ) const;
};

// This templated class represents the fourier coefficients of a real valued, 2D signal
// As in the 1D case, since we assume that the original signal is real, we only store
// half the coefficients.
template< class Real=float >
class FourierKey2D : public InnerProductSpace< Real , FourierKey2D< Real > >
{
	int res , bw;
	Complex<Real>* values;
public:
    /////////////////////////////////
    // Inner product space methods //
    void Add            ( const FourierKey2D& key );
    void Scale          ( Real s );
    Real InnerProduct   ( const FourierKey2D& key ) const;
    /////////////////////////////////

	FourierKey2D( void );
	FourierKey2D( const FourierKey2D& key );
	FourierKey2D( int res );
	~FourierKey2D( void );
	FourierKey2D& operator = ( const FourierKey2D& key );

	// Returns the complex dimension of the array
	int size( void ) const;
	// Returns the resolution of the signal
	int resolution( void ) const;
	// Allocates memory for the array
	int resize( int resolution , bool clr=true );

	// Clears the values of the array to 0
	void clear(void);

	// Returns a reference to the indexed array element.
	// Because we are only storing half the coefficients, we have
	// start() <= i <= end()
	//       0 <= j <= end()
	Complex< Real >& operator() ( int i , int j );
	// start() <= i,j <= end()
	Complex< Real >  operator() ( int i , int j ) const;

	// Reads in an array from the specified file
	int read(const char* fileName);
	int read(FILE* fp);

	// Writes out the array to the specified file
	int write(const char* fileName) const;
	int write(FILE* fp) const;

	static int Entries( int bw );
	int start( void ) const;
	int end( void ) const;
};

// This templated class represents the fourier coefficients of a real valued, signal
// on the surface of the sphers.
// Since we assume that the original signal is real, we only store half the coefficients,
// and for calculations of the dot products we assume that the zonal coefficients are real.
template< class Real=float >
class FourierKeyS2 : public InnerProductSpace< Real , FourierKeyS2< Real > >
{
	int bw;
	Complex< Real >* values;
public:
    /////////////////////////////////
    // Inner product space methods //
    void Add            ( const FourierKeyS2& key );
    void Scale          ( Real s );
    Real InnerProduct   ( const FourierKeyS2& key ) const;
    /////////////////////////////////

	FourierKeyS2( void );
	FourierKeyS2( const FourierKeyS2& key );
	FourierKeyS2( int resolution );
	~FourierKeyS2( void );
	FourierKeyS2& operator = ( const FourierKeyS2& key );

	// Returns the complex dimension of the array
	int bandWidth( void ) const;
	// Returns the resolution of the signal
	int resolution( void ) const;
	// Allocates memory for the array
	int resize( int resolution , bool clr=true );

	// Clears the values of the array to 0
	void clear( void );

	// Returns a reference to the indexed array element
	// In this indexing method, "f" represents the frequency and "i"
	// represents the index of the function within the frequency:
	// 0 <= f < bandWidth()
	// 0 <= i <= f
	Complex<Real>  operator() ( int b , int i ) const;
	Complex<Real>& operator() ( int b , int i );

	// Reads in an array from the specified file
	int read( const char* fileName );
	int read( FILE* fp );

	// Writes out the array to the specified file
	int write( const char* fileName ) const;
	int write( FILE* fp ) const;

	static int Entries( int bw );
};

// This templated class represents the fourier coefficients of a real valued, signal
// on the group of 3D rotations.
// Since we assume that the original signal is real, we only store half the coefficients.
template< class Real=float >
class FourierKeySO3 : public InnerProductSpace< Real , FourierKeySO3< Real > >
{
	int bw;
	Complex<Real>* values;
public:
    /////////////////////////////////
    // Inner product space methods //
    void Add            ( const FourierKeySO3& key );
    void Scale          ( Real s );
    Real InnerProduct   ( const FourierKeySO3& key ) const;
    /////////////////////////////////
	FourierKeySO3( void );
	FourierKeySO3( const FourierKeySO3& key );
	FourierKeySO3( int res );
	~FourierKeySO3( void );
	FourierKeySO3& operator = ( const FourierKeySO3& key );

	// Returns the complex dimension of the array
	int bandWidth( void ) const;
	// Returns the resolution of the signal
	int resolution( void ) const;
	// Allocates memory for the array
	int resize( int resolution , bool clr=true );

	// Clears the values of the array to 0
	void clear( void );

	// Returns a reference to the indexed array element
	// In this indexing method, "f" represents the frequency and the indices "i" and "j"
	// the function within the frequency:
	// 0 <= f < bandWidth()
	//  0 <= i <= f
	// -f <= j <= f
	Complex< Real >  operator() ( int b , int i , int j ) const;
	Complex< Real >& operator() ( int b , int i , int j );

	// Reads in an array from the specified file
	int read( const char* fileName );
	int read( FILE* fp );

	// Writes out the array to the specified file
	int write( const char* fileName ) const;
	int write( FILE* fp ) const;

	static int Entries( int bw );
};

// This templated class is responsible for computing the forward and inverse
// Fourier transforms of periodic functions defined either in 1D or in 2D
template< class Real=float >
class FourierTransform
{
public:
	// This method takes in a real valued function on the circle and computes
	// the Fourier coefficients, writing them into "key"
	int	ForwardFourier( CircularArray<Real>& g , FourierKey1D<Real>& key );

	// This method takes the Fourier coefficients of a real valued function
	// on the circle and returns the originial signal, writing it into "g"
	int InverseFourier( FourierKey1D<Real>& key , CircularArray<Real>& g );

	// This method takes in a real valued function on a 2D grid and computes
	// the Fourier coefficients, writing them into "key"
	int	ForwardFourier( SquareGrid<Real>& g , FourierKey2D<Real>& key );

	// This method takes the Fourier coefficients of a real valued function
	// on a 2D grid and returns the originial signal, writing it into "g"
	// Warning!!! This method will not preserve the input in Fourier key
	int InverseFourier( FourierKey2D<Real>& key , SquareGrid<Real>& g );

	// This static method returns the bandwidth up to which Fourier coefficients
	// will be computed for a signal of resolution "res"
	static int BandWidth( int res ){ return (res>>1)+1; }
};

// This templated class is responsible for computing the forward and inverse
// spherical harmonic transforms of functions defined on the sphere. It allocates
// the appropriate scratch space, based on the resolution of the grid for which
// transforms will be computed.
template<class Real=float>
class HarmonicTransform
{
	class ScratchSpace
	{
	public:
		int bw;
		Real *workSpace , *resultSpace , *transposeResultSpace;
		Real **table , **transposeTable;
		ScratchSpace(void);
		~ScratchSpace(void);
		void resize( const int& bw );
	};
	ScratchSpace scratch;
public:
	HarmonicTransform( void );
	HarmonicTransform( int resolution );
	
	// This method allocates the appropriate amount of scratch space, given
	// the resolution of the signals to be transformed.
	// You do not actually have to call this method, as the transforms will
	// automatically detect if the resolution of the signal doesn't match the
	// resolution of the scratch space, and will call resize if they don't.
	void resize(const int& resolution);

	// This method takes in a real valued function on a sphere and computes
	// the spherical harmonic coefficients, writing them into "key"
	int ForwardFourier(SphericalGrid<Real>& g,FourierKeyS2<Real>& key);

	// This method takes the spherical harmonic coefficients of a real valued function
	// on a sphere and returns the originial signal, writing it into "g"
	int InverseFourier(FourierKeyS2<Real>& key,SphericalGrid<Real>& g);
};

// This templated class is responsible for computing the inverse
// Wigner-D transform of functions defined on the group of 3D rotations. It
// allocates the appropriate scratch space, based on the resolution of the
// grid for which transforms will be computed.
template<class Real=float>
class WignerDTransform{
	class ScratchSpace{
	public:
		int bw;
		fftw_complex *data,*coeffs,*workspace_cx,*workspace_cx2;
		double *workspace_re;
		fftw_plan p;
		ScratchSpace(void);
		~ScratchSpace(void);

		void resize(const int& bw);
	};
	ScratchSpace scratch;
public:
	// This method allocates the appropriate amount of scratch space, given
	// the resolution of the signals to be transformed.
	// You do not actually have to call this method, as the transforms will
	// automatically detect if the resolution of the signal doesn't match the
	// resolution of the scratch space, and will call resize if they don't.
	void resize(const int& resolution);

	// This method takes the spherical harmonic coefficients of a real valued function
	// on a sphere and returns the originial signal, writing it into "g"
	int InverseFourier(FourierKeySO3<Real>& key,RotationGrid<Real>& g);
};

#include "Fourier1D.inl"
#include "Fourier2D.inl"
#include "FourierS2.inl"
#include "FourierSO3.inl"
#endif // FOURIER_INCLUDED
