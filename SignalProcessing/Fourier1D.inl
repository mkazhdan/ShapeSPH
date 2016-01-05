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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

//#include "FFTW/fftw3.h"
#include "fftw3.h"
#include <math.h>

//////////////////
// FourierKey1D //
//////////////////
template< class Real > int FourierKey1D< Real >::Entries( int bw ) { return bw; }
template<class Real> FourierKey1D<Real>::FourierKey1D( void ) : bw(0) , res(0) , values(NULL) { ; }
template<class Real> FourierKey1D<Real>::FourierKey1D( int r ) : bw(0) , res(0) , values(NULL) { resize( r ); }
template<class Real> FourierKey1D<Real>::FourierKey1D( const FourierKey1D< Real >& key ) : bw(0) , res(0) , values(NULL)
{
	resize( key.res );
	memcpy( values , key.values , sizeof( Complex< Real > ) * Entries(bw) );
}
template< class Real >
FourierKey1D< Real >& FourierKey1D< Real >::operator = ( const FourierKey1D< Real >& key )
{
	resize( key.res );
	memcpy( values , key.values , sizeof( Complex< Real > ) * Entries(bw) );
	return *this;
}
template<class Real> FourierKey1D<Real>::~FourierKey1D( void )
{
	if( values ) delete[] values;
	values = NULL;
	bw = res = 0;
}
template<class Real> int FourierKey1D<Real>::read( const char* fileName )
{
	FILE* fp = fopen( fileName , "rb" );
	if( !fp ) return 0;
	int r = read(fp);
	fclose(fp);
	return r;
}
template<class Real> int FourierKey1D<Real>::write(const char* fileName) const{
	FILE* fp=fopen(fileName,"wb");
	if(!fp){return 0;}
	int w=write(fp);
	fclose(fp);
	return w;
}
template<class Real> int FourierKey1D<Real>::read(FILE* fp){
	int resolution,r;
	r = int( fread( &resolution , sizeof(int) , 1 , fp ) );
	if( !r ) return 0; 
	resize( resolution );
	r = int( fread( values , sizeof( Complex<Real> ) , bw , fp ) );
	if( r==bw ) return 1;
	else        return 0;
}
template<class Real> int FourierKey1D<Real>::write( FILE* fp ) const
{
	int w;
	w = int( fwrite( &res , sizeof(int) , 1 , fp ) );
	if( !w ) return 0;
	w = int( fwrite( values , sizeof( Complex<Real> ) , bw , fp ) );
	if( w==bw ) return 1;
	else        return 0;
}
template< class Real > int FourierKey1D< Real >::size( void ) const { return bw; }
template< class Real > int FourierKey1D< Real >::resolution( void ) const { return res; }
template< class Real > int FourierKey1D< Real >::resize( int resolution , bool clr )
{
	int b = FourierTransform<Real>::BandWidth( resolution );
	if( resolution<0 ) return 0;
	else if( b!=bw )
	{
		if(values) delete[] values;
		values = NULL;
		bw = res = 0;
		if(b)
		{
			values = new Complex<Real>[b];
			if( !values ) return 0;
			else bw=b;
		}
	}
	res = resolution;
	if( clr ) clear();
	return 1;
}
template< class Real > void FourierKey1D< Real >::clear( void ){ if(bw) memset( values , 0 , sizeof( Complex<Real> )*bw); }
template< class Real >
void FourierKey1D< Real >::Add( const FourierKey1D< Real >& key )
{
	for( int i=0 ; i<=std::min< int >( end() , key.end() ) ; i++ ) values[i] += key.values[i];
}
template< class Real >
void FourierKey1D< Real >::Scale( Real s )
{
	for( int i=0 ; i<=end() ; i++ ) values[i] *= s;
}
template< class Real >
Real FourierKey1D< Real >::InnerProduct( const FourierKey1D< Real >& key ) const
{
	Real dot = 0;
	for( int i=std::max< int >( start() , key.start() ) ; i<=std::min< int >( end() , key.end() ) ; i++ ) dot += ( (*this)(i)*key(i).conjugate() ).r;
	return dot;
}
template< class Real > Complex< Real >& FourierKey1D<Real>::operator() ( int i ) { return values[i]; }
template< class Real > Complex< Real >  FourierKey1D<Real>::operator() ( int i ) const
{
	if( i<start() || i>end() ) return Complex< Real >( 0 );
	else if( i==end() && !(res&1) ) return Complex< Real >( values[i].r );
	else return i>=0 ? values[i] : values[-i].conjugate();
}
template< class Real > int FourierKey1D< Real >::start( void ) const { return res&1 ? -bw+1 : -bw+2; }
template< class Real > int FourierKey1D< Real >::end  ( void ) const { return bw-1; }

//////////////////////
// FourierTransform //
//////////////////////
template< class Real > void __FFTW_1D__( int r , void* iValues , void* oValues );
template<> void __FFTW_1D__< float >( int r , void* iValues , void* oValues )
{
	fftwf_plan plan = fftwf_plan_dft_r2c_1d( r , (float*)iValues , (fftwf_complex*)oValues , FFTW_PRESERVE_INPUT | FFTW_ESTIMATE );
	fftwf_execute( plan );
	fftwf_destroy_plan( plan );
}
template<> void __FFTW_1D__< double >( int r , void* iValues , void* oValues )
{
	fftw_plan plan = fftw_plan_dft_r2c_1d( r , (double*)iValues , (fftw_complex*)oValues , FFTW_PRESERVE_INPUT | FFTW_ESTIMATE );
	fftw_execute( plan );
	fftw_destroy_plan( plan );
}
template< class Real > void __FFTW_1D__< Real >( int r , void* iValues , void* oValues )
{
	fprintf( stderr , "[ERROR] FFTW_1D only supported for floats and doubles\n" );
	exit( 0 );
}
template< class Real >
int FourierTransform< Real >::ForwardFourier( CircularArray< Real >& g , FourierKey1D< Real >& key )
{
	if( key.resolution()!=g.resolution() ) key.resize( g.resolution() );
	__FFTW_1D__< Real >( g.resolution() , &g(0) , &key(0) );
	// The forward transform takes the constant function f(theta) = 1 => key(0) = res
	// We want it to take the constant function: f(theta) = 1/sqrt(2*Pi) => key(0) = 1 <=> f(theta) = 1 => key(0) = sqrt(2*Pi)
	Real n = Real( sqrt( 2.0*PI ) ) / Real( g.resolution() );
	for( int i=0 ; i<=key.end() ; i++ ) key(i) *= n;
	return 1;
}
template< class Real > void __IFFTW_1D__( int r , void* iValues , void* oValues );
template<> void __IFFTW_1D__< float >( int r , void* iValues , void* oValues )
{
	fftwf_plan plan = fftwf_plan_dft_c2r_1d( r , (fftwf_complex*)iValues , (float*)oValues , FFTW_PRESERVE_INPUT | FFTW_ESTIMATE );
	fftwf_execute( plan );
	fftwf_destroy_plan( plan );
}
template<> void __IFFTW_1D__< double >( int r , void* iValues , void* oValues )
{
	fftw_plan plan = fftw_plan_dft_c2r_1d( r , (fftw_complex*)iValues , (double*)oValues , FFTW_PRESERVE_INPUT | FFTW_ESTIMATE );
	fftw_execute( plan );
	fftw_destroy_plan( plan );
}
template< class Real > void __IFFTW_1D__< Real >( int r , void* iValues , void* oValues )
{
	fprintf( stderr , "[ERROR] IFFTW_1D only supported for floats and doubles\n" );
	exit( 0 );
}
template< class Real >
int FourierTransform< Real >::InverseFourier( FourierKey1D< Real >& key , CircularArray< Real >& g )
{
	if( key.resolution()!=g.resolution() ) g.resize( key.resolution() );
	__IFFTW_1D__< Real >( g.resolution() , &key(0) , &g(0) );
	// The inverse transform takes the coefficient: key(0) = 1  => f(theta) = 1
	// We want it to take the coefficient: key(0) = 1 => f(theta) = 1/sqrt(2*Pi)
	Real n = Real( sqrt( 1./(2.0 * PI) ) );
	for( int i=0 ; i<g.resolution() ; i++ ) g(i) *= n;
	return 1;
}
