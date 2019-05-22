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
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include "SignalProcessing/CubeGrid.h"
#include "SignalProcessing/Fourier.h"
#include "Util/Ply.h"
#include "Util/CmdLineParser.h"
#include "Util/Geometry.h"
#include "Util/Util.h"
#include "Util/Rasterizer.h"
#include "Util/EDT.h"
#include "Util/SphereSampler.h"
#include "Util/lineqn.h"
#include "Util/TriangleMesh.h"
#include "Util/SphericalPolynomials.h"
#include "Util/Signature.h"

cmdLineString In( "in" ) , Out( "out" );
cmdLineInt Resolution( "res" , 64 ) , BandWidth( "bw" , 16 ) , Radii( "radii" , 32 ) , AnisotropicScale( "aScale" , 0 ) , Threads( "threads" , omp_get_num_procs() );
cmdLineFloat MomentRadiusScale( "radius" , 2.f ) , FallOff( "fallOff" , float(sqrt(8.)) );
cmdLineReadable NoCQ( "noCQ" ) , Double( "double" ) , Verbose( "verbose" ) , Binary( "binary" );

cmdLineReadable* params[] = { &In , &Out , &Resolution , &BandWidth , &Radii , &AnisotropicScale , &Threads , &MomentRadiusScale , &FallOff , &NoCQ , &Double , &Verbose , &Binary , NULL };

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input mesh>\n" , In.name );
	printf( "\t[--%s <output signature>]\n" , Out.name );
	printf( "\t[--%s <voxel resolution>=%d]\n" , Resolution.name , Resolution.value );
	printf( "\t[--%s <sph band-width>=%d]\n" , BandWidth.name , BandWidth.value );
	printf( "\t[--%s <sampling radii>=%d]\n" , Radii.name , Radii.value );
	printf( "\t[--%s <threads>=%d]\n" , Threads.name , Threads.value );
	printf( "\t[--%s <anisotropic scale>=%d]\n" , AnisotropicScale.name , AnisotropicScale.value );
	printf( "\t[--%s <moment radius scale>=%f]\n" , MomentRadiusScale.name , MomentRadiusScale.value );
	printf( "\t[--%s <Gaussian EDT fall off>=%f]\n" , FallOff.name , FallOff.value );
	printf( "\t[--%s]\n" , NoCQ.name );
	printf( "\t[--%s]\n" , Double.name );
	printf( "\t[--%s]\n" , Binary.name );
	printf( "\t[--%s]\n" , Verbose.name );
}

template< class Real >
int run( void )
{
	std::vector< Point3D< Real > > vertices;
	std::vector< TriangleIndex > triangles;
	// Read in the mesh
	{
		int fileType;
		std::vector< PlyColorVertex< float > > _vertices;
		PlyReadTriangles( In.value , _vertices , triangles , PlyColorVertex< float >::ReadProperties , NULL , PlyColorVertex< float >::ReadComponents , fileType );
		vertices.resize( _vertices.size() );
		for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] = Point3D< Real >( _vertices[i].point );
		SquareMatrix< Real , 4 > xForm = GetAligningXForm( vertices , triangles , Real(MomentRadiusScale.value) , AnisotropicScale.value );
		XForm( vertices , xForm , Threads.value );
	}

	// Normalize translation and scale
#pragma omp parallel for num_threads( Threads.value )
	for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] = vertices[i] * Real(0.5) + Point3D< Real >( Real(0.5) , Real(0.5) , Real(0.5) );

	CubeGrid< char > grid;
	CubeGrid< Real > gedt;
	std::vector< FourierKeyS2< Real > > sKeys;
	double t;

	// Compute the rasterization
	t = Time();
	{
		grid.resize( Resolution.value );
		Rasterize( vertices , triangles , grid , Real(1.5) , Threads.value );
	}
	if( Verbose.set ) printf( "\tRasterization time: %.2f(s)\n" , Time()-t );

	// Compute the Gaussian EDT
	t = Time();
	{
		GaussianEDT( grid , gedt , Real(FallOff.value) , Threads.value );
	}
	if( Verbose.set ) printf( "\tGaussian EDT time: %.2f(s)\n" , Time()-t );

	// Compute the spherical harmonic decomposition of the concentric spheres (weighted by radius)
	t = Time();
	{
		Real radius = Real(Resolution.value)/2;
		Point3D< Real > center = Point3D< Real >( radius , radius , radius );
		SampleSpheres( gedt , sKeys , center , radius , Radii.value , Resolution.value , Threads.value );
		
		Real norm = 0;
		for( int i=0 ; i<sKeys.size() ; i++ ) norm += sKeys[i].squareNorm();
		norm = Real( sqrt(norm) );
		for( int i=0 ; i<sKeys.size() ; i++ ) sKeys[i] /= norm;
	}
	if( Verbose.set ) printf( "\tSpherical Harmonic time: %.2f(s)\n" , Time()-t );

	// Compute the signature
	if( Out.set )
	{
		int bw = BandWidth.value;
		Signature< Real > sig( !NoCQ.set ? (bw+1) * int( sKeys.size() ) : bw * int( sKeys.size() ) );
		for( int i=0 ; i<sKeys.size() ; i++ )
			if( !NoCQ.set )
			{
				int idx = 0;
				Point3D< Real > cq = ConstantAndQuadratic( sKeys[i] );
				for( int j=0 ; j<3 ; j++ ) sig[i*(bw+1)+idx ] = cq[j] , idx++;
				for( int b=0 ; b<bw ; b++ ) if( b!=0 && b!=2 )
				{
					Real _norm2 = sKeys[i](b,0).squareNorm();
					for( int j=1 ; j<=b ; j++ ) _norm2 += sKeys[i](b,j).squareNorm()*2;
					sig[i*(bw+1)+idx] = Real( sqrt(_norm2) );
					idx++;
				}
			}
			else
			{
				for( int b=0 ; b<bw ; b++ )
				{
					Real _norm2 = sKeys[i](b,0).squareNorm();
					for( int j=1 ; j<=b ; j++ ) _norm2 += sKeys[i](b,j).squareNorm()*2;
					sig[i*bw+b] = Real( sqrt(_norm2) );
				}
			}
		sig.write( Out.value , Binary.set );
	}
	return EXIT_SUCCESS;
}
int main( int argc , char* argv[] )
{
	std::vector<std::string> nonoptArgs;
	cmdLineParse( argc , argv , params , nonoptArgs );
	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}
	if( BandWidth.value>Resolution.value/2 )
	{
		fprintf( stderr , "[WARNING] Resetting band-width: %d -> %d\n" , BandWidth.value , Resolution.value/2 );
		BandWidth.value = Resolution.value/2;
	}
	int ret;
	double t = Time();
	if( Double.set ) ret = run< double >();
	else             ret = run< float  >();
	if( Verbose.set ) printf( "Total Time: %.2f(s)\n" , Time()-t );
	return ret;
}