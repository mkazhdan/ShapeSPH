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

#include "Util/Ply.h"
#include "Util/CmdLineParser.h"
#include "Util/Geometry.h"
#include "Util/Util.h"
#include "SignalProcessing/CubeGrid.h"
#include "SignalProcessing/Fourier.h"
#include "Util/Rasterizer.h"
#include "Util/EDT.h"
#include "Util/TriangleMesh.h"

cmdLineString In1( "in1" ) , In2( "in2" ) , Out( "out" );
cmdLineInt Resolution( "res" , 64 ) , AnisotropicScale( "aScale" , 0 ) , Threads( "threads" , omp_get_num_procs() );
cmdLineFloat MomentRadiusScale( "radius" , 2.f ) , FallOff( "fallOff" , float( sqrt(8.) ) );
cmdLineReadable GEDT( "gedt" ) , Double( "double" ) , Verbose( "verbose" );

cmdLineReadable* params[] = { &In1 , &In2 , &Out , &AnisotropicScale , &Resolution , &Threads , &MomentRadiusScale , &GEDT , &FallOff , &Double , &Verbose , NULL };

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <source mesh>\n" , In1.name );
	printf( "\t --%s <target mesh>\n" , In2.name );
	printf( "\t[--%s <aligned source mesh>]\n" , Out.name );
	printf( "\t[--%s <voxel resolution>=%d]\n" , Resolution.name , Resolution.value );
	printf( "\t[--%s <threads>=%d]\n" , Threads.name , Threads.value );
	printf( "\t[--%s <anisotropic scale>=%d]\n" , AnisotropicScale.name , AnisotropicScale.value );
	printf( "\t[--%s <moment radius scale>=%f]\n" , MomentRadiusScale.name , MomentRadiusScale.value );
	printf( "\t[--%s <Gaussian fall-off>=%f]\n" , FallOff.name , FallOff.value );
	printf( "\t[--%s]\n" , GEDT.name );
	printf( "\t[--%s]\n" , Double.name );
	printf( "\t[--%s]\n" , Verbose.name );
}

template< class Real >
void SetSphereKeys( const CubeGrid< Real >& grid , HarmonicTransform< Real >& xForm , std::vector< FourierKeyS2< Real > >& keys , int threads )
{
	int res = grid.resolution();
	xForm.resize( res );
	keys.resize( res/2 );
	for( int i=0 ; i<res/2 ; i++ ) keys[i].resize( res );
	SphericalGrid< Real > sGrid;
	sGrid.resize( res );

	Real radius = Real(res)/2;
	Point3D< Real > center( radius , radius , radius );
	for( int i=1 ; i<=res/2 ; i++ )
	{
		Real r = radius * Real(i)/(res/2);
		grid.SphereSample( &center[0] , r , sGrid , threads );
		Real scale = Real( sqrt( 4*M_PI*r*r ) );
		Real* _sGrid = sGrid[0];
#pragma omp parallel for num_threads( threads )
		for( int j=0 ; j<sGrid.resolution()*sGrid.resolution() ; j++ ) _sGrid[j] *= scale;
		xForm.ForwardFourier( sGrid , keys[i-1] );
	}
}
template< class Real >
void AddWignerDCoefficients( const std::vector< FourierKeyS2< Real > >& inKeys1 , const std::vector< FourierKeyS2< Real > >& inKeys2 , FourierKeySO3< Real >& outKey , int threads=1 )
{
	for( int b=0 ; b<outKey.bandWidth() ; b++ ) if( inKeys1[0].bandWidth()>b && inKeys2[0].bandWidth()>b )
	{
#pragma omp parallel for num_threads( threads )
		for( int ij=0 ; ij<(b+1)*(b+1) ; ij++ )
		{
			int i = ij / (b+1) , j = ij % (b+1);
			Complex< Real > temp;
			for( int k=0 ; k<inKeys1.size() && k<inKeys2.size() ; k++ ) temp += inKeys1[k](b,i).conjugate() * inKeys2[k](b,j);
			outKey(b,i,j) += temp;

			temp *= 0;
			if( i && j )
			{
				for( int k=0 ; k<inKeys1.size() && k<inKeys2.size() ; k++ ) temp += inKeys1[k](b,i).conjugate() * inKeys2[k](b,j).conjugate();
				outKey(b,i,-j) += temp;
			}
		}
	}
}

template< class Real >
SquareMatrix< Real , 3 > _main_( const std::vector< Point3D< Real > >& vertices1 , const std::vector< TriangleIndex >& triangles1 , const std::vector< Point3D< Real > >& vertices2 , const std::vector< TriangleIndex >& triangles2 )
{
	CubeGrid< char > grid1 , grid2;
	CubeGrid< Real > raster1 , raster2 , sqr_edt1 , sqr_edt2 , gedt1 , gedt2;

	double t;
	t = Time();
	{
		grid1.resize( Resolution.value) , grid2.resize( Resolution.value );
		Rasterize( vertices1 , triangles1 , grid1 , Real(1.5) , Threads.value );
		Rasterize( vertices2 , triangles2 , grid2 , Real(1.5) , Threads.value );

		if( !GEDT.set )
		{
			raster1.resize( Resolution.value ) , raster2.resize( Resolution.value );
			char *_grid1 = grid1[0] , *_grid2 = grid2[0];
			Real *_raster1 = raster1[0] , *_raster2 = raster2[0];
#pragma omp parallel for num_threads( Threads.value )
			for( int i=0 ; i<Resolution.value*Resolution.value*Resolution.value ; i++ ) _raster1[i] = Real( _grid1[i] ) , _raster2[i] = Real( _grid2[i] );
		}
	}
	if( Verbose.set ) printf( "\t\tRasterization Time: %.2f(s)\n" , Time()-t );

	t = Time();
	if( GEDT.set )
	{
		GaussianEDT( grid1 , gedt1 , Real( FallOff.value ) , Threads.value );
		GaussianEDT( grid2 , gedt2 , Real( FallOff.value ) , Threads.value );
	}
	else
	{
		SquaredEDT( grid1 , sqr_edt1 , Threads.value );
		SquaredEDT( grid2 , sqr_edt2 , Threads.value );
	}
	if( Verbose.set ) printf( "\t\tEDT Time: %.2f(s)\n" , Time()-t );


	std::vector< FourierKeyS2< Real > > rasterKey1 , rasterKey2 , edtKey1 , edtKey2 , gedtKey1 , gedtKey2;
	t = Time();
	{
		HarmonicTransform< Real > xForm;
		xForm.resize( Resolution.value );
		if( GEDT.set )
		{
			SetSphereKeys( gedt1 , xForm , gedtKey1 , Threads.value );
			SetSphereKeys( gedt2 , xForm , gedtKey2 , Threads.value );
		}
		else
		{
			SetSphereKeys( raster1 , xForm , rasterKey1 , Threads.value );
			SetSphereKeys( raster2 , xForm , rasterKey2 , Threads.value );
			SetSphereKeys( sqr_edt1 , xForm , edtKey1 , Threads.value );
			SetSphereKeys( sqr_edt2 , xForm , edtKey2 , Threads.value );
		}
	}
	if( Verbose.set ) printf( "\t\tHarmonic Key Time: %.2f(s)\n" , Time()-t );


	RotationGrid< Real > so3Grid;
	Real norm2 = 0;
	t = Time();
	{
		WignerDTransform< Real > xForm;
		FourierKeySO3< Real > key;
		xForm.resize( Resolution.value );
		key.resize( Resolution.value );
		if( GEDT.set ) AddWignerDCoefficients( gedtKey1 , gedtKey2 , key , Threads.value );
		else           AddWignerDCoefficients( rasterKey1 , edtKey2 , key , Threads.value ) , AddWignerDCoefficients( edtKey1 , rasterKey2 , key , Threads.value );
		if( GEDT.set ) for( int r=0 ; r<Resolution.value/2 ; r++ ) norm2 += gedtKey1[r].squareNorm() + gedtKey2[r].squareNorm();
		xForm.InverseFourier( key , so3Grid );
	}
	if( Verbose.set ) printf( "\t\tWigner-D Time: %.2f(s)\n" , Time()-t );

	SquareMatrix< Real , 3 > minRotation=SquareMatrix< Real , 3 >::Identity() , maxRotation=SquareMatrix< Real , 3 >::Identity();
	Real minCorrelation=0 , maxCorrelation=0;
	for( int i=0 ; i<Resolution.value ; i++ ) for( int j=0 ; j<Resolution.value ; j++ ) for( int k=0 ; k<Resolution.value ; k++ )
	{
		if( (!i && !j && !k) || so3Grid(i,j,k)<minCorrelation )
		{
			Real matrix[3][3];
			so3Grid.setCoordinates( i , j , k , matrix );
			for( int ii=0 ; ii<3 ; ii++ ) for( int jj=0 ; jj<3 ; jj++ ) minRotation(ii,jj) = matrix[jj][ii];
			minCorrelation = so3Grid(i,j,k);
		}
		if( (!i && !j && !k) || so3Grid(i,j,k)>maxCorrelation )
		{
			Real matrix[3][3];
			so3Grid.setCoordinates( i , j , k , matrix );
			for( int ii=0 ; ii<3 ; ii++ ) for( int jj=0 ; jj<3 ; jj++ ) maxRotation(ii,jj) = matrix[jj][ii];
			maxCorrelation = so3Grid(i,j,k);
		}
	}
	if( GEDT.set )
	{
		if( Verbose.set ) printf( "\tMin/Max Error: %e / %e\n" , sqrt( norm2 - maxCorrelation*2 ) , sqrt( norm2 - minCorrelation*2 ) );
		return maxRotation;
	}
	else
	{
		if( Verbose.set ) printf( "\tMin/Max Correlation: %e / %e\n" , minCorrelation , maxCorrelation );
		return minRotation;
	}
}
template< class Real >
void WriteTriangles( const std::vector< Point3D< Real > >& vertices , const std::vector< TriangleIndex >& triangles , char* fileName )
{
	std::vector< PlyVertex< float > > _vertices( vertices.size() );
	for( int i=0 ; i<vertices.size() ; i++ ) _vertices[i].point = Point3D< float >( vertices[i] );
	PlyWriteTriangles( fileName , _vertices , triangles , PlyVertex< float >::WriteProperties , PlyVertex< float >::WriteComponents , PLY_BINARY_NATIVE );
}

template< class Real >
int _main_( void )
{
	SquareMatrix< Real , 4 > xForm1 , xForm2;
	std::vector< Point3D< Real > > vertices1 , vertices2;
	std::vector< TriangleIndex > triangles1 , triangles2;
	{
		int fileType;
		std::vector< PlyColorVertex< float > > _vertices;
		PlyReadTriangles( In1.value , _vertices , triangles1 , PlyColorVertex< float >::ReadProperties , NULL , PlyColorVertex< float >::ReadComponents , fileType );
		vertices1.resize( _vertices.size() );
		for( int i=0 ; i<vertices1.size() ; i++ ) vertices1[i] = Point3D< Real >( _vertices[i].point );
		xForm1 = GetAligningXForm( vertices1 , triangles1 , Real(MomentRadiusScale.value) , AnisotropicScale.value );
		XForm( vertices1 , xForm1 , Threads.value );
	}
	{
		int fileType;
		std::vector< PlyColorVertex< float > > _vertices;
		PlyReadTriangles( In2.value , _vertices , triangles2 , PlyColorVertex< float >::ReadProperties , NULL , PlyColorVertex< float >::ReadComponents , fileType );
		vertices2.resize( _vertices.size() );
		for( int i=0 ; i<vertices2.size() ; i++ ) vertices2[i] = Point3D< Real >( _vertices[i].point );
		xForm2 = GetAligningXForm( vertices2 , triangles2 , Real(MomentRadiusScale.value) , AnisotropicScale.value );
		XForm( vertices2 , xForm2 , Threads.value );
	}
	SquareMatrix< Real , 3 > rotation;
	double t = Time();
	{
#pragma omp parallel for num_threads( Threads.value )
		for( int i=0 ; i<vertices1.size() ; i++ ) vertices1[i] = vertices1[i] * Real(0.5) + Point3D< Real >( Real(0.5) , Real(0.5) , Real(0.5) );
#pragma omp parallel for num_threads( Threads.value )
		for( int i=0 ; i<vertices2.size() ; i++ ) vertices2[i] = vertices2[i] * Real(0.5) + Point3D< Real >( Real(0.5) , Real(0.5) , Real(0.5) );
		rotation = _main_( vertices1 , triangles1 , vertices2 , triangles2 );
#pragma omp parallel for num_threads( Threads.value )
		for( int i=0 ; i<vertices1.size() ; i++ ) vertices1[i] = ( vertices1[i] - Point3D< Real >( Real(0.5) , Real(0.5) , Real(0.5) ) ) * Real(2.);
#pragma omp parallel for num_threads( Threads.value )
		for( int i=0 ; i<vertices2.size() ; i++ ) vertices2[i] = ( vertices2[i] - Point3D< Real >( Real(0.5) , Real(0.5) , Real(0.5) ) ) * Real(2.);
	}
	if( Verbose.set ) printf( "\tAlignment Time: %.2f(s)\n" , Time()-t );
	printf( "\t%f %f %f\n\t%f %f %f\n\t%f %f %f\n" , rotation(0,0) , rotation(1,0) , rotation(2,0) , rotation(0,1) , rotation(1,1) , rotation(2,1) , rotation(0,2) , rotation(1,2) , rotation(2,2) );

	if( Out.set )
	{
		SquareMatrix< Real , 4 > _rotation = SquareMatrix< Real , 4 >::Identity();
		for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) _rotation(i,j) = rotation(i,j);
		SquareMatrix< Real , 4 > xForm = xForm2.inverse() * _rotation;
		XForm( vertices1 , xForm );
		WriteTriangles( vertices1 , triangles1 , Out.value );
	}
	return EXIT_SUCCESS;
}
int main( int argc , char* argv[] )
{
	std::vector<std::string> nonoptArgs;
	cmdLineParse( argc , argv , params , nonoptArgs );
	if( !In1.set || !In2.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}

	int ret;
	double t = Time();
	if( Double.set ) ret = _main_< double >( );
	else             ret = _main_< float  >( );
	if( Verbose.set ) printf( "Running Time: %.2f(s)\n" , Time()-t );

	return ret;
}