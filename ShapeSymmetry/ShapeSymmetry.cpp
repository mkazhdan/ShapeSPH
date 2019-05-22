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

cmdLineString In( "in" ) , OutHeader( "out" );
cmdLineInt Resolution( "res" , 64 ) , Threads( "threads" , omp_get_num_procs() ) , MaxRotationalSymmetry( "maxSym" , 6 );
cmdLineFloat MomentRadiusScale( "radius" , 2.f ) , FallOff( "fallOff" , float( sqrt(8.) ) );
cmdLineReadable GEDT( "gedt" ) , Double( "double" ) , Verbose( "verbose" );

cmdLineReadable* params[] = { &In , &OutHeader , &Resolution , &Threads , &MomentRadiusScale , &GEDT , &FallOff , &Double , &MaxRotationalSymmetry , &Verbose , NULL };

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <source mesh>\n" , In.name );
	printf( "\t[--%s <output header>]\n" , OutHeader.name );
	printf( "\t[--%s <voxel resolution>=%d]\n" , Resolution.name , Resolution.value );
	printf( "\t[--%s <threads>=%d]\n" , Threads.name , Threads.value );
	printf( "\t[--%s <moment radius scale>=%f]\n" , MomentRadiusScale.name , MomentRadiusScale.value );
	printf( "\t[--%s <Gaussian fall-off>=%f]\n" , FallOff.name , FallOff.value );
	printf( "\t[--%s <maximal order of rotational symmetry>=%d]\n" , MaxRotationalSymmetry.name , MaxRotationalSymmetry.value );
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
void AddWignerDCoefficients( const std::vector< FourierKeyS2< Real > >& inKeys1 , const std::vector< FourierKeyS2< Real > >& inKeys2 , FourierKeySO3< Real >& outKey1 , FourierKeySO3< Real >& outKey2 , int threads=1 )
{
	for( int b=0 ; b<outKey1.bandWidth() && b<outKey2.bandWidth() ; b++ ) if( inKeys1[0].bandWidth()>b && inKeys2[0].bandWidth()>b )
	{
		Real sign = Real( (b%2) ? -1 : 1 );
#pragma omp parallel for num_threads( threads )
		for( int ij=0 ; ij<(b+1)*(b+1) ; ij++ )
		{
			int i = ij / (b+1) , j = ij % (b+1);
			Complex< Real > temp;
			for( int k=0 ; k<inKeys1.size() && k<inKeys2.size() ; k++ ) temp += inKeys1[k](b,i).conjugate() * inKeys2[k](b,j);
			outKey1(b,i,j) += temp , outKey2(b,i,j) += temp * sign;

			temp *= 0;
			if( i && j )
			{
				for( int k=0 ; k<inKeys1.size() && k<inKeys2.size() ; k++ ) temp += inKeys1[k](b,i).conjugate() * inKeys2[k](b,j).conjugate();
				outKey1(b,i,-j) += temp , outKey2(b,i,-j) += temp * sign;
			}
		}
	}
}
template< class Real >
void _main_( const std::vector< Point3D< Real > >& vertices , const std::vector< TriangleIndex >& triangles , SphericalGrid< Real >& axialSymmetry , SphericalGrid< Real >& refSymmetry , std::vector< SphericalGrid< Real > >& rotSymmetry )
{
	CubeGrid< char > grid;
	CubeGrid< Real > raster , sqr_edt , gedt;

	double t;
	t = Time();
	{
		grid.resize( Resolution.value );
		Rasterize( vertices , triangles , grid , Real(1.5) , Threads.value );

		if( !GEDT.set )
		{
			raster.resize( Resolution.value );
			char *_grid = grid[0];
			Real *_raster = raster[0];
#pragma omp parallel for num_threads( Threads.value )
			for( int i=0 ; i<Resolution.value*Resolution.value*Resolution.value ; i++ ) _raster[i] = Real( _grid[i] );
		}
	}
	if( Verbose.set ) printf( "\t\tRasterization Time: %.2f(s)\n" , Time()-t );

	t = Time();
	if( GEDT.set ) GaussianEDT( grid , gedt , Real( FallOff.value ) , Threads.value );
	else            SquaredEDT( grid , sqr_edt , Threads.value );
	if( Verbose.set ) printf( "\t\tEDT Time: %.2f(s)\n" , Time()-t );

	std::vector< FourierKeyS2< Real > > rasterKey , edtKey , gedtKey;
	t = Time();
	{
		HarmonicTransform< Real > xForm;
		xForm.resize( Resolution.value );
		if( GEDT.set ) SetSphereKeys( gedt , xForm , gedtKey , Threads.value );
		else SetSphereKeys( sqr_edt , xForm , edtKey , Threads.value ) , SetSphereKeys( raster , xForm , rasterKey , Threads.value );
	}
	if( GEDT.set )
	{
		Real norm = 0;
		for( int i=0 ; i<gedtKey.size() ; i++ ) norm += gedtKey[i].squareNorm();
		norm = Real( sqrt(norm) );
		for( int i=0 ; i<gedtKey.size() ; i++ ) gedtKey[i] /= norm;
	}
	else
	{
		Real norm = raster.squareNorm();
		for( int i=0 ; i<rasterKey.size() ; i++ ) rasterKey[i] /= norm , edtKey[i] /= Real(Resolution.value * Resolution.value);
	}
	if( Verbose.set ) printf( "\t\tHarmonic Key Time: %.2f(s)\n" , Time()-t );

	RotationGrid< Real > so3Grid1 , so3Grid2;
	t = Time();
	{
		WignerDTransform< Real > xFormSO3;
		FourierKeySO3< Real > key1 , key2;
		xFormSO3.resize( Resolution.value );
		key1.resize( Resolution.value ) , key2.resize( Resolution.value );
		if( GEDT.set ) AddWignerDCoefficients(   gedtKey , gedtKey , key1 , key2 , Threads.value );
		else           AddWignerDCoefficients( rasterKey ,  edtKey , key1 , key2 , Threads.value );
		// Since we are correlating a shape with itself, the symmetric contribution is the same,
		// as it just reorders the group elements before summing the dot-products
		xFormSO3.InverseFourier( key1 , so3Grid1 ) , xFormSO3.InverseFourier( key2 , so3Grid2 );
	}
	if( Verbose.set ) printf( "\t\tWigner-D Time: %.2f(s)\n" , Time()-t );

	t = Time();
	{
		axialSymmetry.resize( Resolution.value );
		refSymmetry.resize( Resolution.value );
		for( int i=0 ; i<rotSymmetry.size() ; i++ ) rotSymmetry[i].resize( Resolution.value );
#pragma omp parallel for num_threads( Threads.value )
		for( int ij=0 ; ij<(Resolution.value/2)*Resolution.value ; ij++ )
		{
			int i = ij / Resolution.value , j = ij % Resolution.value;
			int _i = i+Resolution.value/2 , _j = Resolution.value-1-j;
			Point3D< Real > p;
			refSymmetry.setCoordinates( i , j , &p[0] );
			{
				int d=2;
				Real temp=0;
				for( int k=0 ; k<d ; k++ )
				{
					SquareMatrix< Real , 3 > rot = RotationMatrix< Real >( p , Real(2.*M_PI/d)*k );
					Real x , y , z , R[3][3];
					for( int ii=0 ; ii<3 ; ii++ ) for( int jj=0 ; jj<3 ; jj++ ) R[ii][jj] = rot(ii,jj);
					so3Grid2.setCoordinates( R , x , y , z );
					if( k%2 ) temp += so3Grid2.sample( x , y , z );
					else      temp += so3Grid1.sample( x , y , z );
				}
				refSymmetry(i,j) = refSymmetry(_i,_j) = temp/d;
			}
			for( int d=2 ; d<rotSymmetry.size()+2 ; d++ )
			{
				Real temp=0;
				for( int k=0 ; k<d ; k++ )
				{
					SquareMatrix< Real , 3 > rot = RotationMatrix< Real >( p , Real(2.*M_PI/d)*k );
					Real x , y , z , R[3][3];
					for( int ii=0 ; ii<3 ; ii++ ) for( int jj=0 ; jj<3 ; jj++ ) R[ii][jj] = rot(ii,jj);
					so3Grid1.setCoordinates( R , x , y , z );
					temp += so3Grid1.sample( x , y , z );
				}
				rotSymmetry[d-2](i,j) = rotSymmetry[d-2](_i,_j) = temp/d;
			}
			{
				int d=Resolution.value;
				Real temp=0;
				for( int k=0 ; k<d ; k++ )
				{
					SquareMatrix< Real , 3 > rot = RotationMatrix< Real >( p , Real(2.*M_PI/d)*k );
					Real x , y , z , R[3][3];
					for( int ii=0 ; ii<3 ; ii++ ) for( int jj=0 ; jj<3 ; jj++ ) R[ii][jj] = rot(ii,jj);
					so3Grid1.setCoordinates( R , x , y , z );
					temp += so3Grid1.sample( x , y , z );
				}
				axialSymmetry(i,j) = axialSymmetry(_i,_j) = temp/d;
			}
		}
	}
	if( Verbose.set ) printf( "\t\tSymmetry Sampling Time: %.2f(s)\n" , Time()-t );
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
	SquareMatrix< Real , 4 > xForm;
	std::vector< Point3D< Real > > vertices;
	std::vector< TriangleIndex > triangles;
	{
		int fileType;
		std::vector< PlyColorVertex< float > > _vertices;
		PlyReadTriangles( In.value , _vertices , triangles , PlyColorVertex< float >::ReadProperties , NULL , PlyColorVertex< float >::ReadComponents , fileType );
		vertices.resize( _vertices.size() );
		for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] = Point3D< Real >( _vertices[i].point );
		XForm( vertices , GetAligningXForm( vertices , triangles , Real(MomentRadiusScale.value) , 0 ) , Threads.value );
	}

	SphericalGrid< Real > axialSymmetry , refSymmetry;
	std::vector< SphericalGrid< Real > > rotSymmetry( std::max< int >( 0 , MaxRotationalSymmetry.value-1 ) );
	double t = Time();
	{
#pragma omp parallel for num_threads( Threads.value )
		for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] = vertices[i] * Real(0.5) + Point3D< Real >( Real(0.5) , Real(0.5) , Real(0.5) );
		_main_( vertices , triangles , axialSymmetry , refSymmetry , rotSymmetry );
#pragma omp parallel for num_threads( Threads.value )
		for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] = ( vertices[i] - Point3D< Real >( Real(0.5) , Real(0.5) , Real(0.5) ) ) * Real(2.);
	}
	if( Verbose.set ) printf( "\tSymmetry Time: %.2f(s)\n" , Time()-t );

	if( OutHeader.set )
	{
		char fileName[512];
		sprintf( fileName ,   "%s.ref.sgrid" , OutHeader.value ) ,   refSymmetry.write( fileName );
		sprintf( fileName , "%s.axial.sgrid" , OutHeader.value ) , axialSymmetry.write( fileName );
		for( int i=0 ; i<rotSymmetry.size() ; i++ ) sprintf( fileName , "%s.rot.%d.sgrid" , OutHeader.value , i+2 ) , rotSymmetry[i].write( fileName );
	}
	if( Verbose.set )
	{
		Real min , max;
		{
			min = max = axialSymmetry(0,0);
			for( int i=0 ; i<Resolution.value ; i++ ) for( int j=0 ; j<Resolution.value ; j++ ) min = std::min< Real >( min , axialSymmetry(i,j) ) , max = std::max< Real >( max , axialSymmetry(i,j) );
			printf( "        Axial: [%f,%f]\n" , min , max );
		}
		{
			min = max = refSymmetry(0,0);
			for( int i=0 ; i<Resolution.value ; i++ ) for( int j=0 ; j<Resolution.value ; j++ ) min = std::min< Real >( min , refSymmetry(i,j) ) , max = std::max< Real >( max , refSymmetry(i,j) );
			printf( "   Reflective: [%f,%f]\n" , min , max );
		}
		for( int d=0 ; d<rotSymmetry.size() ; d++ )
		{
			min = max = rotSymmetry[d](0,0);
			for( int i=0 ; i<Resolution.value ; i++ ) for( int j=0 ; j<Resolution.value ; j++ ) min = std::min< Real >( min , rotSymmetry[d](i,j) ) , max = std::max< Real >( max , rotSymmetry[d](i,j) );
			printf( "Rotational[%d]: [%f,%f]\n" , d+2 , min , max );
		}
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

	int ret;
	double t = Time();
	if( Double.set ) ret = _main_< double >( );
	else             ret = _main_< float  >( );
	if( Verbose.set ) printf( "Running Time: %.2f(s)\n" , Time()-t );

	return ret;
}