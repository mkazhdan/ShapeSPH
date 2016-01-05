#ifndef RASTERIZER_INCLUDED
#define RASTERIZER_INCLUDED

#include <vector>
#include <omp.h>
#include "Util/Geometry.h"
#include "SignalProcessing/CubeGrid.h"


template< class Real >
void Rasterize( Point3D< Real > v1 ,                                           CubeGrid< char >& grid , Real scale=Real(1.) );
template< class Real >
void Rasterize( Point3D< Real > v1 , Point3D< Real > v2 ,                      CubeGrid< char >& grid , Real scale=Real(1.) );
template< class Real >
void Rasterize( Point3D< Real > v1 , Point3D< Real > v2 , Point3D< Real > v3 , CubeGrid< char >& grid , Real scale=Real(1.) );

template< class Real >
void Rasterize( const std::vector< Point3D< Real > >& vertices , const std::vector< TriangleIndex >& triangles , CubeGrid< char >& grid , Real scale=Real(1.) , int threads=1 );
template< class Real >
void RasterizeEdges( const std::vector< Point3D< Real > >& vertices , const std::vector< TriangleIndex >& triangles , CubeGrid< char >& grid , Real scale=Real(1.) , int threads=1 );

///////////////////////////////
// Rasterization definitions //
///////////////////////////////
template< class Real >
void Rasterize( const std::vector< Point3D< Real > >& vertices , const std::vector< TriangleIndex >& triangles , CubeGrid< char >& grid , Real scale , int threads )
{
	if( !grid.resolution() )
	{
		fprintf( stderr , "[WARNING] Cannot rasterize triangles to grid of resolution zero\n" );
		return;
	}
#pragma omp parallel for num_threads( threads )
	for( int i=0 ; i<triangles.size() ; i++ )
		Rasterize< Real >( vertices[ triangles[i][0] ] * Real( grid.resolution() ) , vertices[ triangles[i][1] ] * Real( grid.resolution() ) , vertices[ triangles[i][2] ] * Real( grid.resolution() ) , grid , scale );
}
template< class Real >
void RasterizeEdges( const std::vector< Point3D< Real > >& vertices , const std::vector< TriangleIndex >& triangles , CubeGrid< char >& grid , Real scale , int threads )
{
	if( !grid.resolution() )
	{
		fprintf( stderr , "[WARNING] Cannot rasterize triangle edges to grid of resolution zero\n" );
		return;
	}
#pragma omp parallel for num_threads( threads )
	for( int i=0 ; i<triangles.size() ; i++ ) for( int j=0 ; j<3 ; j++ )
		Rasterize< Real >( vertices[ triangles[i][j] ] * Real( grid.resolution() ) , vertices[ triangles[i][(j+1)%3] ] * Real( grid.resolution() ) , grid , scale );
}
template< class Real >
void Rasterize( Point3D< Real > v1 , Point3D< Real > v2 , Point3D< Real > v3 , CubeGrid< char >& grid , Real scale )
{
	if( !grid.resolution() )
	{
		fprintf( stderr , "[WARNING] Cannot rasterize triangleto grid of resolution zero\n" );
		return;
	}

	Point3D< Real > w1 = v2-v1 , w2 = v3-v1;
	Real l1 = Point3D< Real >::SquareNorm( w1 ) , l2 = Point3D< Real >::SquareNorm( w2 );
	if( !l1 && !l2 ) // The triangle is a point
		Rasterize< Real >( v1 ,      grid , scale );
	else if( !l1 ) // The triangle is the line from v1 to v3
		Rasterize< Real >( v1 , v3 , grid , scale );
	else if( !l2 ) // The triangle is the line from v1 to v2
		Rasterize< Real >( v1 , v2 , grid , scale );
	else
	{
		w1 /= Real( sqrt( double( l1 ) ) );
		Real dot = Point3D< Real >::Dot( w1 , w2 );
		Real l = Real( sqrt( l2 - dot*dot ) );
		int steps = (int)( l + 3 );
		steps = (int)( steps*scale );
		for( int i=0 ; i<steps ; i++ )
		{
			Real t = Real(i)/(steps-1);
			Rasterize< Real >( v2*(1-t) + v1*t , v2*(1-t) + v3*t , grid , scale );
		}
	}
}
template< class Real >
void Rasterize( Point3D< Real > v1 , Point3D< Real > v2 , CubeGrid< char >& grid , Real scale )
{
	Point3D< Real > w = v2-v1;
	Real l = Real( sqrt( double( Point3D< Real >::SquareNorm( w ) ) ) );
	int steps = (int)( l+3 );
	steps = (int)( steps*scale );
	for( int i=0 ; i<steps ; i++ )
	{
		Real t = Real(i)/(steps-1);
		Rasterize< Real >( v1*t + v2*(1-t) , grid );
	}
}
template< class Real >
void Rasterize( Point3D< Real > v , CubeGrid< char >& grid , Real scale )
{
	int x = (int)( v[0]+0.5 ) , y = (int)( v[1]+0.5 ) , z = (int)( v[2]+0.5 );
	if( x>=0 && y>=0 && z>=0 && x<grid.resolution() && y<grid.resolution() && z<grid.resolution() ) grid(x,y,z)=1;
}

#endif // RASTERIZER_INCLUDED