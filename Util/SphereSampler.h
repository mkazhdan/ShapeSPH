#ifndef SPHERE_SAMPLER_INCLUDED
#define SPHERE_SAMPLER_INCLUDED

#include <vector>
#include <omp.h>
#include "Util/Geometry.h"
#include "SignalProcessing/CubeGrid.h"
#include "SignalProcessing/SphericalGrid.h"
#include "SignalProcessing/Fourier.h"


template< class Real >
void SampleSpheres( const CubeGrid< Real >& grid , std::vector< SphericalGrid< Real > >& spheres , Point3D< Real > center , Real maxRadius , int radii , int sphereResolution , int threads=1 );
template< class Real >
void SampleSpheres( const CubeGrid< Real >& grid , std::vector< FourierKeyS2< Real > >& sphericalHarmonics , Point3D< Real > center , Real maxRadius , int radii , int sphereResolution , int threads=1 );

template< class Real >
void SubSampleSpheres( const CubeGrid< Real >& grid , std::vector< SphericalGrid< Real > >& spheres , Point3D< Real > center , Real maxRadius , int radii , int sphereResolution , int subSphereResolution , int threads=1 );
template< class Real >
void SubSampleSpheres( const CubeGrid< Real >& grid , std::vector< FourierKeyS2< Real > >& sphericalHarmonics , Point3D< Real > center , Real maxRadius , int radii , int sphereResolution , int subSphereResolution , int threads=1 );


template< class Real >
void SampleSpheres( const std::vector< SphericalGrid< Real > >& spheres , CubeGrid< Real >& grid , Point3D< Real > center , Real maxRadius , int gridResolution , int threads=1 );

///////////////////////////////
// SampleSpheres definitions //
///////////////////////////////
template< class Real >
void SampleSpheres( const CubeGrid< Real >& grid , std::vector< SphericalGrid< Real > >& spheres , Point3D< Real > center , Real maxRadius , int radii , int sphereResolution , int threads )
{
	spheres.resize( radii );
	for( int i=0 ; i<radii ; i++ )
	{
		Real radius = ( Real(i+0.5)/radii ) * maxRadius;
		spheres[i].resize( sphereResolution );
		grid.SphereSample( &center[0] , radius , spheres[i] , threads );
		Real scale = Real( sqrt( 4*M_PI*radius*radius ) );
		Real* _sphere = spheres[i][0];
#pragma omp parallel for num_threads( threads )
		for( int j=0 ; j<sphereResolution*sphereResolution ; j++ ) _sphere[j] *= scale;
	}
}
template< class Real >
void SampleSpheres( const CubeGrid< Real >& grid , std::vector< FourierKeyS2< Real > >& sphericalHarmonics , Point3D< Real > center , Real maxRadius , int radii , int sphereResolution , int threads )
{
	HarmonicTransform< Real > xForm( sphereResolution );
	SphericalGrid< Real > sphere( sphereResolution );
	sphericalHarmonics.resize( radii );
	for( int i=0 ; i<radii ; i++ )
	{
		Real radius = ( Real(i+0.5)/radii ) * maxRadius;
		sphericalHarmonics[i].resize( sphereResolution );
		grid.SphereSample( &center[0] , radius , sphere , threads );
		Real scale = Real( sqrt( 4*M_PI*radius*radius ) );
		Real* _sphere = sphere[0];
#pragma omp parallel for num_threads( threads )
		for( int j=0 ; j<sphereResolution*sphereResolution ; j++ ) _sphere[j] *= scale;
		xForm.ForwardFourier( sphere , sphericalHarmonics[i] );
	}
}
template< class Real >
void SubSampleSpheres( const CubeGrid< Real >& grid , std::vector< SphericalGrid< Real > >& spheres , Point3D< Real > center , Real maxRadius , int radii , int sphereResolution , int subSphereResolution , int threads )
{
	Real thickness = Real( maxRadius/radii );
	spheres.resize( radii );
	for( int i=0 ; i<radii ; i++ )
	{
		Real radius = ( Real(i+0.5)/radii ) * maxRadius;
		spheres[i].resize( sphereResolution );
		grid.SphereSample( &center[0] , radius , spheres[i] , subSphereResolution , thickness , threads );
		Real scale = Real( sqrt( 4*M_PI*radius*radius ) );
		Real* _sphere = spheres[i][0];
#pragma omp parallel for num_threads( threads )
		for( int j=0 ; j<sphereResolution*sphereResolution ; j++ ) _sphere[j] *= scale;
	}
}
template< class Real >
void SubSampleSpheres( const CubeGrid< Real >& grid , std::vector< FourierKeyS2< Real > >& sphericalHarmonics , Point3D< Real > center , Real maxRadius , int radii , int sphereResolution , int subSphereResolution , int threads )
{
	Real thickness = Real( maxRadius/radii );
	HarmonicTransform< Real > xForm( sphereResolution );
	SphericalGrid< Real > sphere( sphereResolution );
	sphericalHarmonics.resize( radii );
	for( int i=0 ; i<radii ; i++ )
	{
		Real radius = ( Real(i+0.5)/radii ) * maxRadius;
		sphericalHarmonics[i].resize( sphereResolution );
		grid.SphereSample( &center[0] , radius , sphere , subSphereResolution , thickness , threads );
		Real scale = Real( sqrt( 4*M_PI*radius*radius ) );
		Real* _sphere = sphere[0];
#pragma omp parallel for num_threads( threads )
		for( int j=0 ; j<sphereResolution*sphereResolution ; j++ ) _sphere[j] *= scale;
		xForm.ForwardFourier( sphere , sphericalHarmonics[i] );
	}
}
template< class Real >
void SampleSpheres( const std::vector< SphericalGrid< Real > >& spheres , CubeGrid< Real >& grid , Point3D< Real > center , Real maxRadius , int gridResolution , int threads )
{
	grid.resize( gridResolution );
#pragma omp parallel for num_threads( threads )
	for( int i=0 ; i<gridResolution ; i++ ) for( int j=0 ; j<gridResolution ; j++ ) for( int k=0 ; k<gridResolution ; k++ )
	{
		Point3D< Real > p = ( Point3D< Real >( Real(i) , Real(j) , Real(k) ) - center ) / maxRadius;
		Real r = Real( sqrt( Point3D< Real >::SquareNorm(p) ) );
		p /= r , r *= spheres.size();
		Real x , y;
		Real scale = Real( sqrt( 4*M_PI*r*r ) );

		spheres[0].setCoordinates( &p[0] , x , y );
		int r1 = int( floor( r ) ) , r2 = r1+1;
		Real dr = Real( 1. - (r-r1) );
		Real temp=0;
		if( r1>=0 && r1<spheres.size() ) temp += spheres[r1](x,y) * Real(   dr);
		if( r2>=0 && r2<spheres.size() ) temp += spheres[r2](x,y) * Real(1.-dr);
		grid( i , j , k ) = temp / scale;
	}
}

#endif // SPHERE_SAMPLER_INCLUDED