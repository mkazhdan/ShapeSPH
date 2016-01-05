#ifndef TRIANGLE_MESH_INCLUDED
#define TRIANGLE_MESH_INCLUDED

template< class Real > void XForm( std::vector< Point3D< Real > >& vertices , const SquareMatrix< Real , 3 >& xForm , int threads=1 );
template< class Real > void XForm( std::vector< Point3D< Real > >& vertices , const SquareMatrix< Real , 4 >& xForm , int threads=1 );

template< class Real >
Real Area( const std::vector< Point3D< Real > >& vertices , const std::vector< TriangleIndex >& triangles );

template< class Real >
Point3D< Real > Center( const std::vector< Point3D< Real > >& vertices , const std::vector< TriangleIndex >& triangles );

template< class Real >
SquareMatrix< Real , 3 > CovarianceMatrix( Point3D< Real > v1 , Point3D< Real > v2 , Point3D< Real > v3 , Point3D< Real > center );
template< class Real >
SquareMatrix< Real , 3 > CovarianceMatrix( const std::vector< Point3D< Real > >& vertices , const std::vector< TriangleIndex >& triangles , Point3D< Real > center );

template< class Real >
Real BoundingRadius( const std::vector< Point3D< Real > >& vertices , Point3D< Real > center );
template< class Real >
Real MomentRadius( const std::vector< Point3D< Real > >& vertices , const std::vector< TriangleIndex >& triangles , Point3D< Real > center );

template< class Real >
SquareMatrix< Real , 4 > GetAligningXForm( const std::vector< Point3D< Real > >& vertices , const std::vector< TriangleIndex >& triangles , Real radiusScale , int anisotropic=0 );

template< class Real > std::pair< Point3D< Real > , Point3D< Real > > GetBoundingBox( const std::vector< Point3D< Real > >& vertices );

//////////////////////////////
// TriangleMesh definitions //
//////////////////////////////
#include "Util/lineqn.h"

template< class Real >
void XForm( std::vector< Point3D< Real > >& vertices , const SquareMatrix< Real , 3 >& xForm , int threads )
{
#pragma omp parallel for num_threads( threads )
	for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] = xForm * vertices[i];
}
template< class Real >
void XForm( std::vector< Point3D< Real > >& vertices , const SquareMatrix< Real , 4 >& xForm , int threads )
{
#pragma omp parallel for num_threads( threads )
	for( int i=0 ; i<vertices.size() ; i++ )
	{
		Point< Real , 4 > p;
		for( int j=0 ; j<3 ; j++ ) p[j] = vertices[i][j];
		p[3] = Real(1.);
		p = xForm * p;
		for( int j=0 ; j<3 ; j++ ) vertices[i][j] = p[j];
	}
}
template< class Real >
Real Area( const std::vector< Point3D< Real > >& vertices , const std::vector< TriangleIndex >& triangles )
{
	Real area = 0;
	for( int i=0 ; i<triangles.size() ; i++ )
	{
		Point3D< Real > v[] = { vertices[ triangles[i][0] ] , vertices[ triangles[i][1] ] , vertices[ triangles[i][2] ] };
		area += Real( Point3D< Real >::Length( Point3D< Real >::CrossProduct( v[1]-v[0] , v[2]-v[0] ) ) / 2. );
	}
	return area;
}
template< class Real >
Point3D< Real > Center( const std::vector< Point3D< Real > >& vertices , const std::vector< TriangleIndex >& triangles )
{
	Real area = 0;
	Point3D< Real > center;
	for( int i=0 ; i<triangles.size() ; i++ )
	{
		Point3D< Real > v[] = { vertices[ triangles[i][0] ] , vertices[ triangles[i][1] ] , vertices[ triangles[i][2] ] };
		Point3D< Real > c = ( v[0] + v[1] + v[2] ) / Real(3.);
		Real a = Real( Point3D< Real >::Length( Point3D< Real >::CrossProduct( v[1]-v[0] , v[2]-v[0] ) ) / 2. );
		center += c * a;
		area += a;
	}
	return center / area;
}
template< class Real >
SquareMatrix< Real , 3 > CovarianceMatrix( Point3D< Real > v1 , Point3D< Real > v2 , Point3D< Real > v3 , Point3D< Real > center )
{
	SquareMatrix< Real , 3 > covariance;
	Real a = Point3D< Real >::Length( Point3D< Real >::CrossProduct( v2-v1 , v3-v1 ) );
	Point3D< Real > v[] = { v1-center , v3-v1 , v2-v3 };
	const Real Factors[3][3] =
	{
		{ Real(1./2) , Real(1./3) , Real(1./ 6) } ,
		{ Real(1./3) , Real(1./4) , Real(1./ 8) } ,
		{ Real(1./6) , Real(1./8) , Real(1./12) }
	};

	for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ )
	{
		Real dot = Real(0);
		for( int k=0 ; k<3 ; k++ ) for( int l=0 ; l<3 ; l++ ) dot += ( v[k][i] * v[l][j] ) * Factors[k][l];
		covariance( i , j ) = dot * a;
	}
	return covariance;

}
template< class Real >
SquareMatrix< Real , 3 > CovarianceMatrix( const std::vector< Point3D< Real > >& vertices , const std::vector< TriangleIndex >& triangles , Point3D< Real > center )
{
	SquareMatrix< Real , 3 > covariance;
	for( int i=0 ; i<triangles.size() ; i++ ) covariance += CovarianceMatrix( vertices[ triangles[i][0] ] , vertices[ triangles[i][1] ] , vertices[ triangles[i][2] ] , center );
	return covariance;
}


template< class Real >
Real BoundingRadius( const std::vector< Point3D< Real > >& vertices , Point3D< Real > center )
{
	Real radius2 = 0;
	for( int i=0 ; i<vertices.size() ; i++ ) radius2 = std::max< Real >( radius2 , Point3D< Real >::SquareNorm( vertices[i]-center ) );
	return Real( sqrt( radius2 ) );
}
template< class Real >
Real MomentRadius( const std::vector< Point3D< Real > >& vertices , const std::vector< TriangleIndex >& triangles , Point3D< Real > center )
{
	SquareMatrix< Real , 3 > c = CovarianceMatrix( vertices , triangles , center ) / Area( vertices , triangles );
	return Real( sqrt( c.trace() ) );
}

template< class Real , int Dim >
void PrintMatrix( SquareMatrix< Real , Dim > M )
{
	for( int i=0 ; i<Dim ; i++ )
	{
		for( int j=0 ; j<Dim ; j++ ) printf( " %f" , M(j,i) );
		printf( "\n" );
	}
	printf( "\n" );
}
template< class Real >
SquareMatrix< Real , 4 > GetAligningXForm( const std::vector< Point3D< Real > >& vertices , const std::vector< TriangleIndex >& triangles , Real radiusScale , int anisotropic )
{
	SquareMatrix< Real , 4 > translate=SquareMatrix< Real , 4 >::Identity() , scale=SquareMatrix< Real , 4 >::Identity() , aScale=SquareMatrix< Real , 4 >::Identity();
	std::vector< Point3D< Real > > _vertices = vertices;
	for( int i=0 ; i<anisotropic ; i++ )
	{
		Point3D< Real > center = Center( _vertices , triangles );
		// Rescale initially so that the covariance matrix has reasonable size
		Real radius = BoundingRadius( _vertices , center );
		for( int j=0 ; j<_vertices.size() ; j++ ) _vertices[j] = (_vertices[j]-center ) / radius;

		SquareMatrix< Real , 3 > cov = CovarianceMatrix( _vertices , triangles , Point3D< Real >( Real(0) , Real(0) , Real(0) ) ) , covRoot;
		Real A[3][3] , d[3];
		for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ ) A[j][k] = cov(j,k);
		eigdc< Real , 3 >( A , d );
		{
			SquareMatrix< Real , 3 > cRotate , cScale;
			for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ ) cRotate(j,k) = A[j][k];
			cScale(0,0) = Real( sqrt( d[0] ) ) , cScale(1,1) = Real( sqrt( d[1] ) ) , cScale(2,2) = Real( sqrt( d[2] ) );
			covRoot = SquareMatrix< Real , 3 >( cRotate.transpose() ) * cScale * cRotate;
			covRoot = covRoot.inverse();
		}
		SquareMatrix< Real , 4 > subTranslate , subScale;
		subScale = subTranslate = SquareMatrix< Real , 4 >::Identity();
		for( int j=0 ; j<3 ; j++ ) subTranslate(3,j) = -center[j];
		for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ ) subScale(j,k) = covRoot(j,k)/radius;
		for( int j=0 ; j<_vertices.size() ; j++ ) _vertices[j] = covRoot * _vertices[j];

		aScale = subScale * subTranslate * aScale;
	}
	Point3D< Real > center = Center( _vertices , triangles );
	Real radius;
	if( radiusScale>0 ) radius = MomentRadius( _vertices , triangles , center ) * radiusScale;
	else                radius = BoundingRadius( _vertices , center );
	for( int i=0 ; i<3 ; i++ ) translate(3,i) = -center[i] , scale(i,i) = Real(1./radius);
	return scale * translate * aScale;
}
template< class Real >
std::pair< Point3D< Real > , Point3D< Real > > GetBoundingBox( const std::vector< Point3D< Real > >& vertices )
{
	Point3D< Real > min , max;
	min = max = vertices[0];
	for( int i=0 ; i<vertices.size() ; i++ ) for( int j=0 ; j<3 ; j++ ) min[j] = std::min< Real >( min[j] , vertices[i][j] ) , max[j] = std::max< Real >( max[j] , vertices[i][j] );
	return std::pair< Point3D< Real > , Point3D< Real > >( min , max );
}

#endif // TRIANGLE_MESH_INCLUDED