/* -*- C++ -*-
Copyright (c) 2006, Michael Kazhdan and Matthew Bolitho
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
#ifndef GEOMETRY_INCLUDED
#define GEOMETRY_INCLUDED
#include <cmath>
#include <cassert>
#include <complex>
#include <vector>
#include <cstring>
#include "Algebra.h"

#define PAN_FIX 1

template<class Real>
Real Random(void);

template< class Real , int Dim >
class Point : public InnerProductSpace< Real , Point< Real , Dim > >
{
public:
    /////////////////////////////////
    // Inner product space methods //
    void Add            ( const Point& p );
    void Scale          ( Real s );
    Real InnerProduct   ( const Point< Real , Dim >& p ) const;
    /////////////////////////////////

    Real coords[Dim];
    Point( void ) { memset( coords , 0 , sizeof(Real)*Dim ); }
	Point( Point< Real , Dim+1 > p ){ for( int i=0 ; i<Dim ; i++ ) coords[i] = p.coords[i]; }
	Point( Point< Real , Dim-1 > p , Real pCoord=Real(1) ){ for( int i=0 ; i<Dim-1 ; i++ ) coords[i] = p.coords[i] ; p[Dim-1] = Real(pCoord); }
    template< class Real2 >
    operator Point< Real2, Dim > ( void ) const
    {
        Point< Real2, Dim > p;
        for( int d=0 ; d<Dim ; d++ ) p.coords[d] = Real2( coords[d] ); 
        return p;
    }
    Real& operator [] (int idx) { return coords[idx]; }
    const Real& operator [] (int idx) const { return coords[idx]; }
};

template<class Real,int Cols,int Rows>
class Matrix : public InnerProductSpace<Real,Matrix<Real,Cols,Rows> >
{
public:
    //////////////////////////
    // Vector space methods //
    void Add            ( const Matrix& m );
    void Scale          ( Real s );
    Real InnerProduct   ( const Matrix& m ) const;
    //////////////////////////

    Real coords[Cols][Rows];
    Matrix ( void ) { memset( coords , 0 , sizeof( Real ) * Cols * Rows ); }
    template<class Real2>
    operator Matrix< Real2 , Cols , Rows > ( void ) const
    {
        Matrix< Real2, Cols , Rows > m;
        for( int c=0 ; c<Cols ; c++ ) for ( int r=0 ; r<Rows ; r++ ) m.coords[c][r] = Real2( coords[c][r] ); 
        return m;
    }
    template<int C,int R>
    Matrix(const Matrix<Real,C,R>& m)
    {
        for(int i=0;i<Cols && i<C;i++)
            for(int j=0;j<Rows && j<R;j++)
                coords[i][j]=m.coords[i][j];
    }
    Real& operator () (int c,int r) { return coords[c][r]; }
    const Real& operator () (int c,int r) const { return coords[c][r]; }

    template<int Cols1>
    Matrix<Real,Cols1,Rows> operator * ( const Matrix< Real , Cols1 , Cols >& m ) const;

    Matrix<Real,Rows,Cols> transpose( void )    const;

    template<class Real2>
    Point<Real2,Rows> operator * ( const Point< Real2 , Cols >& v ) const;
    template<class Real2>
    Point<Real2,Rows> operator () ( const Point< Real2 , Cols >& v ) const;
};

template< class Real , int Dim >
class SquareMatrix : public Algebra< Real , SquareMatrix< Real , Dim > > , public Matrix< Real , Dim , Dim >
{
public:
    ////////////////////////////////
    // Additional algebra methods //
    void Multiply (const SquareMatrix& m);
    void SetIdentity(void);
    ////////////////////////////////

    SquareMatrix( const Matrix< Real , Dim , Dim > & m ) { memcpy( coords , m.coords , sizeof(Real)*Dim*Dim );}
    SquareMatrix( void )                                 { memset( coords , 0        , sizeof(Real)*Dim*Dim );}
    Real subDeterminant(int c,int r) const;
    Real determinant(void) const;
    Real trace(void) const;
    SquareMatrix inverse(void) const;

    using Matrix< Real , Dim , Dim >::operator *;
    using Matrix< Real , Dim , Dim >::operator ();
	using Matrix< Real , Dim , Dim >::coords;
};
template< class V , int Dim , class _R = typename V::R >
class Gradient : public VectorSpace< _R , Gradient< V , Dim , _R > >
{
public:
	//////////////////////////
	// Vector space methods //
    void Add            ( const Gradient& g ) { for( int c=0  ; c<Dim ; c++ ) gradients[c] += g.gradients[c]; }
    void Scale          ( _R s ) { for( int c=0 ; c<Dim ; c++ ) gradients[c] *= s; }
	//                      //
	//////////////////////////

    V gradients[Dim];
    Gradient( void ) { for( int d=0 ; d<Dim ;  d++ ) gradients[d] *= 0; }
    V& operator[] ( int idx ) { return gradients[idx]; }
    const V& operator[] ( int idx ) const { return gradients[idx]; }

    template< class V2 , class _R2>
    operator Gradient< V2, Dim , _R2 > ( void ) const
    {
        Gradient< V2 , Dim , _R2 > g;
        for( int d=0 ; d<Dim ; d++ ) g.gradients[d] = V2( gradients[d] ); 
        return g;
    }

    template< class Real >
    Gradient Project( const Point< Real , Dim >& dir ) const
    {
        V dot;
        Gradient g;
        g *= 0;
        dot *= 0;
        Real len = Real( sqrt( Point< Real , Dim >::SquareNorm( dir ) ) );
        if( !len ) return g;
        Point< Real , Dim > _dir = dir / len;
        for( int d=0 ; d<Dim ; d++ ) dot += gradients[d] * _dir[d];
        for( int d=0 ; d<Dim ; d++ ) g.gradients[d] = dot * _dir[d];
        return g;
    }
};

template< class V , int Dim , class _R = typename V::R >
class ConstantFunction : public VectorSpace< _R , ConstantFunction< V , Dim , _R > >
{
public:
    V value;
    Gradient< V , Dim , _R > gradients;
    ConstantFunction( void ) { value *= 0 , gradients *= 0;}

    template< class Real > V operator( ) ( const Point< Real , Dim >& p ) const { return value; }
    template< class Real > Gradient< V , Dim , _R > gradient( const Point< Real , Dim >& p ) const { return gradients; }

    //////////////////////////
    // Vector space methods //
    void Add            ( const ConstantFunction& cf ) { value += cf.value; }
    void Scale          ( _R s ) { value *= s , this->offset *= s; }
    //////////////////////////
};

template< class V , int Dim , class _R = typename V::R >
class LinearFunction : public VectorSpace< _R , LinearFunction< V , Dim , _R > >
{
public:
    Gradient< V , Dim , _R > gradients;
    V offset;
    LinearFunction( void ) { offset *= 0 ; }
    template< class Real >
    V operator( ) ( const Point< Real , Dim >& p ) const
    {
        V v;
        v *= 0;
        for( int d=0 ; d<Dim ; d++ ) v += gradients[d] * p[d];
        v -= offset;
        return v;
    }
    template< class Real >
    LinearFunction fitToHyperplane( const Point< Real , Dim >& p , const Point< Real , Dim >& n ) const
    {
        LinearFunction f;
        Real len = Point< Real , Dim >::SquareNorm( n );
        if( !len )
        {
            f.gradients *= 0;
            f.offset = -(*this)( p );
        }
        else
        {
            Point< Real , Dim > normal = n / Real( sqrt( double( len ) ) );
            V dot;
            dot *= 0;
            for( int d=0 ; d<Dim ; d++ ) dot += gradients[d] * normal[d];
            for( int d=0 ; d<Dim ; d++ ) f.gradients[d] = gradients[d] - dot * normal[d];
            f.offset *= 0;
            f.offset = -(*this)( p ) + f( p );
        }
        return f;
    }
    template< class V2 , class _R2 >
    operator LinearFunction< V2 , Dim , _R2 > ( void ) const
    {
        LinearFunction< V2 , Dim , _R2 > lf;
        lf.offset = V2 ( offset );
        lf.gradients = Gradient< V2 , Dim , _R2 >( gradients );
        return lf;
    }
    template< class Real >
    Gradient< V , Dim , _R > gradient( const Point< Real , Dim >& p ) const { return gradients; }

    // Warning, this function requires the taking of an inverse, which may fail...
    template< class Real >
    static LinearFunction BestFit( const Point< Real , Dim >* points , const V* values , int count )
    {
        LinearFunction lf;
        V constraint[Dim];
        SquareMatrix< Real , Dim > M , Minv;
        M *= 0;
        for( int d=0 ; d<Dim ; d++ ) constraint[d] *= 0;
        for( int i=0 ; i<count ; i++ )
        {
            for( int k=0 ; k<Dim ; k++ ) for( int l=0 ; l<Dim ; l++ ) M( k , l ) += points[i][k] * points[i][l];
            for( int j=0 ; j<count ; j++ ) for( int k=0 ; k<Dim ; k++ ) for( int l=0 ; l<Dim ; l++ ) M( k , l ) -= points[i][k] * points[j][l] / Real( count ); 

            for( int d=0 ; d<Dim ; d++ ) constraint[d] += values[i] * points[i][d];
            for( int j=0 ; j<count ; j++ ) for( int d=0 ; d<Dim ; d++ ) constraint[d] -= values[j] * points[i][d] / Real( count );
        }
        Minv = M.inverse();

        lf *= 0;
        for( int c=0 ; c<Dim ; c++ ) for( int r=0 ; r<Dim ; r++ ) lf.gradients[r] += constraint[c] * Minv( c , r );
        for( int i=0 ; i<count ; i++ )
        {
            for( int d=0 ; d<Dim ; d++ ) lf.offset += lf.gradients[d] * points[i][d];
            lf.offset -= values[i];
        }
        lf.offset /= Real( count );
        return lf;
    }


    //////////////////////////
    // Vector space methods //
    void Add            ( const LinearFunction& lf ) { this->gradient += lf.gradient , offset += lf.offset; }
    void Scale          ( _R s ) { gradients *= s , offset *= s; }
    //////////////////////////
};

template< class Real , int Dim >
struct OrientedPoint
{
    Point< Real , Dim > position , normal;
    template< class Real2 > operator Point< Real2, Dim > ( void ) const { return Point< Real2 , Dim >( position ); }
};

template< class Real >
class Point2D : public Point<Real,2>
{
public:
    Point2D( void ) : Point< Real , 2 >(){;}
    Point2D( const Point< Real , 2 >& p) { memcpy(this->coords,p.coords,sizeof(Real)*2); }
    Point2D( Real v1 , Real v2 ) { coords[0] = v1 , coords[1] = v2; }

	using Point< Real , 2 >::coords;
};

template< class Real >
class Point3D : public Point< Real , 3 >
{
public:
    using Point<Real, 3>::coords;
    Point3D( void ) : Point<Real,3>(){;}
    Point3D( const Point< Real , 3 >& p) { memcpy( coords , p.coords , sizeof(Real)*3 ); }
    Point3D( Real v1 , Real v2 , Real v3 ) { coords[0] = v1 , coords[1] = v2 , coords[2] = v3; }
    Point3D( const Real *p ){ memcpy( coords , p.coords , sizeof(Real)*3 ); }

    static Point3D CrossProduct( const Point3D& p1 , const Point3D & p2 );
};


template< class V , class _R = typename V::R >
class Gradient3D : public Gradient< V , 3 , _R >
{
public:
    Gradient3D( void ) : Gradient< V , 3 , _R >( ){ ; }
    Gradient3D( const Gradient< V , 3 , _R >& lf ) { for( int d=0 ; d<3 ; d++ ) this->gradients[d] = lf.gradients[d]; }
};

template< class V , class _R = typename V::R >
class LinearFunction3D : public LinearFunction< V , 3 , _R >
{
public:
    LinearFunction3D( void ) : LinearFunction< V , 3 , _R >( ){ ; }
    LinearFunction3D( const LinearFunction< V , 3 , _R >& lf )
    {
        for( int d=0 ; d<3 ; d++ ) this->gradients[d] = lf.gradients[d];
        this->offset = lf.offset;
    }
    template< class Real >
    static LinearFunction3D GetInterpolant( const Point3D< Real >* vertices , const V* values , Point3D< Real > n )
    {
        Point3D< Real > p1 = Point3D< Real >( vertices[1] - vertices[0] );
        Point3D< Real > p2 = Point3D< Real >( vertices[2] - vertices[0] );
        Point3D< Real > v1 = Point3D< Real >::CrossProduct( p2 , n  );
        Point3D< Real > v2 = Point3D< Real >::CrossProduct( p1 , n  );

        Real d1 = Point3D< Real >::Dot( v1 , p1 );
        Real d2 = Point3D< Real >::Dot( v2 , p2 );
        if( !d1 || !d2 )
        {
            LinearFunction3D< V , _R > lf;
            lf.gradients *= 0;
#if PAN_FIX
            lf.offset = -( values[0] + values[1] + values[2] ) / 3;
#else // !PAN_FIX
            lf.offset =  ( values[0] + values[1] + values[2] ) / 3;
#endif // PAN_FIX
            return lf;
        }
        else
        {
            v1 /= d1;
            v2 /= d2;
        }

        LinearFunction3D< V , _R > lf;
        lf.offset =
            - ( values[0] * ( 1.0 + Point3D< Real >::Dot( vertices[0] , v1 ) + Point3D< Real >::Dot( vertices[0] , v2 ) )
            - values[1] * Point3D< Real >::Dot( vertices[0] , v1 )
            - values[2] * Point3D< Real >::Dot( vertices[0] , v2 ) );
        for( int d=0 ; d<3 ; d++ ) lf.gradients[d] = - values[0]*v1[d] - values[0]*v2[d] + values[1]*v1[d] + values[2]*v2[d];
        return lf;
    }
    template< class Real >
    static LinearFunction3D GetInterpolant( const Point3D< Real >& v1 ,
                                          const Point3D< Real >& v2 ,
                                          const Point3D< Real >& v3 ,
                                          const V& s1 , const V& s2 , const V s3 ,
                                          Point3D< Real > n )
    {
        Point3D< Real > p1 = Point3D< Real >( v2 - v1 );
        Point3D< Real > p2 = Point3D< Real >( v3 - v1 );
        Point3D< Real > _v1 = Point3D< Real >::CrossProduct( p2 , n  );
        Point3D< Real > _v2 = Point3D< Real >::CrossProduct( p1 , n  );

        Real d1 = Point3D< Real >::Dot( _v1 , p1 );
        Real d2 = Point3D< Real >::Dot( _v2 , p2 );
        if( !d1 || !d2 )
        {
            LinearFunction3D< V , _R > lf;
            lf.gradients *= 0;
#if PAN_FIX
            lf.offset = -( s1 + s2 + s3 ) / _R(3);
#else // !PAN_FIX
            lf.offset = ( s1 + s2 + s3 ) / _R(3);
#endif // PAN_FIX
            return lf;
        }
        else
        {
            _v1 /= d1;
            _v2 /= d2;
        }

        LinearFunction3D< V , _R > lf;
        lf.offset = - ( s1 * ( Real(1.0)
                        + Point3D< Real >::Dot( v1 , _v1 )
                        + Point3D< Real >::Dot( v1 , _v2 ))
                    - s2 * Point3D< Real >::Dot( v1 , _v1 )
                    - s3 * Point3D< Real >::Dot( v1 , _v2 ) );
        for( int d=0 ; d<3 ; d++ ) lf.gradients[d] = - s1*_v1[d] - s1*_v2[d] + s2*_v1[d] + s3*_v2[d];
        return lf;
    }
    template< class Real >
    static LinearFunction3D GetInterpolant( const Point3D< Real >* vertices , const V* values )
    {
        Point3D< Real > p1 = Point3D< Real >( vertices[1] - vertices[0] );
        Point3D< Real > p2 = Point3D< Real >( vertices[2] - vertices[0] );
        Point3D< Real > n  = Point3D< Real >::CrossProduct( p1 , p2 );
        return GetInterpolant( vertices , values , n );
    }
    template< class Real >
    static LinearFunction3D GetInterpolant( const Point3D< Real >& v1 , const Point3D< Real >& v2 , const Point3D< Real >& v3 , const V& s1 , const V& s2 , const V& s3 )
    {
        Point3D< Real > p1 = Point3D< Real >( v2 - v1 );
        Point3D< Real > p2 = Point3D< Real >( v3 - v1 );
        Point3D< Real > n  = Point3D< Real >::CrossProduct( p1 , p2 );
        return GetInterpolant( v1 ,v2 , v3 , s1 , s2 , s3 , n );
    }
};

template<class Real>
class OrientedPoint2D : public OrientedPoint<Real,2>{;};
template<class Real>
class OrientedPoint3D : public OrientedPoint<Real,3>{;};

template<class Real>
Point3D<Real> RandomBallPoint(void);

template<class Real>
Point3D<Real> RandomSpherePoint(void);

template<class Real>
SquareMatrix< Real , 3 > RotationMatrix( Real a , Real b , Real c , Real d );

template<class Real>
SquareMatrix< Real , 3 > RotationMatrix( const Point3D<Real>& axis , const Real& angle );

template<class Real>
SquareMatrix< Real , 3 > RandomRotationMatrix( void );

template< class Real >
void BarycentricCoordinates( const Point3D< Real >& p , const Point3D< Real >& v1 , const Point3D< Real >& v2, const Point3D< Real >& v3 , Real& a0 , Real& a1 , Real& a2 );

class Edge
{
public:
    double p[2][2];
    double Length( void ) const
	{
        double d[2];
        d[0]=p[0][0]-p[1][0];
        d[1]=p[0][1]-p[1][1];

        return sqrt(d[0]*d[0]+d[1]*d[1]);
    }
};
class Triangle
{
public:
    double p[3][3];
    double Area(void) const{
        double v1[3],v2[3],v[3];
        for(int d=0;d<3;d++){
            v1[d]=p[1][d]-p[0][d];
            v2[d]=p[2][d]-p[0][d];
        }
        v[0]= v1[1]*v2[2]-v1[2]*v2[1];
        v[1]=-v1[0]*v2[2]+v1[2]*v2[0];
        v[2]= v1[0]*v2[1]-v1[1]*v2[0];
        return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])/2;
    }
    double AspectRatio(void) const{
        double d=0;
        int i,j;
        for(i=0;i<3;i++){
      for(i=0;i<3;i++)
            for(j=0;j<3;j++){d+=(p[(i+1)%3][j]-p[i][j])*(p[(i+1)%3][j]-p[i][j]);}
        }
        return Area()/d;
    }
};
class EdgeIndex
{
public:
    int v[2];
    int& operator[] ( int idx ) { return v[idx]; }
    const int& operator[] ( int idx ) const { return v[idx]; }
};

class TriangleIndex
{
protected:
    unsigned int v[3];
public:
    TriangleIndex() { v[0] = v[1] = v[2] = 0; }
    TriangleIndex(unsigned int v0, unsigned int v1, unsigned int v2)
        { v[0] = v0; v[1] = v1; v[2] = v2; }
    unsigned int &operator[](unsigned int idx) { return v[idx]; }
    unsigned int  operator[](unsigned int idx) const { return v[idx]; }
};



template< class Real >
class MinimalAreaTriangulation
{
    double* bestTriangulation;
    int* midPoint;
    double GetArea( const int& i , const int& j , const std::vector< Point3D< Real > >& vertices );
    void GetTriangulation( const int& i , const int& j , const std::vector<Point3D<Real> >& vertices,std::vector<TriangleIndex>& triangles , int& idx);
public:
    MinimalAreaTriangulation(void);
    ~MinimalAreaTriangulation(void);
    double GetArea(const std::vector<Point3D<Real> >& vertices);
    void GetTriangulation( const std::vector<Point3D<Real> >& vertices , std::vector<TriangleIndex>& triangles );
};

template<class Vertex>
class Mesh
{
public:
    std::vector<Vertex> vertices;
    std::vector<std::vector<int> > polygons;
};

template< class Real >
void SplitTriangle( std::vector< Point3D< Real > >& vertices , TriangleIndex triangle , Point3D< Real > pNormal , Real pOffset ,
                    std::vector< TriangleIndex >& backTriangles , std::vector< TriangleIndex >& frontTriangles );
template< class Real >
void SplitTriangle( std::vector< Point3D< Real > >& vertices , TriangleIndex triangle , int direction , Real offset ,
                    std::vector< TriangleIndex >& backTriangles , std::vector< TriangleIndex >& frontTriangles );
#include "Geometry.inl"

#endif // GEOMETRY_INCLUDED
