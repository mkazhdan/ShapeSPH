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

#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif // M_PI

#ifndef _WIN32
#include <stdlib.h>
#include <string.h>
#endif
#include <float.h>
#include <hash_map>

inline long long HalfEdgeKey( int i1 , int i2 )
{
	return ( ( (long long) i1 )<<32 ) | ( (long long) i2 );
}
inline long long EdgeKey( int i1 , int i2 )
{
	if( i1>i2 ) return HalfEdgeKey( i1 , i2 );
	else		return HalfEdgeKey( i2 , i1 );
}
inline void FactorEdgeKey( long long key , int& idx1 , int& idx2 )
{
    long long i1 , i2;
    i1 = key>>32;
    i2 = (key<<32)>>32;
    idx1 = int( i1 );
    idx2 = int( i2 );
}


///////////
// Point //
///////////
template<class Real,int Dim>
void Point<Real,Dim>::Add		(const Point<Real,Dim>& p)	{	for(int d=0;d<Dim;d++)	coords[d]+=p.coords[d];	}
template<class Real,int Dim>
void Point<Real,Dim>::Scale		(Real s)					{	for(int d=0;d<Dim;d++)	coords[d]*=s;	}
template<class Real,int Dim>
Real Point<Real,Dim>::InnerProduct(const Point<Real,Dim>& p)	const
{
	Real dot=0;
	for(int i=0;i<Dim;i++)	dot+=p.coords[i]*coords[i];
	return dot;
}


/////////////
// Point3D //
/////////////
template< class Real >
Point3D<Real> Point3D<Real>::CrossProduct( const Point3D<Real>& p1 , const Point3D<Real> & p2 )
{
	Point3D<Real> p;
	p.coords[0]= p1.coords[1]*p2.coords[2]-p1.coords[2]*p2.coords[1];
	p.coords[1]=-p1.coords[0]*p2.coords[2]+p1.coords[2]*p2.coords[0];
	p.coords[2]= p1.coords[0]*p2.coords[1]-p1.coords[1]*p2.coords[0];
	return p;
}

////////////
// Matrix //
////////////
template<class Real,int Cols,int Rows>
void Matrix<Real,Cols,Rows>::Add(const Matrix<Real,Cols,Rows>& m)
{
	for(int i=0;i<Cols;i++)	for(int j=0;j<Rows;j++)	coords[i][j]+=m.coords[i][j];
}
template<class Real,int Cols,int Rows>
void Matrix<Real,Cols,Rows>::Scale(Real s)
{
	for(int i=0;i<Cols;i++)	for(int j=0;j<Rows;j++)	coords[i][j]*=s;
}
template<class Real,int Cols,int Rows>
Real Matrix<Real,Cols,Rows>::InnerProduct(const Matrix<Real,Cols,Rows>& m) const
{
	Real dot=0;
	for(int i=0;i<Cols;i++)
		for(int j=0;j<Rows;j++)
			dot+=m.coords[i][j]*coords[i][j];
	return dot;
}
template<class Real,int Cols,int Rows>
template<int Cols1>
Matrix<Real,Cols1,Rows> Matrix<Real,Cols,Rows>::operator * (const Matrix<Real,Cols1,Cols>& m) const
{
	Matrix<Real,Cols1,Rows> n;
	for(int i=0;i<Cols1;i++)
		for(int j=0;j<Rows;j++)
			for(int k=0;k<Cols;k++)
				n.coords[i][j]+=m.coords[i][k]*coords[k][j];
	return n;
}
template<class Real,int Cols,int Rows>
template<class Real2>
Point<Real2,Rows> Matrix<Real,Cols,Rows>::operator () (const Point<Real2,Cols>& v) const	{	return (*this)*v;	}
template<class Real,int Cols,int Rows>
template<class Real2>
Point<Real2,Rows> Matrix<Real,Cols,Rows>::operator * (const Point<Real2,Cols>& v) const
{
	Point<Real2,Rows> out;
	for(int i=0;i<Rows;i++)
		for(int j=0;j<Cols;j++)
			out.coords[i] += Real2( coords[j][i] ) * v.coords[j];
	return out;
}
template<class Real,int Cols,int Rows>
Matrix<Real,Rows,Cols> Matrix<Real,Cols,Rows>::transpose(void) const
{
	Matrix<Real,Rows,Cols> out;
	for(int i=0;i<Cols;i++)
		for(int j=0;j<Rows;j++)
			out.coords[j][i]=coords[i][j];
	return out;
}

//////////////////
// SquareMatrix //
//////////////////
template<>
inline double SquareMatrix<double,1>::determinant(void)	const	{return coords[0][0];}

template<class Real,int Dim>
Real SquareMatrix<Real,Dim>::subDeterminant(int c,int r)	const
{
	SquareMatrix<double,Dim-1> temp;
	int ii=0;
	for(int i=0;i<Dim;i++)
	{
		if(i==c)	continue;
		int jj=0;
		for(int j=0;j<Dim;j++)
		{
			if(j==r)	continue;
			temp.coords[ii][jj]=this->coords[i][j];
			jj++;
		}
		ii++;
	}
	return Real( temp.determinant() );
}

template<class Real,int Dim>
Real SquareMatrix<Real,Dim>::determinant(void)	const
{
	Real det=0;
	for(int i=0;i<Dim;i++)
		if((i&1)==0)	det+=this->coords[i][0]*subDeterminant(i,0);
		else			det-=this->coords[i][0]*subDeterminant(i,0);
	return det;
}
template< class Real , int  Dim >
Real SquareMatrix< Real , Dim >::trace( void ) const
{
	Real tr = 0;
	for( int i=0 ; i<Dim ; i++ ) tr += coords[i][i];
	return tr;
}
template<class Real,int Dim>
SquareMatrix<Real,Dim> SquareMatrix<Real,Dim>::inverse(void)	const
{
	SquareMatrix iXForm;
	Real d=determinant();
	for(int i=0;i<Dim;i++)
		for(int j=0;j<Dim;j++)
			if(((i+j)&1)==0)	iXForm.coords[j][i]= subDeterminant(i,j)/d;
			else				iXForm.coords[i][j]=-subDeterminant(j,i)/d;
	return iXForm;
}
template<class Real,int Dim>
void SquareMatrix<Real,Dim>::Multiply (const SquareMatrix<Real,Dim>& m)
{
	SquareMatrix temp=*this;
	for(int i=0;i<Dim;i++)
		for(int j=0;j<Dim;j++)
		{
			this->coords[i][j]=0;
			for(int k=0;k<Dim;k++)	this->coords[i][j]+=temp.coords[k][j]*m.coords[i][k];
		}
}
template<class Real,int Dim>
void SquareMatrix<Real,Dim>::SetIdentity(void)
{
	memset(this->coords,0,sizeof(Real)*Dim*Dim);
	for(int i=0;i<Dim;i++)	this->coords[i][i]=1;
}

template<class Real>
Real Random( void )
{
	long long foo = ( (long long) rand() ) * RAND_MAX + rand();
	return Real( ( double(foo) / RAND_MAX) / RAND_MAX );
}

template<class Real>
Real Random2( void )
{
	long long temp = (long long) ( rand() )*RAND_MAX+rand();
	return Real( (double(temp)/RAND_MAX)/RAND_MAX );
}

template< class Real >
Point3D<Real> RandomBallPoint(void){
	Point3D<Real> p;
	while( 1 )
	{
		p.coords[0]=Real(1.0-2.0*Random2<Real>());
		p.coords[1]=Real(1.0-2.0*Random2<Real>());
		p.coords[2]=Real(1.0-2.0*Random2<Real>());
		double l=SquareLength(p);
		if( l<=1 ){return p;}
	}
}
template< class Real >
Point3D<Real> RandomSpherePoint( void )
{
	Point3D<Real> p = RandomBallPoint<Real>();
	Real l = Point3D< Real >::Length(p);
	p.coords[0] /= l;
	p.coords[1] /= l;
	p.coords[2] /= l;
	return p;
}

template< class Real >
SquareMatrix< Real , 3 > RotationMatrix( const Point3D<Real>& axis , const Real& angle )
{
	double a = cos( angle / 2 );
	double b , c , d;
	Point3D< Real > ax = axis * Real(sin( angle / 2 ) / Point3D< Real >::Length( axis ));
	b = ax[0] , c = ax[1] , d = ax[2];
	return RotationMatrix< Real >( Real( a ) , Real( b ) , Real( c ) , Real( d ) );
}
template<class Real>
SquareMatrix< Real , 3 > RotationMatrix( Real a , Real b , Real c , Real d )
{
	SquareMatrix< Real , 3 > rot;
	rot( 0 , 0 ) = 1 - 2*c*c - 2*d*d;
	rot( 1 , 0 ) = 2*b*c - 2*a*d;
	rot( 2 , 0 ) = 2*b*d + 2*a*c;
	rot( 0 , 1 ) = 2*b*c + 2*a*d;
	rot( 1 , 1 ) = 1 - 2*b*b - 2*d*d;
	rot( 2 , 1 ) = 2*c*d - 2*a*b;
	rot( 0 , 2 ) = 2*b*d - 2*a*c;
	rot( 1 , 2 ) = 2*c*d + 2*a*b;
	rot( 2 , 2 ) = 1 - 2*b*b - 2*c*c;
	return rot;
}

template<class Real>
SquareMatrix< Real , 3 > RandomRotationMatrix( void )
{
	Point3D< Real > axis = RandomSpherePoint< Real > ( );
	Real angle = Real( 2.0 * M_PI * Random< Real > ( ) );
	return RotationMatrix( axis , angle );
}

template<class Real>
Real SubDeterminant( const SquareMatrix< Real , 4 >& xForm , int c1 , int r1 , int c2 , int r2 )
{
	return xForm.coords[c1][r1]*xForm.coords[c2][r2]-xForm.coords[c1][r2]*xForm.coords[c2][r1];
}
template<class Real>
Real SubDeterminant( const SquareMatrix< Real , 4 >& xForm , int c , int r )
{
	int c1,r1,c2,r2,row;
	Real d=0,sgn=1.0;
	row=0;
	if(row==r){row++;}
	for(int i=0;i<4;i++)
	{
		if(i==c){continue;}
		c1=0;
		while(c1==i || c1==c){c1++;}
		c2=c1+1;
		while(c2==i || c2==c){c2++;}
		r1=0;
		while(r1==row || r1==r){r1++;}
		r2=r1+1;
		while(r2==row || r2==r){r2++;}
		
		d+=sgn*xForm.coords[i][row]*SubDeterminant(xForm,c1,r1,c2,r2);
		sgn*=-1.0;
	}
	return d;
}


template<class Real>
Real Determinant( const SquareMatrix< Real , 4 >& xForm )
{
	Real d=Real(0);
	for(int i=0;i<4;i++)
		if((i&1)==0)	d+=SubDeterminant(xForm,i,0)*xForm.coords[i][0];
		else			d-=SubDeterminant(xForm,i,0)*xForm.coords[i][0];
	return d;
}

//////////////////////////////
// MinimalAreaTriangulation //
//////////////////////////////
template <class Real>
MinimalAreaTriangulation<Real>::MinimalAreaTriangulation(void)
{
	bestTriangulation=NULL;
	midPoint=NULL;
}
template <class Real>
MinimalAreaTriangulation<Real>::~MinimalAreaTriangulation(void)
{
	if(bestTriangulation)
		delete[] bestTriangulation;
	bestTriangulation=NULL;
	if(midPoint)
		delete[] midPoint;
	midPoint=NULL;
}
template <class Real>
void MinimalAreaTriangulation<Real>::GetTriangulation( const std::vector<Point3D<Real> >& vertices , std::vector<TriangleIndex>& triangles )
{
	triangles.resize( vertices.size() - 2 );
	if( vertices.size()==3 )
	{
		triangles[0][0]=0;
		triangles[0][1]=1;
		triangles[0][2]=2;
		return;
	}
	else if( vertices.size()==4 )
	{
		TriangleIndex tIndex[2][2];
		Real area[2];

		area[0]=area[1]=0;

		tIndex[0][0][0]=0;
		tIndex[0][0][1]=1;
		tIndex[0][0][2]=2;
		tIndex[0][1][0]=2;
		tIndex[0][1][1]=3;
		tIndex[0][1][2]=0;

		tIndex[1][0][0]=0;
		tIndex[1][0][1]=1;
		tIndex[1][0][2]=3;
		tIndex[1][1][0]=3;
		tIndex[1][1][1]=1;
		tIndex[1][1][2]=2;

		Point3D<Real> n,p1,p2;
		for(int i=0;i<2;i++)
			for(int j=0;j<2;j++)
				for(int k=0;k<3;k++)
				{
					p1.coords[k]=vertices[tIndex[i][j][1]].coords[k]-vertices[tIndex[i][j][0]].coords[k];
					p2.coords[k]=vertices[tIndex[i][j][2]].coords[k]-vertices[tIndex[i][j][0]].coords[k];
					n = Point3D< Real >::CrossProduct( p1 , p2 );
					area[i] += Point3D< Real >::Length(n);
				}
		if(area[0]>area[1])
		{
			triangles[0]=tIndex[1][0];
			triangles[1]=tIndex[1][1];
		}
		else
		{
			triangles[0]=tIndex[0][0];
			triangles[1]=tIndex[0][1];
		}
		return;
	}

	if(bestTriangulation) delete[] bestTriangulation;
	if(midPoint) delete[] midPoint;
	bestTriangulation=NULL;
	midPoint=NULL;
	size_t eCount=vertices.size();
	bestTriangulation=new double[eCount*eCount];
	midPoint=new int[eCount*eCount];
	for (unsigned int i = 0; i < eCount * eCount; i++)
        bestTriangulation[i] = -1;
	memset(midPoint,-1,sizeof(int)*eCount*eCount);
	GetArea(0,1,vertices);
//	triangles.clear();
	int idx = 0;
//	GetTriangulation(0,1,vertices,triangles);
	GetTriangulation( 0 , 1 , vertices , triangles , idx );
}
template <class Real>
double MinimalAreaTriangulation<Real>::GetArea(const std::vector<Point3D<Real> >& vertices)
{
	if(bestTriangulation)
		delete[] bestTriangulation;
	if(midPoint)
		delete[] midPoint;
	bestTriangulation=NULL;
	midPoint=NULL;
	size_t eCount=vertices.size();
	bestTriangulation=new double[eCount*eCount];
	midPoint=new int[eCount*eCount];
	for(int i=0;i<eCount*eCount;i++)
		bestTriangulation[i]=-1;
	memset(midPoint,-1,sizeof(int)*eCount*eCount);
	return GetArea(0,1,vertices);
}
template<class Real>
void MinimalAreaTriangulation<Real>::GetTriangulation( const int& i , const int& j , const std::vector<Point3D<Real> >& vertices , std::vector<TriangleIndex>& triangles , int& idx )
{
	TriangleIndex tIndex;
	size_t eCount=vertices.size();
	int ii=i;
	if( i<j ) ii+=(int)eCount;
	if( j+1>=ii ) return;
	ii=midPoint[i*eCount+j];
	if(ii>=0)
	{
		tIndex[0]=i;
		tIndex[1]=j;
		tIndex[2]=ii;
		triangles[idx++] = tIndex;
		GetTriangulation( i , ii , vertices , triangles , idx );
		GetTriangulation( ii , j , vertices , triangles , idx );
	}
}

template<class Real>
double MinimalAreaTriangulation<Real>::GetArea(const int& i,const int& j,const std::vector<Point3D<Real> >& vertices)
{
	double a=FLT_MAX,temp;
	size_t eCount=vertices.size();
	size_t idx=i*eCount+j;
	size_t ii=i;
	if(i<j) ii+=eCount;
	if(j+1>=(int) ii)
	{
		bestTriangulation[idx]=0;
		return 0;
	}
	int mid=-1;
	for(unsigned int r=j+1;r<ii;r++)
	{
		int rr=r%eCount;
		size_t idx1=i*eCount+rr,idx2=rr*eCount+j;
		Point3D<Real> p,p1,p2;
		for(int k=0;k<3;k++)
		{
			p1.coords[k]=vertices[i].coords[k]-vertices[rr].coords[k];
			p2.coords[k]=vertices[j].coords[k]-vertices[rr].coords[k];
		}
		p = Point3D< Real >::CrossProduct( p1 , p2 );
		temp=Point3D< Real >::Length(p);

		if(bestTriangulation[idx1]>0)
		{
			temp+=bestTriangulation[idx1];
			if(temp>a)
				continue;
			if(bestTriangulation[idx2]>0)
				temp+=bestTriangulation[idx2];
			else
				temp+=GetArea(rr,j,vertices);
		}
		else
		{
			if(bestTriangulation[idx2]>0)
				temp+=bestTriangulation[idx2];
			else
				temp+=GetArea(rr,j,vertices);
			if(temp>a)
				continue;
			temp+=GetArea(i,rr,vertices);
		}

		if(temp<a)
		{
			a=temp;
			mid=rr;
		}
	}
	bestTriangulation[idx]=a;
	midPoint[idx]=mid;

	return a;
}

template< class Real >
void SplitTriangle( std::vector< Point3D< Real > >& vertices , TriangleIndex triangle , Point3D< Real > pNormal , Real pOffset ,
				    std::vector< TriangleIndex >& backTriangles , std::vector< TriangleIndex >& frontTriangles )
{
	int bVerts[4] , fVerts[4] , bCount = 0 , fCount = 0;
	Real values[3];
	for( int i=0 ; i<3 ; i++ ) values[i] = Point3D< Real >::Dot( vertices[ triangle[i] ] , pNormal ) - pOffset;
	for( int i=0 ; i<3 ; i++ )
	{
		int i1 = (i+2)%3;
		int i2 = (i+1)%3;
		if( values[i]*values[i1]<0 )
		{
			Real t = values[i] / ( values[i] - values[i1] );
			Point3D< Real > newVert = vertices[ triangle[i1] ]*t + vertices[ triangle[ i ] ]*Real(1.0-t);
			vertices.push_back( newVert );
			bVerts[ bCount++ ] = vertices.size()-1;
			fVerts[ fCount++ ] = vertices.size()-1;
		}
		if( values[i]==0 )
		{
			if( values[i1]<0 || values[i2]<0 ) bVerts[ bCount++ ] = triangle[i];
			if( values[i1]>0 || values[i2]>0 ) fVerts[ fCount++ ] = triangle[i];
		}
		else
			if( values[i]<0 ) bVerts[ bCount++ ] = triangle[i];
			else			  fVerts[ fCount++ ] = triangle[i];
	}
	if( bCount==3 )
	{
		backTriangles.resize( 1 );
		backTriangles[0][0] = bVerts[0];
		backTriangles[0][1] = bVerts[1];
		backTriangles[0][2] = bVerts[2];
	}
	if( bCount==4 )
	{
		backTriangles.resize( 2 );
		backTriangles[0][0] = bVerts[0];
		backTriangles[0][1] = bVerts[1];
		backTriangles[0][2] = bVerts[2];
		backTriangles[1][0] = bVerts[2];
		backTriangles[1][1] = bVerts[3];
		backTriangles[1][2] = bVerts[0];
	}
	if( fCount==3 )
	{
		frontTriangles.resize( 1 );
		frontTriangles[0][0] = fVerts[0];
		frontTriangles[0][1] = fVerts[1];
		frontTriangles[0][2] = fVerts[2];
	}
	if( fCount==4 )
	{
		frontTriangles.resize( 2 );
		frontTriangles[0][0] = fVerts[0];
		frontTriangles[0][1] = fVerts[1];
		frontTriangles[0][2] = fVerts[2];
		frontTriangles[1][0] = fVerts[2];
		frontTriangles[1][1] = fVerts[3];
		frontTriangles[1][2] = fVerts[0];
	}
}
template< class Real >
void SplitTriangle( std::vector< Point3D< Real > >& vertices , TriangleIndex triangle , int direction , Real offset ,
				    std::vector< TriangleIndex >& backTriangles , std::vector< TriangleIndex >& frontTriangles )
{
	int bVerts[4] , fVerts[4] , bCount = 0 , fCount = 0;
	Real values[3];
	for( int i=0 ; i<3 ; i++ ) values[i] = vertices[ triangle[i] ][ direction ] - offset;
	for( int i=0 ; i<3 ; i++ )
	{
		int i1 = (i+2)%3;
		int i2 = (i+1)%3;
		if( values[i]*values[i1]<0 )
		{
			Real t = values[i] / ( values[i] - values[i1] );
			Point3D< Real > newVert = vertices[ triangle[i1] ]*t + vertices[ triangle[ i ] ]*Real(1.0-t);
			vertices.push_back( newVert );
			bVerts[ bCount++ ] = vertices.size()-1;
			fVerts[ fCount++ ] = vertices.size()-1;
		}
		if( values[i]==0 )
		{
			if( values[i1]<0 || values[i2]<0 ) bVerts[ bCount++ ] = triangle[i];
			if( values[i1]>0 || values[i2]>0 ) fVerts[ fCount++ ] = triangle[i];
		}
		else
			if( values[i]<0 ) bVerts[ bCount++ ] = triangle[i];
			else			  fVerts[ fCount++ ] = triangle[i];
	}
	if( bCount==3 )
	{
		backTriangles.resize( 1 );
		backTriangles[0][0] = bVerts[0];
		backTriangles[0][1] = bVerts[1];
		backTriangles[0][2] = bVerts[2];
	}
	if( bCount==4 )
	{
		backTriangles.resize( 2 );
		backTriangles[0][0] = bVerts[0];
		backTriangles[0][1] = bVerts[1];
		backTriangles[0][2] = bVerts[2];
		backTriangles[1][0] = bVerts[2];
		backTriangles[1][1] = bVerts[3];
		backTriangles[1][2] = bVerts[0];
	}
	if( fCount==3 )
	{
		frontTriangles.resize( 1 );
		frontTriangles[0][0] = fVerts[0];
		frontTriangles[0][1] = fVerts[1];
		frontTriangles[0][2] = fVerts[2];
	}
	if( fCount==4 )
	{
		frontTriangles.resize( 2 );
		frontTriangles[0][0] = fVerts[0];
		frontTriangles[0][1] = fVerts[1];
		frontTriangles[0][2] = fVerts[2];
		frontTriangles[1][0] = fVerts[2];
		frontTriangles[1][1] = fVerts[3];
		frontTriangles[1][2] = fVerts[0];
	}
}

template< class Real >
void BarycentricCoordinates( const Point3D< Real >& p , const Point3D< Real >& v1 , const Point3D< Real >& v2, const Point3D< Real >& v3 , Real& a0 , Real& a1 , Real& a2 )
{
	Point3D< Real > p0 =  p - v1;
	Point3D< Real > p1 = v2 - v1;
	Point3D< Real > p2 = v3 - v1;
	Point3D< Real >  n  = Point3D<Real>::CrossProduct( p1 , p2 );
	Point3D< Real > _v1 = Point3D<Real>::CrossProduct( p2 , n  );
	Point3D< Real > _v2 = Point3D<Real>::CrossProduct( p1 , n  );

	_v1 /= Point3D<Real>::Dot( _v1 , p1 );
	_v2 /= Point3D<Real>::Dot( _v2 , p2 );

	a1 = Point3D<Real>::Dot( _v1 , p0 );
	a2 = Point3D<Real>::Dot( _v2 , p0 );
	a0 = Real(1.0) - a1 - a2;
}

