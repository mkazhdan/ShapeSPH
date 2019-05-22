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
#include <Include/fftw3.h>
#include <SignalProcessing/Complex.h>



template< class Real >
CubeGrid<Real>::CubeGrid( void ) : res(0) , values(NULL) { ; }
template< class Real >
CubeGrid<Real>::CubeGrid( int res ) : res(0) , values(NULL) { resize(res); }
template<class Real>
CubeGrid<Real>::~CubeGrid(void){if(values){resize(0);}}
template<class Real>
int CubeGrid<Real>::read(const char* fileName){
	FILE* fp=fopen(fileName,"rb");
	if(!fp){return 0;}
	int r=read(fp);
	fclose(fp);
	return r;
}
template<class Real>
int CubeGrid<Real>::write(const char* fileName) const{
	FILE* fp=fopen(fileName,"wb");
	if(!fp){return 0;}
	int w=write(fp);
	fclose(fp);
	return w;
}
template<class Real>
int CubeGrid<Real>::read(FILE* fp){
	int io,r;
	io=int(fread(&r,sizeof(int),1,fp));
	if(!io){return 0;}
	resize(r);
	io=int(fread(values,sizeof(Real),res*res*res,fp));
	if(io==res*res*res){return 1;}
	else{return 0;}
}
template<class Real>
int CubeGrid<Real>::write(FILE* fp) const {
	int io;
	io=int(fwrite(&res,sizeof(int),1,fp));
	if(!io){return 0;}
	io=int(fwrite(values,sizeof(Real),res*res*res,fp));
	if(io==res*res*res){return 1;}
	else{return 0;}
}
template<class Real>
int CubeGrid<Real>::resolution(void) const{return res;}

template<>
int CubeGrid<float>::resize(const int& r){
	if(r<0){return 0;}
	else{
		if(values){fftwf_free(values);}
		values=NULL;
		res=0;
		if(r){
			values=(float*)fftwf_malloc(sizeof(float)*r*r*r);
			if(!values){return 0;}
			else{res=r;}
			clear();
		}
		return 1;
	}
}
template<>
int CubeGrid<double>::resize(const int& r){
	if(r<0){return 0;}
	else{
		if(values){fftw_free(values);}
		values=NULL;
		res=0;
		if(r){
			values=(double*)fftw_malloc(sizeof(double)*r*r*r);
			if(!values){return 0;}
			else{res=r;}
			clear();
		}
		return 1;
	}
}
template<>
int CubeGrid<Complex<float> >::resize(const int& r){
	if(r<0){return 0;}
	else{
		if(values){fftwf_free(values);}
		values=NULL;
		res=0;
		if(r){
			values=(Complex<float>*)fftwf_malloc(sizeof(fftwf_complex)*r*r*r);
			if(!values){return 0;}
			else{res=r;}
			clear();
		}
		return 1;
	}
}
template<>
int CubeGrid<Complex<double> >::resize(const int& r){
	if(r<0){return 0;}
	else{
		if(values){fftw_free(values);}
		values=NULL;
		res=0;
		if(r){
			values=(Complex<double>*)fftw_malloc(sizeof(fftw_complex)*r*r*r);
			if(!values){return 0;}
			else{res=r;}
			clear();
		}
		return 1;
	}
}

template<class Real>
int CubeGrid<Real>::resize(const int& r){
	if(r<0){return 0;}
	else{
		if(values){delete[] values;}
		values=NULL;
		res=0;
		if(r){
			values=new Real[r*r*r];
			if(!values){return 0;}
			else{res=r;}
			clear();
		}
		return 1;
	}
}
template<class Real>
void CubeGrid<Real>::clear(void){if(res){memset(values,0,sizeof(Real)*res*res*res);}}

template< class Real > Real* CubeGrid< Real >::operator[] ( int x ){ return values + res*res*x; }
template< class Real > const Real* CubeGrid< Real >::operator[] ( int x ) const { return values + res*res*x; }
template< class Real >
Real& CubeGrid<Real>::operator() (const int& i,const int& j,const int& k){
	int x=i,y=j,z=k;
	if( x<0 ) x=res-((-x)%res);
	x%=res;
	if (y<0 ) y=res-((-y)%res);
	y%=res;
	if( z<0 ) z=res-((-z)%res);
	z%=res;
	return values[x*res*res+y*res+z];
}
template<class Real>
Real CubeGrid<Real>::operator() (const int& i,const int& j,const int& k) const {
	int x=i,y=j,z=k;
	if(x<0){x=res-((-x)%res);}
	x%=res;
	if(y<0){y=res-((-y)%res);}
	y%=res;
	if(z<0){z=res-((-z)%res);}
	z%=res;
	return values[x*res*res+y*res+z];
}
template<class Real>
Real CubeGrid<Real>::operator() ( const double& x , const double& y , const double& z ) const
{
#if 1
	Real temp = 0;
	
	const int xx = int( floor(x) ) , yy = int( floor(y) ) , zz = int( floor(z) );
#if 1
	const double ex2 = x-(double)xx , ex1 = 1-ex2 , ey2 = y-(double)yy , ey1 = 1-ey2 , ez2 = z-(double)zz , ez1 = 1-ez2;
	const int x1 = xx*res*res , x2 = (xx+1)*res*res;
	const int y1 = yy*res , y2 = (yy+1)*res;
	const int z1 = zz , z2 = zz+1;
	const Real *values1 = values + x1 , *values2 = values + x2;
	const bool inz1 = (zz>=0 && zz<res) , inz2 = (zz+1>=0 && zz+1<res);
	if( xx>=0 && xx<res )
	{
		Real t=0;
		if( yy>=0 && yy<res )
		{
			if( inz1 ) t += Real( values1[y1+z1] * ey1*ez1 );
			if( inz2 ) t += Real( values1[y1+z2] * ey1*ez2 );
		}
		if( yy+1>=0 && yy+1<res )
		{
			if( inz2 ) t += Real( values1[y2+z2] * ey2*ez2 );
			if( inz1 ) t += Real( values1[y2+z1] * ey2*ez1 );
		}
		temp += Real( t*ex1 );
	}
	if( xx+1>=0 && xx+1<res )
	{
		Real t=0;
		if( yy>=0 && yy<res )
		{
			if( inz1 ) t += Real( values2[y1+z1] * ey1*ez1 );
			if( inz2 ) t += Real( values2[y1+z2] * ey1*ez2 );
		}
		if( yy+1>=0 && yy+1<res )
		{
			if( inz2 ) t += Real( values2[y2+z2] * ey2*ez2 );
			if( inz1 ) t += Real( values2[y2+z1] * ey2*ez1 );
		}
		temp += Real( t*ex2 );
	}
#else
	double ex = x-(double)xx , ey = y-(double)yy , ez = z-(double)zz;
	if( xx>=0   && xx<res   && yy>=0   && yy<res   && zz>=0   && zz<res   ) temp += Real( values[(xx  )*res*res+(yy  )*res+(zz  )] * (1-(ex))*(1-(ey))*(1-(ez)) );
	if( xx>=0   && xx<res   && yy>=0   && yy<res   && zz+1>=0 && zz+1<res ) temp += Real( values[(xx  )*res*res+(yy  )*res+(zz+1)] * (1-(ex))*(1-(ey))*(  (ez)) );
	if( xx>=0   && xx<res   && yy+1>=0 && yy+1<res && zz+1>=0 && zz+1<res ) temp += Real( values[(xx  )*res*res+(yy+1)*res+(zz+1)] * (1-(ex))*(  (ey))*(  (ez)) );
	if( xx>=0   && xx<res   && yy+1>=0 && yy+1<res && zz>=0   && zz<res   ) temp += Real( values[(xx  )*res*res+(yy+1)*res+(zz  )] * (1-(ex))*(  (ey))*(1-(ez)) );
	if( xx+1>=0 && xx+1<res && yy>=0   && yy<res   && zz>=0   && zz<res   ) temp += Real( values[(xx+1)*res*res+(yy  )*res+(zz  )] * (  (ex))*(1-(ey))*(1-(ez)) );
	if( xx+1>=0 && xx+1<res && yy>=0   && yy<res   && zz+1>=0 && zz+1<res ) temp += Real( values[(xx+1)*res*res+(yy  )*res+(zz+1)] * (  (ex))*(1-(ey))*(  (ez)) );
	if( xx+1>=0 && xx+1<res && yy+1>=0 && yy+1<res && zz+1>=0 && zz+1<res ) temp += Real( values[(xx+1)*res*res+(yy+1)*res+(zz+1)] * (  (ex))*(  (ey))*(  (ez)) );
	if( xx+1>=0 && xx+1<res && yy+1>=0 && yy+1<res && zz>=0   && zz<res   ) temp += Real( values[(xx+1)*res*res+(yy+1)*res+(zz  )] * (  (ex))*(  (ey))*(1-(ez)) );
#endif
	return temp;
#else

	int x1,y1,z1;
	double dx,dy,dz;

	if(x<0){x1=-int(-x)-1;}
	else{x1=int(x);}
	if(y<0){y1=-int(-y)-1;}
	else{y1=int(y);}
	if(z<0){z1=-int(-z)-1;}
	else{z1=int(z);}

	dx=x-x1 , dy=y-y1 , dz=z-z1;
	return
		Real
		(
			(
				(
					(*this)(x1  ,y1,z1)*(Real(1.0)-dx)+
					(*this)(x1+1,y1,z1)*(          dx)
				) * (Real(1.0)-dy) +
			(
					(*this)(x1  ,y1+1,z1)*(Real(1.0)-dx)+
					(*this)(x1+1,y1+1,z1)*(          dx)
				) * dy
			) * (Real(1.0)-dz) +
			(
				(
					(*this)(x1  ,y1,z1+1)*(Real(1.0)-dx)+
					(*this)(x1+1,y1,z1+1)*(          dx)
				) * (Real(1.0)-dy) +
			(
					(*this)(x1  ,y1+1,z1+1)*(Real(1.0)-dx)+
					(*this)(x1+1,y1+1,z1+1)*(          dx)
				) * dy
			) * dz
		)
		;
#endif
}
template<class Real>
Real CubeGrid<Real>::squareNorm(void) const{return Dot(*this,*this);}
template<class Real>
Real CubeGrid<Real>::SquareDifference(const CubeGrid& g1,const CubeGrid& g2){ return g1.squareNorm()+g2.squareNorm()-2*Dot(g1,g2); }
template<class Real>
Real CubeGrid<Real>::Dot( const CubeGrid& g1,const CubeGrid& g2 )
{
	Real d=0;
	if(g1.res != g2.res){
		fprintf(stderr,"Could not compare arrays of different sizes: %d != %d\n",g1.res,g2.res);
		exit(0);
	}
	for(int i=0;i<g1.res*g1.res*g1.res;i++) d+=g1.values[i]*g2.values[i];
	return Real(d/(g1.res*g1.res*g1.res)*8*PI*PI*PI);
}
template< class Real >
void CubeGrid< Real >::SphereSample( const Real* center , Real radius , SphericalGrid< Real >& sGrid , int threads ) const
{
#pragma omp parallel for num_threads( threads )
	for( int i=0 ; i<sGrid.resolution() ; i++ ) for( int j=0 ; j<sGrid.resolution() ; j++ )
	{
		Real coords[3];
		sGrid.setCoordinates( i , j , coords );
		sGrid( i , j ) = (*this)( center[0] + coords[0]*radius , center[1] + coords[1]*radius , center[2] + coords[2]*radius );
	}
}
template< class Real >
void CubeGrid< Real >::SphereSample( const Real* center , Real radius , SphericalGrid< Real >& sGrid , int subRes , Real thickness , int threads ) const
{
#pragma omp parallel for num_threads( threads )
	for( int i=0 ; i<sGrid.resolution() ; i++ ) for( int j=0 ; j<sGrid.resolution() ; j++ )
	{
		Real value=0 , weight=0;
		for( int ii=0 ; ii<subRes ; ii++ ) for( int jj=0 ; jj<subRes ; jj++ )
		{
			Real _i = Real(i-0.5) + Real(ii+0.5)/subRes , _j = Real(j-0.5) + Real(jj+0.5)/subRes;
			Real c[3];
			sGrid.setCoordinates( _i , _j , c );
			for( int kk=0 ; kk<subRes ; kk++ )
			{
				Real r = radius + thickness * ( Real(kk+0.5)/subRes - Real(0.5) );
				// Note that this weighting is slightly off because it doesn't take into account the collapse near the poles
				value += (*this)( center[0] + c[0]*r , center[1] + c[1]*r , center[2] + c[2]*r ) * r * r , weight += r*r;
			}
		}
		sGrid( i , j ) = value/weight;
	}
}
