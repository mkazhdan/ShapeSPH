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
#include <math.h>
#include <float.h>

template<class Real>
Real ArcTan2(const Real& y,const Real& x){
	/* This first case should never happen */
	if( y==0 && x==0 ) return 0;
	if( x==0 ) return y>0 ? Real(  PI/2.0 ) : Real( -PI/2.0 );
	if( x>=0 ) return Real( atan(y/x) );
	else       return y>=0 ? Real( atan(y/x)+PI ) : Real( atan(y/x)-PI );
}

template<class Real>
SphericalGrid<Real>::SphericalGrid(void){
	res=0;
	values=NULL;
}
template<class Real>
SphericalGrid<Real>::SphericalGrid( int r )
{
	res = 0;
	values = NULL;
	if( r ) resize( r );
}
template<class Real>
SphericalGrid<Real>::~SphericalGrid(void){
	if(values){delete[] values;}
	values=NULL;
	res=0;
}
template<class Real>
int SphericalGrid<Real>::read(const char* fileName){
	FILE* fp=fopen(fileName,"rb");
	if(!fp){return 0;}
	int r=read(fp);
	fclose(fp);
	return r;
}
template<class Real>
int SphericalGrid<Real>::write(const char* fileName) const{
	FILE* fp=fopen(fileName,"wb");
	if(!fp){return 0;}
	int w=write(fp);
	fclose(fp);
	return w;
}
template<class Real>
int SphericalGrid<Real>::read(FILE* fp){
	int r,io;
	io=int(fread(&r,sizeof(int),1,fp));
	if(!io){return 0;}
	resize(r);
	io=int(fread(values,sizeof(Real),res*res,fp));
	if(io==res*res){return 1;}
	else{return 0;}
}
template<class Real>
int SphericalGrid<Real>::write(FILE* fp) const {
	int io;
	io=int(fwrite(&res,sizeof(int),1,fp));
	if(!io){return 0;}
	io=int(fwrite(values,sizeof(Real),res*res,fp));
	if(io==res*res){return 1;}
	else{return 0;}
}
template<class Real>
int SphericalGrid<Real>::resolution(void) const{return res;}
template<class Real>
int SphericalGrid<Real>::resize(const int& r){
	if(r<0){return 0;}
	else{
		if(values){delete[] values;}
		values=NULL;
		res=0;
		if(r){
			values=new Real[r*r];
			if(!values){return 0;}
			else{res=r;}
		}
		clear();
		return 1;
	}
}
template<class Real>
void SphericalGrid<Real>::clear(void){if(res){memset(values,0,sizeof(Real)*res*res);}}
template<class Real>
void SphericalGrid<Real>::SetCoordinates(const Real& theta,const Real& phi,Real coords[3]){
	coords[0]=Real(sin(phi)*cos(theta));
	coords[1]=Real(cos(phi));
	coords[2]=Real(sin(phi)*sin(theta));
}
template<class Real>
void SphericalGrid<Real>::SetCoordinates(const Real coords[3],Real& theta,Real& phi){
	if(fabs(coords[0])<16*FLT_EPSILON && fabs(coords[2])<16*FLT_EPSILON){
		if(coords[1]<0){phi=Real(PI);theta=0;}
		else{theta=phi=0;}
	}
	else{
		theta=Real(ArcTan2(coords[2],coords[0]));
		if(theta<0){theta+=Real(2.0*PI);}
		phi=Real(acos(coords[1]));
	}
}


template<class Real>
void SphericalGrid<Real>::setCoordinates( const int& i , const int& j , Real coords[3] ) const
{
	Real theta , phi;
	theta = Real( 2.0*PI*i/res );
	phi = Real( PI*(2.0*j+1)/(2.0*res) );
	SetCoordinates( theta , phi , coords );
}
template<class Real>
void SphericalGrid<Real>::setCoordinates( const Real& i , const Real& j , Real coords[3] ) const
{
	Real theta,phi;
	theta = Real( 2.0*PI*i/res );
	phi = Real( PI*(2.0*j+1)/(2.0*res) );
	SetCoordinates( theta , phi , coords );
}
template<class Real>
void SphericalGrid<Real>::setCoordinates( const Real coords[3] , Real& i , Real& j ) const{
	Real theta,phi;
	SetCoordinates(coords,theta,phi);
	i = Real( theta*res/Real(2.0*PI) );
	j = Real( (phi*Real(2.0*res)/PI-1)/2 );
}

template<class Real>
Real& SphericalGrid<Real>::operator() (const int& i,const int& j){
	int x=j,y=i;

	if(x<0){x=2*res-((-i)%(2*res));}
	x%=2*res;
	if(x>=res){
		x=2*res-x-1;
		y+=res/2;
	}
	if(y<0){y=res-((-y)%res);}
	y=y%res;
	return values[x*res+y];
}
template<class Real>
Real SphericalGrid<Real>::operator() (const int& i,const int& j) const {
	int x=j,y=i;

	if(x<0){x=2*res-((-i)%(2*res));}
	x%=2*res;
	if(x>=res){
		x=2*res-x-1;
		y+=res/2;
	}
	if(y<0){y=res-((-y)%res);}
	y=y%res;
	return values[x*res+y];
}
template<class Real>
Real SphericalGrid<Real>::operator() (const double& x,const double& y){
	int x1,y1;
	Real dx,dy;

	if(x<0){x1=-int(-x)-1;}
	else{x1=int(x);}
	if(y<0){y1=-int(-y)-1;}
	else{y1=int(y);}

	dx = Real( x-x1 );
	dy = Real( y-y1 );
	return (*this)(x1,y1)*(Real(1.0)-dx)*(Real(1.0)-dy)+(*this)(x1+1,y1)*dx*(Real(1.0)-dy)+(*this)(x1,y1+1)*(Real(1.0)-dx)*dy+(*this)(x1+1,y1+1)*dx*dy;
}
template< class Real >
Real SphericalGrid<Real>::operator() (const double& x,const double& y) const
{
	int x1,y1;
	Real dx,dy;

	if(x<0){x1=-int(-x)-1;}
	else{x1=int(x);}
	if(y<0){y1=-int(-y)-1;}
	else{y1=int(y);}

	dx = Real(x-x1);
	dy = Real(y-y1);
	return (*this)(x1,y1)*(Real(1.0)-dx)*(Real(1.0)-dy)+(*this)(x1+1,y1)*dx*(Real(1.0)-dy)+(*this)(x1,y1+1)*(Real(1.0)-dx)*dy+(*this)(x1+1,y1+1)*dx*dy;
}

template<class Real>
Real* SphericalGrid<Real>::operator[] (const int& i){return &values[i*res];}
template<class Real>
Real SphericalGrid<Real>::squareNorm(void) const{return Dot(*this,*this);}
template<class Real>
Real SphericalGrid<Real>::SquareDifference(const SphericalGrid& g1,const SphericalGrid& g2){return g1.squareNorm()+g2.squareNorm()-2*Dot(g1,g2);}
template<class Real>
Real SphericalGrid<Real>::Dot(const SphericalGrid& g1,const SphericalGrid& g2){
	double d=0,c=1;
	if(g1.res != g2.res){
		fprintf(stderr,"Could not compare arrays of different sizes: %d != %d\n",g1.res,g2.res);
		exit(0);
	}
	for(int i=0;i<g1.res;i++){
		double t1=cos(PI*(2.0*i+2)/(2.0*g1.res));
		double t2=(c-t1)/(2*g1.res);
		c=t1;
		for(int j=0;j<g1.res;j++){d+=g1.values[i*g1.res+j]*g2.values[i*g1.res+j]*t2;}
	}
	return Real(d*4*PI);
}

template<class Real>
void SphericalGrid<Real>::Transform(const Real matrix[3][3],const Real* in,Real* out){
	for(int i=0;i<3;i++){
		out[i]=0;
		for(int j=0;j<3;j++){out[i]+=matrix[j][i]*in[j];}
	}
}

template<class Real>
void SphericalGrid<Real>::Rotate(const SphericalGrid& in,const Real rotation[3][3],SphericalGrid& out){
	int res=in.resolution();
	Real rotationTranspose[3][3];
	Real inC[3],outC[3];
	Real ii,jj;
	out.resize(res);

	for(int i=0;i<3;i++){for(int j=0;j<3;j++){rotationTranspose[i][j]=rotation[j][i];}}
	for(int i=0;i<res;i++){
		for(int j=0;j<res;j++){
			out.setCoordinates(i,j,inC);
			Transform(rotationTranspose,inC,outC);
			out.setCoordinates(outC,ii,jj);
			out(i,j)=in(ii,jj);
		}
	}
}
