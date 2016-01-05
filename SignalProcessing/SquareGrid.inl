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
#include <fftw3.h>
#include "Complex.h"



template<class Real> SquareGrid< Real >::SquareGrid( void ) : res(0) , values(NULL) { ; }
template<class Real> SquareGrid< Real >::SquareGrid( int r ) : res(0) , values(NULL) { resize(r); }
template<class Real>
SquareGrid<Real>::~SquareGrid(void){ if( values ) resize(0); }
template<class Real>
int SquareGrid<Real>::read(const char* fileName){
	FILE* fp=fopen(fileName,"rb");
	if(!fp){return 0;}
	int r=read(fp);
	fclose(fp);
	return r;
}
template<class Real>
int SquareGrid<Real>::write(const char* fileName) const{
	FILE* fp=fopen(fileName,"wb");
	if(!fp){return 0;}
	int w=write(fp);
	fclose(fp);
	return w;
}
template<class Real>
int SquareGrid<Real>::read(FILE* fp){
	int io,r;
	io=int(fread(&r,sizeof(int),1,fp));
	if(!io){return 0;}
	resize(r);
	io=int(fread(values,sizeof(Real),res*res,fp));
	if(io==res*res){return 1;}
	else{return 0;}
}
template<class Real>
int SquareGrid<Real>::write(FILE* fp) const {
	int io;
	io=int(fwrite(&res,sizeof(int),1,fp));
	if(!io){return 0;}
	io=int(fwrite(values,sizeof(Real),res*res,fp));
	if(io==res*res){return 1;}
	else{return 0;}
}
template<class Real>
int SquareGrid<Real>::resolution(void) const{return res;}

int SquareGrid<float>::resize( int  r)
{
	if( r<0 ) return 0;
	else
	{
		if( values ) fftwf_free( values );
		values = NULL , res = 0;
		if( r )
		{
			values = (float*)fftwf_malloc( sizeof(float)*r*r );
			if(!values) return 0;
			else res=r;
			clear();
		}
		return 1;
	}
}
int SquareGrid<double>::resize( int r )
{
	if( r<0 ) return 0;
	else
	{
		if( values ) fftw_free(values);
		values = NULL , res = 0;
		if( r )
		{
			values = (double*)fftw_malloc( sizeof(double)*r*r );
			if( !values ) return 0;
			else res=r;
			clear();
		}
		return 1;
	}
}
int SquareGrid< Complex<float> >::resize( int r)
{
	if( r<0 ) return 0;
	else
	{
		if( values ) fftwf_free(values);
		values = NULL , res=0;
		if( r )
		{
			values = (Complex<float>*)fftwf_malloc( sizeof(fftwf_complex)*r*r );
			if( !values ) return 0;
			else res=r;
			clear();
		}
		return 1;
	}
}
int SquareGrid< Complex<double> >::resize( int r )
{
	if( r<0 ) return 0;
	else
	{
		if( values ) fftw_free(values);
		values = NULL , res=0;
		if( r )
		{
			values = (Complex<double>*)fftw_malloc( sizeof(fftw_complex)*r*r );
			if( !values ) return 0;
			else res=r;
			clear();
		}
		return 1;
	}
}

template<class Real>
int SquareGrid< Real >::resize( int r )
{
	if( r<0 ) return 0;
	else
	{
		if( values ) delete[] values;
		values = NULL;
		res=0;
		if( r )
		{
			values = new Real[r*r];
			if( !values ) return 0;
			else res=r;
			clear();
		}
		return 1;
	}
}
template<class Real>
void SquareGrid<Real>::clear(void){if(res){memset(values,0,sizeof(Real)*res*res);}}

template<class Real>
Real* SquareGrid<Real>::operator[] (const int& i){return &values[i*res];}
template<class Real>
Real& SquareGrid<Real>::operator() (const int& i,const int& j){
	int x=i,y=j;
	if(x<0){x=res-((-x)%res);}
	x%=res;
	if(y<0){y=res-((-y)%res);}
	y%=res;
	return values[x*res+y];
}
template<class Real>
Real SquareGrid<Real>::operator() (const int& i,const int& j) const {
	int x=i,y=j;
	if(x<0){x=res-((-x)%res);}
	x%=res;
	if(y<0){y=res-((-y)%res);}
	y%=res;
	return values[x*res+y];
}
template<class Real>
Real SquareGrid<Real>::operator() (const double& x,const double& y){
	int x1,y1;
	double dx,dy;

	if(x<0){x1=-int(-x)-1;}
	else{x1=int(x);}
	if(y<0){y1=-int(-y)-1;}
	else{y1=int(y);}

	dx=x-x1;
	dy=y-y1;
	return (*this)(x1,y1)*(Real(1.0)-dx)*(Real(1.0)-dy)+(*this)(x1+1,y1)*dx*(Real(1.0)-dy)+(*this)(x1,y1+1)*(Real(1.0)-dx)*dy+(*this)(x1+1,y1+1)*dx*dy;
}
template<class Real>
Real SquareGrid<Real>::operator() (const double& x,const double& y) const {
	int x1,y1;
	double dx,dy;

	if(x<0){x1=-int(-x)-1;}
	else{x1=int(x);}
	if(y<0){y1=-int(-y)-1;}
	else{y1=int(y);}

	dx=x-x1;
	dy=y-y1;
	return (*this)(x1,y1)*(Real(1.0)-dx)*(Real(1.0)-dy)+(*this)(x1+1,y1)*dx*(Real(1.0)-dy)+(*this)(x1,y1+1)*(Real(1.0)-dx)*dy+(*this)(x1+1,y1+1)*dx*dy;
}
template<class Real>
Real SquareGrid<Real>::squareNorm(void) const{return Dot(*this,*this);}
template<class Real>
Real SquareGrid<Real>::SquareDifference(const SquareGrid& g1,const SquareGrid& g2){return g1.squareNorm()+g2.squareNorm()-2*Dot(g1,g2);}
template<class Real>
Real SquareGrid<Real>::Dot(const SquareGrid& g1,const SquareGrid& g2)
{
	Real d = 0;
	if( g1.res != g2.res ) fprintf(stderr,"Could not compare arrays of different sizes: %d != %d\n",g1.res,g2.res) , exit(0);
	for( int i=0 ; i<g1.res*g1.res ; i++ ) d += g1.values[i]*g2.values[i];
	d /= Real( g1.res * g1.res );
	return Real( 4.* PI * PI * d );
}