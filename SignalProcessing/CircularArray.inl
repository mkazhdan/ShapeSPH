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


template< class Real > CircularArray< Real >::CircularArray( void ) : res(0) , values(NULL) { ; }
template< class Real > CircularArray< Real >::CircularArray( int r ) : res(0) , values(NULL) { resize(r); }
template<class Real>
CircularArray<Real>::~CircularArray(void){
	if(values){delete[] values;}
	values=NULL;
	res=0;
}
template<class Real>
int CircularArray<Real>::read(const char* fileName){
	FILE* fp=fopen(fileName,"rb");
	if(!fp){return 0;}
	int r=read(fp);
	fclose(fp);
	return r;
}
template<class Real>
int CircularArray<Real>::write(const char* fileName) const{
	FILE* fp=fopen(fileName,"wb");
	if(!fp){return 0;}
	int w=write(fp);
	fclose(fp);
	return w;
}
template<class Real>
int CircularArray<Real>::read(FILE* fp){
	int io,r;
	io=int(fread(&r,sizeof(int),1,fp));
	if(!io){return 0;}
	resize(r);
	io=int(fread(values,sizeof(Real),res,fp));
	if(io==res){return 1;}
	else{return 0;}
}
template<class Real>
int CircularArray<Real>::write(FILE* fp) const {
	int io;
	io=int(fwrite(&res,sizeof(int),1,fp));
	if(!io){return 0;}
	io=int(fwrite(values,sizeof(Real),res,fp));
	if(io==res){return 1;}
	else{return 0;}
}
template<class Real>
int CircularArray<Real>::resolution(void) const{return res;}
template<class Real>
int CircularArray<Real>::resize( int r ){
	if(r<0){return 0;}
	else{
		if(values){delete[] values;}
		values=NULL;
		res=0;
		if(r){
			values=new Real[r];
			if(!values){return 0;}
			else{res=r;}
		}
		clear();
		return 1;
	}
}
template<class Real>
void CircularArray<Real>::clear(void){if(res){memset(values,0,sizeof(Real)*res);}}

template<class Real>
Real& CircularArray<Real>::operator() (const int& i){
	int idx=i;
	if(idx<0){idx=res-((-idx)%res);}
	return values[idx%res];
}
template<class Real>
Real CircularArray<Real>::operator() (const int& i) const {
	int idx=i;
	if(idx<0){idx=res-((-idx)%res);}
	return values[idx%res];
}
template<class Real>
Real CircularArray<Real>::operator() (const double& x) const {
	int x1;
	double dx;
	if(x<0){x1=-int(-x)-1;}
	else{x1=int(x);}
	dx=x-x1;
	return (*this)(x1)*(Real(1.0)-dx)+(*this)(x1+1)*dx;
}
template<class Real>
Real CircularArray<Real>::operator() (const double& x){
	int x1;
	double dx;
	if(x<0){x1=-int(-x)-1;}
	else{x1=int(x);}
	dx=x-x1;
	return (*this)(x1)*(Real(1.0)-dx)+(*this)(x1+1)*dx;
}
template<class Real>
Real CircularArray<Real>::squareNorm(void) const{return Dot(*this,*this);}
template<class Real>
Real CircularArray<Real>::SquareDifference(const CircularArray& g1,const CircularArray& g2){return g1.squareNorm()+g2.squareNorm()-2*Dot(g1,g2);}
template<class Real>
Real CircularArray<Real>::Dot(const CircularArray& g1,const CircularArray& g2){
	Real d=0;
	if(g1.res != g2.res){
		fprintf(stderr,"Could not compare arrays of different sizes: %d != %d\n",g1.res,g2.res);
		exit(0);
	}
	for(int i=0;i<g1.res;i++) d+=g1.values[i]*g2.values[i];
	return Real(d/g1.res*2.0*PI);
}
