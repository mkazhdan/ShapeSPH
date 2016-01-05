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

template< class Real >
RotationGrid< Real >::RotationGrid(void){;}
template< class Real >
RotationGrid< Real >::RotationGrid( int res ){ resize(res); }
template< class Real >
RotationGrid<Real>::~RotationGrid(void){ ; }
template<class Real>
void RotationGrid<Real>::SetCoordinates(const Real& theta,const Real& phi,const Real& psi,Real matrix[3][3]){
	Real ci, cj, ch, si, sj, sh, cc, cs, sc, ss;
    ci = cos(theta);	cj = cos(phi);	ch = cos(psi);
    si = sin(theta);	sj = sin(phi);	sh = sin(psi);

    cc = ci*ch;		cs = ci*sh;		sc = si*ch;		ss = si*sh;

	matrix[1][1] =  cj;		matrix[1][2] =  sj*si;		matrix[1][0] =  sj*ci;
	matrix[2][1] =  sj*sh;	matrix[2][2] = -cj*ss+cc;	matrix[2][0] = -cj*cs-sc;
	matrix[0][1] = -sj*ch;	matrix[0][2] =  cj*sc+cs;	matrix[0][0] =  cj*cc-ss;
}
template<class Real>
void RotationGrid<Real>::SetCoordinates(const Real matrix[3][3],Real& theta,Real& phi,Real& psi){
	Real sy;
	sy=Real(sqrt(matrix[1][2]*matrix[1][2]+matrix[1][0]*matrix[1][0]));
	if(sy>16*FLT_EPSILON){
		theta=Real( atan2(matrix[1][2],matrix[1][0]));
		phi=Real( atan2(sy,matrix[1][1]));
		psi=Real( atan2(matrix[2][1],-matrix[0][1]));
	}
	else{
		theta=Real( atan2(-matrix[2][0],matrix[2][2]));
		phi=Real(atan2(sy,matrix[1][1]));
		psi=0;
	}
}

template<class Real>
void RotationGrid<Real>::SetCoordinates(const Real axis[3],const Real& angle,Real matrix[3][3]){
	Real m1[3][3],m2[3][3];
	Real theta,phi;
	SphericalGrid<>::SetCoordinates(axis,theta,phi);
	SetCoordinates(0,-phi,-theta,m1);
	SetCoordinates(theta,phi,angle,m2);

	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			matrix[i][j]=0;
			for(int k=0;k<3;k++){
				matrix[i][j]+=m2[k][j]*m1[i][k];
			}
		}
	}
}

template<class Real>
void RotationGrid<Real>::setCoordinates(const int& i,const int& j,const int& k,Real matrix[3][3]) const {
	Real theta,phi,psi;
	theta=Real(2.0*i*PI/res);
	phi=Real(PI*(2.0*j+1)/(2.0*res));
	psi=Real(2.0*k*PI/res);

	SetCoordinates(theta,phi,psi,matrix);
}
template<class Real>
void RotationGrid<Real>::setCoordinates(const Real& i,const Real& j,const Real& k,Real matrix[3][3]) const {
	Real theta,phi,psi;
	theta=Real(2.0*i*PI/res);
	phi=Real(PI*(2.0*j+1)/(2.0*res));
	psi=Real(2.0*k*PI/res);

	SetCoordinates(theta,phi,psi,matrix);
}
template<class Real>
void RotationGrid<Real>::setCoordinates(const Real matrix[3][3],Real& i,Real& j,Real& k) const{
	Real theta,phi,psi;
	SetCoordinates(matrix,theta,phi,psi);
	i=(theta*res)/Real(2.0*PI);
	j=Real(2.0*res*phi-PI)/Real(2*PI);
	k=(psi*res)/Real(2.0*PI);
}
template<class Real>
Real& RotationGrid<Real>::operator() (const int& i,const int& j,const int& k){
	int x=i+res/2;
	int y=j;
	int z=k+res/2;

	if(y<0){y=2*res-((-y)%(2*res));}
	y%=2*res;
	if(y>=res){
		y=2*res-y-1;
		x+=res/2;
		z+=res/2;
	}
	if(x<0){x=res-((-x)%res);}
	x=x%res;
	if(z<0){z=res-((-z)%res);}
	z=z%res;
	return values[x*res*res+y*res+z];
}
template<class Real>
Real RotationGrid<Real>::operator() (const int& i,const int& j,const int& k) const {
	int x=i+res/2;
	int y=j;
	int z=k+res/2;

	if(y<0){y=2*res-((-y)%(2*res));}
	y%=2*res;
	if(y>=res){
		y=2*res-y-1;
		x+=res/2;
		z+=res/2;
	}
	if(x<0){x=res-((-x)%res);}
	x=x%res;
	if(z<0){z=res-((-z)%res);}
	z=z%res;
	return values[x*res*res+y*res+z];
}
#if 0
template<class Real>
Real RotationGrid<Real>::operator() (const double& x,const double& y,const double& z){
	int x1,y1,z1;
	double dx,dy,dz;

	if(x<0){x1=-int(-x)-1;}
	else{x1=int(x);}
	if(y<0){y1=-int(-y)-1;}
	else{y1=int(y);}
	if(z<0){z1=-int(-z)-1;}
	else{z1=int(z);}

	dx=x-x1;
	dy=y-y1;
	dz=z-z1;
	return 
		(*this)(x1  ,y1  ,z1  )*(Real(1.0)-dx)*(Real(1.0)-dy)*(Real(1.0)-dz)+
		(*this)(x1+1,y1  ,z1  )*(          dx)*(Real(1.0)-dy)*(Real(1.0)-dz)+
		(*this)(x1  ,y1+1,z1  )*(Real(1.0)-dx)*(          dy)*(Real(1.0)-dz)+
		(*this)(x1+1,y1+1,z1  )*(          dx)*(          dy)*(Real(1.0)-dz)+
		(*this)(x1  ,y1  ,z1+1)*(Real(1.0)-dx)*(Real(1.0)-dy)*(          dz)+
		(*this)(x1+1,y1  ,z1+1)*(          dx)*(Real(1.0)-dy)*(          dz)+
		(*this)(x1  ,y1+1,z1+1)*(Real(1.0)-dx)*(          dy)*(          dz)+
		(*this)(x1+1,y1+1,z1+1)*(          dx)*(          dy)*(          dz);
}
#endif
template<class Real>
Real RotationGrid<Real>::sample( Real x , Real y , Real z ) const {
	int x1,y1,z1;
	Real dx,dy,dz;

	if(x<0){x1=-int(-x)-1;}
	else{x1=int(x);}
	if(y<0){y1=-int(-y)-1;}
	else{y1=int(y);}
	if(z<0){z1=-int(-z)-1;}
	else{z1=int(z);}

	dx=x-x1;
	dy=y-y1;
	dz=z-z1;
	return 
		(*this)(x1  ,y1  ,z1  )*(Real(1.0)-dx)*(Real(1.0)-dy)*(Real(1.0)-dz)+
		(*this)(x1+1,y1  ,z1  )*(          dx)*(Real(1.0)-dy)*(Real(1.0)-dz)+
		(*this)(x1  ,y1+1,z1  )*(Real(1.0)-dx)*(          dy)*(Real(1.0)-dz)+
		(*this)(x1+1,y1+1,z1  )*(          dx)*(          dy)*(Real(1.0)-dz)+
		(*this)(x1  ,y1  ,z1+1)*(Real(1.0)-dx)*(Real(1.0)-dy)*(          dz)+
		(*this)(x1+1,y1  ,z1+1)*(          dx)*(Real(1.0)-dy)*(          dz)+
		(*this)(x1  ,y1+1,z1+1)*(Real(1.0)-dx)*(          dy)*(          dz)+
		(*this)(x1+1,y1+1,z1+1)*(          dx)*(          dy)*(          dz);
}
template<class Real>
Real RotationGrid<Real>::squareNorm(void) const{return Dot(*this,*this);}
template<class Real>
Real RotationGrid<Real>::SquareDifference(const RotationGrid& g1,const RotationGrid& g2){return g1.squareNorm()+g2.squareNorm()-2*Dot(g1,g2);}
template<class Real>
Real RotationGrid<Real>::Dot(const RotationGrid& g1,const RotationGrid& g2){
	double d=0,c=1;
	if(g1.res != g2.res){
		fprintf(stderr,"Could not compare arrays of different sizes: %d != %d\n",g1.res,g2.res);
		exit(0);
	}
	for(int i=0;i<g1.res;i++){
		double t1 =cos(PI*(2.0*i+2)/(2.0*g1.res));
		double t2=((c-t1)/(2.0*g1.res)/(2.0*g1.res));
		c=t1;
		for(int j=0;j<g1.res;j++){for(int k=0;k<g1.res;k++){d+=g1.values[j*g1.res*g1.res+i*g1.res+k]*g2.values[j*g1.res*g1.res+i*g1.res+k]*t2;}}
	}
	return Real(d*4.0/3.0*PI);
}
