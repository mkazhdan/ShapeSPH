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

#include <FST_semi_memo_fftw.h>
#include <cospmls.h>
#include "fftw3.h"
#include <math.h>

//////////////////
// FourierKeyS2 //
//////////////////
template< class Real > FourierKeyS2<Real>::FourierKeyS2( void ) : bw(0) , values(NULL) { ; }
template< class Real > FourierKeyS2< Real >::FourierKeyS2( int resolution ) : bw(0) , values(NULL) { resize( resolution ); }
template<class Real> FourierKeyS2<Real>::FourierKeyS2( const FourierKeyS2< Real >& key ) : bw(0) , values(NULL)
{
	resize( key.resolution() );
	memcpy( values , key.values , sizeof( Complex< Real > ) * Entries(bw) );
}
template< class Real >
FourierKeyS2< Real >& FourierKeyS2< Real >::operator = ( const FourierKeyS2< Real >& key )
{
	resize( key.resolution() );
	memcpy( values , key.values , sizeof( Complex< Real > ) * Entries(bw) );
	return *this;
}

template<class Real> FourierKeyS2<Real>::~FourierKeyS2(void){
	if(values){delete[] values;}
	values=NULL;
	bw=0;
}
template<class Real> int FourierKeyS2<Real>::read(const char* fileName){
	FILE* fp=fopen(fileName,"rb");
	if(!fp){return 0;}
	int r=read(fp);
	fclose(fp);
	return r;
}
template<class Real> int FourierKeyS2<Real>::write(const char* fileName) const{
	FILE* fp=fopen(fileName,"wb");
	if(!fp){return 0;}
	int w=write(fp);
	fclose(fp);
	return w;
}
template<class Real> int FourierKeyS2<Real>::read(FILE* fp){
	int b,r;
	r=int(fread(&b,sizeof(int),1,fp));
	if(!r){return 0;}
	resize(b);
	r=int(fread(values,sizeof(Complex<Real>),((bw*bw+bw)>>1),fp));
	if(r==((bw*bw+bw)>>1)){return 1;}
	else{return 0;}
}
template<class Real> int FourierKeyS2<Real>::write(FILE* fp) const {
	int w;
	w=int(fwrite(&bw,sizeof(int),1,fp));
	if(!w){return 0;}
	w=int(fwrite(values,sizeof(Complex<Real>),((bw*bw+bw)>>1),fp));
	if(w==((bw*bw+bw)>>1)){return 1;}
	else{return 0;}
}
template<class Real> int FourierKeyS2<Real>::bandWidth( void ) const{return bw;}
template<class Real> int FourierKeyS2<Real>::resolution( void ) const {return bw*2;}
template<class Real> int FourierKeyS2<Real>::resize( int resolution , bool clr )
{
	int b=resolution>>1;
	if( b<0 ) return 0;
	else if( b!=bw )
	{
		if( values ) delete[] values;
		values = NULL;
		bw = 0;
		if( b )
		{
			values = new Complex<Real>[(b*b+b)>>1];
			if( !values ) return 0;
			else bw=b;
		}
	}
	if(clr) clear();
	return 1;
}
template<class Real> Complex<Real>& FourierKeyS2<Real>::operator() ( int b , int i)       { return values[(b-i)+(i*bw)-(i*i-i)/2];}
template<class Real> Complex<Real>  FourierKeyS2<Real>::operator() ( int b , int i) const { return values[(b-i)+(i*bw)-(i*i-i)/2];}
template<class Real> void FourierKeyS2<Real>::clear( void ){ if(bw) memset(values,0,sizeof(Complex<Real>)*((bw*bw+bw)>>1));}
template< class Real >
void FourierKeyS2< Real >::Add( const FourierKeyS2< Real >& key )
{
	for( int b=0 ; b<bw && b<=key.bw ; b++ ) for( int i=0 ; i<=b ; i++ ) (*this)(b,i) += key(b,i);
}
template< class Real >
void FourierKeyS2< Real >::Scale( Real s )
{
	for( int b=0 ; b<bw ; b++ ) for( int i=0 ; i<=b ; i++ ) (*this)(b,i) *= s;
}
template< class Real >
Real FourierKeyS2< Real >::InnerProduct( const FourierKeyS2< Real >& key ) const
{
	Real dot = 0;
	// For real data, we have
	// f-hat(l,-m) = (-1)^m * conjugate(f-hat(l,m))
	// Furthermore, we will assume that f-hat(l,0) is strictly real

	// The zonal harmonics
	for( int b=0 ; b<bw && b<key.bw ; b++ )
	{
		dot += (*this)(b,0).r * key(b,0).r;
		for( int i=1 ; i<=b ; i++ ) dot += ( (*this)(b,i) * key(b,i).conjugate() ).r * 2;
	}
	return dot;
}
template<class Real> int FourierKeyS2<Real>::Entries( int bw ){return (bw*bw+bw)>>1;}
/////////////////////////////////////
// HarmonicTransform::ScratchSpace //
/////////////////////////////////////
template<class Real>
HarmonicTransform<Real>::ScratchSpace::ScratchSpace( void )
{
	bw=0;
	workSpace=resultSpace=transposeResultSpace=NULL;
	table=transposeTable=NULL;
#if NEW_HARMONIC
	weights=NULL;
#endif // NEW_HARMONIC
}
template<class Real>
HarmonicTransform<Real>::ScratchSpace::~ScratchSpace(void){resize(0);}
template<class Real>
void HarmonicTransform<Real>::ScratchSpace::resize( const int& b )
{
	if( b!=bw )
	{
		int size=b*2;
		if(workSpace)				{delete[] workSpace;}
		if(resultSpace)				{delete[] resultSpace;}
		if(transposeResultSpace)	{delete[] transposeResultSpace;}
		if(table)					{delete[] table;}
		if(transposeTable)			{delete[] transposeTable;}
#if NEW_HARMONIC
		if( weights ) delete[] weights;
#endif // NEW_HARMONIC
		bw=0;
		workSpace=resultSpace=transposeResultSpace=NULL;
		table=transposeTable=NULL;
#if NEW_HARMONIC
		weights = NULL;
#endif // NEW_HARMONIC
		if( b>0 )
		{
			bw = b;
			workSpace = new Real[4*bw*bw+36*bw];
			resultSpace = new Real[Spharmonic_TableSize(bw)];
			transposeResultSpace = new Real[Spharmonic_TableSize(bw)];
#if NEW_HARMONIC
			weights = new double*[4*bw];
#endif // NEW_HARMONIC
			table			=           Spharmonic_Pml_Table(bw,resultSpace,workSpace);
			transposeTable	= Transpose_Spharmonic_Pml_Table(table,bw,transposeResultSpace,workSpace);
		}
	}
}
///////////////////////
// HarmonicTransform //
///////////////////////
template< class Real > HarmonicTransform< Real >::HarmonicTransform( void ){ ; }
template< class Real > HarmonicTransform< Real >::HarmonicTransform( int resolution ){ resize( resolution ); }
template<class Real>
void HarmonicTransform<Real>::resize( const int& resolution ){ scratch.resize(resolution>>1); }
template<>
int HarmonicTransform< double >::ForwardFourier( SphericalGrid< double >& g , FourierKeyS2< double >& key )
{
	int sz,bw;
	sz=g.resolution();
	bw=sz>>1;
	if(key.resolution()!=sz){key.resize(sz);}
	scratch.resize(bw);
	FST_semi_memo_fftw( g[0] , (fftw_complex*)&key(0,0) , sz , scratch.table , scratch.workSpace );
	return 1;
}
template<>
int HarmonicTransform<float>::ForwardFourier( SphericalGrid<float>& g , FourierKeyS2<float>& key )
{
	int sz = g.resolution() , bw = sz>>1;
	if( key.resolution()!=sz ) key.resize(sz);
	scratch.resize( bw );
	FST_semi_memo_fftw( g[0] , (fftwf_complex*)&key(0,0) , sz , scratch.table , scratch.workSpace );
	return 1;
}
template<class Real>
int HarmonicTransform<Real>::ForwardFourier(SphericalGrid<Real>&,FourierKeyS2<Real>&){
	fprintf(stderr,"Harmonic Transform only supported for floats and doubles\n");
	return 0;
}
template<>
int HarmonicTransform<double>::InverseFourier(FourierKeyS2<double>& key,SphericalGrid<double>& g){
	if(key.resolution()!=g.resolution()){g.resize(key.resolution());}
	int bw=key.bandWidth(),sz=g.resolution();
	scratch.resize(bw);

	InvFST_semi_memo_fftw((fftw_complex*)&key(0,0),g[0],sz,scratch.transposeTable,scratch.workSpace);
	return 1;
}
template<>
int HarmonicTransform<float>::InverseFourier(FourierKeyS2<float>& key,SphericalGrid<float>& g)
{
	if(key.resolution()!=g.resolution()){g.resize(key.resolution());}
	int bw=key.bandWidth(),sz=g.resolution();
	scratch.resize(bw);

	InvFST_semi_memo_fftw((fftwf_complex*)&key(0,0),g[0],sz,scratch.transposeTable,scratch.workSpace);
	return 1;
}
template<class Real>
int HarmonicTransform<Real>::InverseFourier(FourierKeyS2<Real>&,SphericalGrid<Real>&){
	fprintf(stderr,"Harmonic Transform only supported for floats and doubles\n");
	return 0;
}