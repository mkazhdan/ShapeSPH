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
template<class Real> Complex<Real>::Complex(const Real& r,const Real& i){this->r=r;this->i=i;}
template<class Real> Complex<Real> Complex<Real>::operator - (void) const{
	Complex<Real> c;
	c.r=-r;
	c.i=-i;
	return c;
}
template<class Real> Complex<Real> Complex<Real>::conjugate(void) const{
	Complex<Real> c;
	c.r=r;
	c.i=-i;
	return c;
}
template<class Real> Real Complex<Real>::squareNorm(void) const {return r*r+i*i;}
template<class Real> Complex<Real> Complex<Real>::operator + (const Real& r) const {
	Complex<Real> c;
	c.r=this->r+r;
	c.i=i;
	return c;
}
template<class Real> Complex<Real> Complex<Real>::operator - (const Real& r) const {
	Complex<Real> c;
	c.r=this->r-r;
	c.i=i;
	return c;
}
template<class Real> Complex<Real> Complex<Real>::operator * (const Real& r) const {
	Complex<Real> c;
	c.r=this->r*r;
	c.i=this->i*r;
	return c;
}
template<class Real> Complex<Real> Complex<Real>::operator / (const Real& r) const {
	Complex<Real> c;
	c.r=this->r/r;
	c.i=this->i/r;
	return c;
}
template<class Real> Complex<Real>& Complex<Real>::operator = (const Real& r){
	this->r=r;
	this->i=0;
	return *this;
}
template<class Real> Complex<Real>& Complex<Real>::operator += (const Real& r){
	this->r+=r;
	return *this;
}
template<class Real> Complex<Real>& Complex<Real>::operator -= (const Real& r){
	this->r-=r;
	return *this;
}
template<class Real> Complex<Real>& Complex<Real>::operator *= (const Real& r){
	this->r*=r;
	this->i*=r;
	return *this;
}
template<class Real> Complex<Real>& Complex<Real>::operator /= (const Real& r){
	this->r/=r;
	this->i/=r;
	return *this;
}
template<class Real> Complex<Real> Complex<Real>::operator + (const Complex<Real>& c) const {
	Complex<Real> out;
	out.r=r+c.r;
	out.i=i+c.i;
	return out;
}
template<class Real> Complex<Real> Complex<Real>::operator - (const Complex<Real>& c) const {
	Complex<Real> out;
	out.r=r-c.r;
	out.i=i-c.i;
	return out;
}
template<class Real> Complex<Real> Complex<Real>::operator * (const Complex<Real>& c) const {
	Complex<Real> out;
	out.r=r*c.r-i*c.i;
	out.i=r*c.i+i*c.r;
	return out;
}
template<class Real> Complex<Real> Complex<Real>::operator / (const Complex<Real>& c) const {
	Complex<Real> recip=c.conjugate()/c.squareNorm();
	return (*this)*recip;
}
template<class Real> Complex<Real>& Complex<Real>::operator = (const Complex<Real>& c){
	r=c.r;
	i=c.i;
	return *this;
}
template<class Real> Complex<Real>& Complex<Real>::operator += (const Complex<Real>& c){
	r+=c.r;
	i+=c.i;
	return *this;
}
template<class Real> Complex<Real>& Complex<Real>::operator -= (const Complex<Real>& c){
	r-=c.r;
	i-=c.i;
	return *this;
}
template<class Real> Complex<Real>& Complex<Real>::operator *= (const Complex<Real>& c){return ((*this)=((*this)*c));}
template<class Real> Complex<Real>& Complex<Real>::operator /= (const Complex<Real>& c){return ((*this)=((*this)/c));}
