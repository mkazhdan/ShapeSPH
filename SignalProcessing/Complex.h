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
#ifndef COMPLEX_INCLUDED
#define COMPLEX_INCLUDED

// This templated class represents a complex number, supporting the standard
// arithmetic operations
template<class Real=float>
class Complex{
public:
	Real r,i;

	Complex(const Real& r=0,const Real& i=0);

	Complex operator + (const Real& r) const;
	Complex operator - (const Real& r) const;
	Complex operator * (const Real& r) const;
	Complex operator / (const Real& r) const;
	Complex& operator = (const Real& r);
	Complex& operator += (const Real& r);
	Complex& operator -= (const Real& r);
	Complex& operator *= (const Real& r);
	Complex& operator /= (const Real& r);

	
	Complex operator + (const Complex& c) const;
	Complex operator - (const Complex& c) const;
	Complex operator * (const Complex& c) const;
	Complex operator / (const Complex& c) const;
	Complex& operator = (const Complex& r);
	Complex& operator += (const Complex& c);
	Complex& operator -= (const Complex& c);
	Complex& operator *= (const Complex& c);
	Complex& operator /= (const Complex& c);

	Complex operator - (void) const;
	Complex conjugate(void) const;
	Real squareNorm(void) const;
};
#include "Complex.inl"
#endif // COMPLEX_INCLUDED
