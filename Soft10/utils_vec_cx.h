/***************************************************************************
  **************************************************************************
  
                SOFT: SO(3) Fourier transform code

                Version 1.0

  
   Peter Kostelec, Dan Rockmore
   {geelong,rockmore}@cs.dartmouth.edu
  
   Contact: Peter Kostelec
            geelong@cs.dartmouth.edu
  
  
   Copyright 2003 Peter Kostelec, Dan Rockmore
  
  
     This program is free software; you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation; either version 2 of the License, or
     (at your option) any later version.
  
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.
  
     You should have received a copy of the GNU General Public License
     along with this program; if not, write to the Free Software
     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
  
  
   Commercial use is absolutely prohibited.
  
   See the accompanying LICENSE file for details.
  
  ************************************************************************
  ************************************************************************/

/*

  header file for basic utility functions


  vec_add_cx() - add two COMPLEX vectors
  vec_mul_cx() - multiply COMPLEX vector by a scalar
  vec_pt_mul_cx() - pointwise-multiply REAL vector with a COMPLEX vector

  transpose_cx() - transpose a COMPLEX matrix

*/


#ifndef _UTILS_VEC_CX_H
#define _UTILS_VEC_CX_H 1

extern void vec_add_cx( fftw_complex * ,
			fftw_complex * ,
			fftw_complex * ,
			int ) ;

extern void vec_mul_cx( REAL ,
			fftw_complex * ,
			fftw_complex * ,
			int ) ;

extern void vec_pt_mul_cx( REAL * ,
			   fftw_complex * ,
			   fftw_complex * ,
			   int ) ;

extern void transpose_cx(  fftw_complex * ,
			   fftw_complex * ,
			   int ,
			   int ) ;

#endif /* _UTILS_VEC_CX_H */
