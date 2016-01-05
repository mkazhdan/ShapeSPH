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
  header file for full SO3 transform routines using fftw

  - the Wigner-ds are precomputed

  Forward_SO3_Naive_fftw_pc() - do forward full SO(3) transform
  Inverse_SO3_Naive_fftw_pc() - do inverse full SO(3) transform

  Both routines require fftw!!!


*/

#ifndef _SOFT_FFTW_PC_H
#define _SOFT_FFTW_PC_H 1

extern void Forward_SO3_Naive_fftw_pc( int ,
				       fftw_complex * ,
				       fftw_complex * ,
				       fftw_complex * ,
				       fftw_complex * ,
				       REAL * ,
				       fftw_plan * ,
				       fftw_plan * ,
				       REAL * ,
				       int ) ;

extern void Inverse_SO3_Naive_fftw_pc( int ,
				       fftw_complex * ,
				       fftw_complex * ,
				       fftw_complex * ,
				       fftw_complex * ,
				       REAL * ,
				       fftw_plan * ,
				       REAL * ,
				       int ) ;

#endif /* _SOFT_FFTW_PC_H */

