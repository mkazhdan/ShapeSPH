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
  routine for combining spherical (S^2) Fourier coefficients of
  two functions to get the SO(3) Fourier coefficients of their
  correlation.


  - will use Wigner-d symmetries and FFTW in the inverse SO(3) transform

*/

#include "fftw3.h"

#ifndef _SO3_CORRELATE_FFTW_H
#define _SO3_CORRELATE_FFTW_H 1

extern void so3CombineCoef_fftw( int ,
				 int ,
				 int ,
				 REAL * , REAL * ,
				 REAL * , REAL * ,
				 fftw_complex * ) ;

#endif
