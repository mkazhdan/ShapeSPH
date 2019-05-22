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

  header file for the migration to interleaved real/imaginary
  arrays: defining the structure and some operations

*/


#ifndef _COMPLEX_FCTS_H
#define _COMPLEX_FCTS_H 1

/*
  heaven knows what would happen if I set REAL to
  be of type float -> too much code was written more
  or less ignoring this
*/
#define REAL double

struct COMPLEX
{
  REAL re ;
  REAL im ;
} ;


/*
  For the following, I am assuming that these complex
  numbers are in fftw format: double[2] where
  double[0] contains the real part and
  double[1] contains the imaginary part
*/

/* add two complex numbers */
#define CXADDRE( CX1, CX2 ) ( CX1[0] + CX2[0] )
#define CXADDIM( CX1, CX2 ) ( CX1[1] + CX2[1] )

/* subtract two complex numbers */
#define CXSUBRE( CX1, CX2 ) ( CX1[0] - CX2[0] )
#define CXSUBIM( CX1, CX2 ) ( CX1[1] - CX2[1] )

/* multiply two complex numbers */
#define CXMULRE( CX1, CX2 ) ( CX1[0] * CX2[0] - CX1[1] * CX2[1] )
#define CXMULIM( CX1, CX2 ) ( CX1[0] * CX2[1] + CX1[1] * CX2[0] )

/* multiply a complex number by a scalar */
#define CXSMULRE( SC, CX )  (CX[0] * SC)
#define CXSMULIM( SC, CX )  (CX[1] * SC)


#endif /* _COMPLEX_FCTS_H */

