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

  header file for rotate.c -> functions having to do with
  rotating bandlimited functions defined on the sphere

*/


#ifndef _ROTATESO3_H
#define _ROTATESO3_H 1

/* #include "complex.h"  /* needed to define REAL */

extern void genExp( int ,
		    REAL ,
		    REAL * ,
		    REAL * ) ;

extern void wignerdmat( int ,
			REAL * ,
			REAL * ,
			REAL * ,
			REAL * ,
			REAL * ) ;

extern void rotateCoefDegree( int ,
			      REAL * , REAL * ,
			      REAL * , REAL * ,
			      REAL * , REAL * ,
			      REAL * , REAL * ,
			      int ,
			      REAL * ) ;

extern void rotateCoefAll( int ,
			   int ,
			   int ,
			   REAL ,
			   REAL ,
			   REAL ,
			   REAL * , REAL * ,
			   REAL * , REAL * ,
			   REAL * ) ;

extern void rotateFct( int , int , int ,
		       REAL * , REAL * ,
		       REAL * , REAL * ,
		       REAL , REAL , REAL ,
		       REAL * , 
		       REAL ** ,
		       REAL ** ) ;

#endif /* #ifndef _ROTATESO3_H */
