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


/********************************************************************

  fft_grids_so3.c - routines to perform 1-d ffts on grids, expected
                    to be used in SO(3) transforms!!!


  The FFTcode-based routines in this file are:

  grid_fourier_so3
  grid_invfourier_so3


*/

#include <math.h>
#include "FFTcode.h"


/************************************************************************
  Computes the fourier transform of each row of the grid.  This
  is NOT the same as a 2-D Fourier transform.

  Used by FST_semi procedure

  Since this will be input to an associated legendre transform,
  the lines of longitude, or zones, are loaded into the rows
  in a transposed fashion for easy access by the Legendre
  transform rpocedure.  The grid is expected to
  be rows x cols, which is probably (2*bw)^2 x (2*bw)
  
  realgrid, imaggrid - (rows x cols) arrays of real and imag
                       input
  rmatrix, imatrix - (rows x cols) arrays of real and imag
                     output
  
  workspace - double pointer to array of (4 * cols)

  *********************************************************/


void grid_fourier_so3(double *realgrid,
		      double *imaggrid,
		      double *rmatrix,
		      double *imatrix,
		      int rows,
		      int cols,
		      double *scratchpad)
{
  int i;

  for (i = 0 ; i < rows ; i++) 
    {
      FFTInterp(realgrid+(i*cols), imaggrid+(i*cols),
		rmatrix+(i*cols),
		imatrix+(i*cols),
		cols, cols, scratchpad, 1);
    }

}

/***********************************************************************

  Same as above except for inverse Fourier transform is used
  used by InvFST_semi procedure 

  workspace = (24 * bw)

  **********************************************************************/

void grid_invfourier_so3(double *realgrid,
			 double *imaggrid,
			 double *rmatrix,
			 double *imatrix,
			 int rows,
			 int cols,
			 double *scratchpad)
{

  int i;

  for (i=0 ; i<rows ; i++)
    {
      FFTEval(realgrid+(i*cols), imaggrid+(i*cols),
	      rmatrix+(i*cols),
	      imatrix+(i*cols), cols, cols, scratchpad, 1);
    }
  
}
