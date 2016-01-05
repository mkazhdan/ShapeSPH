/***************************************************************************
  **************************************************************************
  
                Spherical Harmonic Transform Kit 2.6
  
   Sean Moore, Dennis Healy, Dan Rockmore, Peter Kostelec
   smoore@bbn.com, {healy,rockmore,geelong}@cs.dartmouth.edu
  
   Contact: Peter Kostelec
            geelong@cs.dartmouth.edu
  
  
   Copyright 1997-2003  Sean Moore, Dennis Healy,
                        Dan Rockmore, Peter Kostelec
  
  
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

  fft_grids.c - routines to perform 1-d ffts on grids, expected
                to be used in spherical transforms!!!


  Compilation flag information: 

  The FFTcode-based routines in this file are:

  grid_fourier
  grid_invfourier

*/

#include <math.h>
#include<string.h>
#include "primitive_FST.h"

#include "FFTcode.h"


/************************************************************************/
/************************************************************************/

/************************************************************************
  Computes the fourier transform of each row of the grid.  This
  is NOT the same as a 2-D Fourier transform.

  Used by FST_semi procedure

  Since this will be input to an associated legendre transform,
  the lines of longitude, or zones, are loaded into the rows
  in a transposed fashion for easy access by the Legendre
  transform rpocedure.  The grid is expected to
  be size * size, which is probably (2*bw) * (2*bw).
  
  realgrid, imaggrid - (size x size) arrays of real and imag
                       input
  rmatrix, imatrix - (size x size) arrays of real and imag
                     output
  size = 2 * bw
  
  workspace - double pointer to array of (6 * size) = (24 * bw)

  *********************************************************/


void grid_fourier(double *realgrid,
		  double *imaggrid,
		  double *rmatrix,
		  double *imatrix,
		  int size,
		  double *workspace)
{

  double *rout, *iout, *scratchpad;
  int i;

  /* need to assign workspace - need 6*size total */

  rout = workspace; /* needs size space */
  iout = rout + size; /* needs size space */
  scratchpad = iout + size; /* needs 4 * size space */

  if(1)
    {
      for (i=0; i<size; i++) 
	{
	  FFTInterp(realgrid+(i*size), imaggrid+(i*size),
		    rmatrix+(i*size),
		    imatrix+(i*size),
		    size, size, scratchpad, 1);
	}
    }
  else
    {
      for (i=0; i<size; i++) 
	{
	  FFTInterp(realgrid+(i*size), imaggrid+(i*size),
		    rout,
		    iout,
		    size, size, scratchpad, 1);

	  memcpy(rmatrix+(i*size),rout, sizeof(double) * size);
	  memcpy(imatrix+(i*size),iout, sizeof(double) * size);

	}
    }

  
  /* now transpose the results */
  transpose(rmatrix,size);
  transpose(imatrix,size);

}

/***********************************************************************

  Same as above except for inverse Fourier transform is used
  used by InvFST_semi procedure 

  workspace = (24 * bw)

  **********************************************************************/

void grid_invfourier(double *realgrid, double *imaggrid,
		     double *rmatrix, double *imatrix,
		     int size, double *workspace)
{

  double *rout, *iout, *scratchpad;
  int i;

  /* need to assign workspace - need 6*size total */
  
  rout = workspace; /* needs size space */
  iout = rout + size; /* needs size space */
  scratchpad = iout + size; /* needs 4 * size space */
  
  for (i=0; i<size; i++) {
    FFTEval(realgrid+(i*size), imaggrid+(i*size),
	    rmatrix+(i*size),
	    imatrix+(i*size), size, size, scratchpad, 1);
  }
  
}



/******************************************/
