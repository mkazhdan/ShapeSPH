/***************************************************************************
  **************************************************************************
  
                           S2kit 1.0

          A lite version of Spherical Harmonic Transform Kit

   Peter Kostelec, Dan Rockmore
   {geelong,rockmore}@cs.dartmouth.edu
  
   Contact: Peter Kostelec
            geelong@cs.dartmouth.edu
  
   Copyright 2004 Peter Kostelec, Dan Rockmore

   This file is part of S2kit.

   S2kit is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   S2kit is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with S2kit; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

   See the accompanying LICENSE file for details.
  
  ************************************************************************
  ************************************************************************/


/*

Source code to test *inverse* spherical harmonic transform
using the seminaive and naive algorithms coded up during October, 1995.

In its current state, assuming that will seminaive at ALL orders.
If you wish to change this, modify the CUTOFF variable in this file, i.e.
cutoff = at what order to switch from semi-naive to naive algorithm

WILL PRECOMPUTE IN MEMORY EVERYTHING BEFORE DOING TRANSFORM!

Sample call

 test_s2_semi_memo_inv coeffsFile outputFile bw

 test_s2_semi_memo_inv coeff_bw8.dat samples_bw8.dat 8

The format of the input coefficient file will be an
interleaved real/imaginary parts of the coefficients, where
the coefficients are given in "code" order, as defined
in, e.g., test_s2_semi_memo:

          f(0,0) f(0,1) f(0,2) ... f(0,bw-1)
          f(1,1) f(1,2) ... f(1,bw-1)
          etc.
          f(bw-2,bw-2), f(bw-2,bw-1)
          f(bw-1,bw-1)
          f(-(bw-1),bw-1)
          f(-(bw-2),bw-2) f(-(bw-2),bw-1)
          etc.
          f(-2,2) ... f(-2,bw-1)
          f(-1,1) f(-1,2) ... f(-1,bw-1)


To help you out, the function

 seanindex(m, l, bw)

defined in FST_semi_memo.c, returns the array index of
the coefficient f_{l,m} (so you know where it is). Since we're
talking in C here, the indexing begins at 0.


The DEFAULT format of the output file of the forward transform test routine

 test_s2_semi_memo_for sampleFile outputFile bw [output_format]

is suitable for input for test_s2_semi_memo_inv.


The format of the output sample file will be an interleaved
real/imaginary parts of the function samples arranged in
"latitude-major" format, i.e. the function will be sampled
in this order:

  (theta_0, phi_0)
  (theta_0, phi_1)
  (theta_0, phi_2)
  ...
  (theta_0, phi_{bw-1})
  (theta_1, phi_0)  
  (theta_1, phi_1)
  ...
  (theta_{bw-1}, phi_{bw-1})

  where theta_k = pi*(2*j+1)/(4*bw)
        phi_j = 2*pi*k/(2*bw)

*/

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "fftw3.h"
#include "makeweights.h"
#include "cospmls.h"
#include "FST_semi_memo.h"
#include "csecond.h"

#define max(A, B) ((A) > (B) ? (A) : (B))

/**************************************************/
/**************************************************/



int main(int argc, char **argv)
{
  FILE *fp ;
  int i, bw, size, cutoff ;
  int rank, howmany_rank ;
  double *rcoeffs, *icoeffs, *rdata, *idata ;
  double *workspace, *weights ;
  double *seminaive_naive_tablespace, *trans_seminaive_naive_tablespace ;
  double **seminaive_naive_table, **trans_seminaive_naive_table;
  double tstart, tstop;
  fftw_plan idctPlan, ifftPlan ;
  fftw_iodim dims[1], howmany_dims[1];

  if (argc < 4)
    {
      fprintf(stdout,
	      "Usage: test_s2_semi_memo_inv coeffsFile outputFile bw\n");
      exit(0);
    }

  bw = atoi(argv[3]);
  size = 2*bw;

  /*** ASSUMING WILL SEMINAIVE ALL ORDERS ***/
  cutoff = bw ;

  /* allocate lots of memory */
  rcoeffs = (double *) malloc(sizeof(double) * (bw * bw));
  icoeffs = (double *) malloc(sizeof(double) * (bw * bw));
  rdata = (double *) malloc(sizeof(double) * (size * size));
  idata = (double *) malloc(sizeof(double) * (size * size));
  weights = (double *) malloc(sizeof(double) * 4 * bw);
  
  workspace = (double *) malloc(sizeof(double) * 
				((8 * (bw*bw)) + 
				 (10 * bw)));

  seminaive_naive_tablespace =
    (double *) malloc(sizeof(double) *
		      (Reduced_Naive_TableSize(bw,cutoff) +
		       Reduced_SpharmonicTableSize(bw,cutoff)));

  trans_seminaive_naive_tablespace =
    (double *) malloc(sizeof(double) *
		      (Reduced_Naive_TableSize(bw,cutoff) +
		       Reduced_SpharmonicTableSize(bw,cutoff)));

  /****
       At this point, check to see if all the memory has been
       allocated. If it has not, there's no point in going further.
  ****/

  if ( (rdata == NULL) || (idata == NULL) ||
       (rcoeffs == NULL) || (icoeffs == NULL) ||
       (weights == NULL) ||
       (seminaive_naive_tablespace == NULL) ||
       (trans_seminaive_naive_tablespace == NULL) ||
       (workspace == NULL) )
    {
      perror("Error in allocating memory");
      exit( 1 ) ;
    }

  /* now precompute the Legendres */
  fprintf(stdout,"Generating seminaive_naive tables...\n");
  seminaive_naive_table = SemiNaive_Naive_Pml_Table(bw, cutoff,
						    seminaive_naive_tablespace,
						    workspace);

  fprintf(stdout,"Generating trans_seminaive_naive tables...\n");
  trans_seminaive_naive_table =
    Transpose_SemiNaive_Naive_Pml_Table(seminaive_naive_table,
					bw, cutoff,
					trans_seminaive_naive_tablespace,
					workspace);

  /* construct fftw plans */

  /* make iDCT plan -> note that I will be using the GURU
     interface to execute this plan within the routine*/
      
  /* inverse DCT */
  idctPlan = fftw_plan_r2r_1d( 2*bw, weights, rdata,
			       FFTW_REDFT01, FFTW_ESTIMATE );

  /*
    now plan for inverse fft - note that this plans assumes
    that I'm working with a transposed array, e.g. the inputs
    for a length 2*bw transform are placed every 2*bw apart,
    the output will be consecutive entries in the array
  */
  rank = 1 ;
  dims[0].n = 2*bw ;
  dims[0].is = 2*bw ;
  dims[0].os = 1 ;
  howmany_rank = 1 ;
  howmany_dims[0].n = 2*bw ;
  howmany_dims[0].is = 1 ;
  howmany_dims[0].os = 2*bw ;

  /* inverse fft */
  ifftPlan = fftw_plan_guru_split_dft( rank, dims,
				       howmany_rank, howmany_dims,
				       rdata, idata,
				       workspace, workspace+(4*bw*bw),
				       FFTW_ESTIMATE );

  /* now make the weights */
  makeweights( bw, weights );

  /* now read in coefficients */
  fp = fopen(argv[1], "r");
  for( i = 0 ; i < bw*bw ; i++ )
    {
      /* first the real part of the coefficient */
      fscanf(fp, "%lf", rcoeffs + i );
      /* now the imaginary part */
      fscanf(fp, "%lf", icoeffs + i );
    }
  fclose( fp ) ;

  
  /* do the inverse spherical transform */
  tstart = csecond();
  InvFST_semi_memo(rcoeffs,icoeffs,
		   rdata, idata,
		   bw,
		   trans_seminaive_naive_table,
		   workspace,
		   0,
		   cutoff,
		   &idctPlan,
		   &ifftPlan );
  tstop = csecond();
  
  fprintf(stderr,"inv time \t = %.4e\n", tstop - tstart);
  fprintf(stdout,"about to write out samples\n");

  fp = fopen(argv[2], "w");
  for( i = 0 ; i < size*size ; i++ )
    fprintf(fp, "%.15f\n%.15f\n", 
	    rdata[i], idata[i]);
  fclose(fp);

  fprintf(stdout,"finished writing samples\n");
  
  /* now clean up */

  fftw_destroy_plan( ifftPlan );
  fftw_destroy_plan( idctPlan );

  free(trans_seminaive_naive_table);
  free(seminaive_naive_table);
  free(trans_seminaive_naive_tablespace);
  free(seminaive_naive_tablespace);
  free(workspace);
  free(weights);
  free(idata);
  free(rdata);
  free(icoeffs);
  free(rcoeffs);

  return 0 ;
  
}
