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

Source code to test FORWARD spherical harmonic transform
using the seminaive and naive algorithms coded up during October, 1995.
It will compute the spherical coefficients of the function whose
sample values are given in the named input file, and will write
those coefficients in the named output file.

In its current state, assuming that will seminaive at ALL orders.
If you wish to change this, modify the CUTOFF variable in this file, i.e.
cutoff = at what order to switch from semi-naive to naive algorithm

WILL PRECOMPUTE IN MEMORY EVERYTHING BEFORE DOING TRANSFORM!

Sample call

 test_s2_semi_memo_for sampleFile outputFile bw [output_format]

 test_s2_semi_memo_for y31_bw8.dat y31_coef.dat 8

 output_format is an optional argument:
 = 0 -> coefficients in "code" order,
        i.e. suitable for test_s2_semi_memo_inv

 = 1 -> coefficients in prettier "human" order,
        i.e. if f_{l,m} is the coefficient of degree l, order m,
	then the coefficients will be arranged this way:
	f_{0,0},
	f_{1,-1}, f_{1,0}, f_{1,1},
	f_{2,-2}, f_{2,-1}, f_{2,0}, f_{2,1}, f_{2,2},
	...

The default output format is "code" order. To help you out,
the function

 seanindex(m, l, bw)

defined in FST_semi_memo.c, returns the array index of
the coefficient f_{l,m} (so you know where it goes). Since we're
talking in C here, the indexing begins at 0.

The format of the input sample file will be an interleaved
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

The format of the output depends on whether or not it's human
or code-ordered. For code-ordered, it's interleaved real/imaginary.
For human-ordered, it's verbose and it's one coefficient per line,
e.g. l = 2  m = 1  2.3 + 6 I


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

/**************************************************/
/**************************************************/



int main(int argc, char **argv)
{
  FILE *fp ;
  int i, bw, size ;
  int l, m, dummy;
  int cutoff, order ;
  int rank, howmany_rank ;
  double *rdata, *idata ;
  double *rcoeffs, *icoeffs ;
  double *weights ;
  double *seminaive_naive_tablespace, *workspace;
  double **seminaive_naive_table ;
  double tstart, tstop;
  fftw_plan dctPlan, fftPlan ;
  fftw_iodim dims[1], howmany_dims[1];

  if (argc < 4)
    {
      fprintf(stdout,
	      "Usage: test_s2_semi_memo_for sampleFile outputFile bw [output_format]\n");
      exit(0);
    }


  bw = atoi(argv[3]);

  /*** ASSUMING WILL SEMINAIVE ALL ORDERS ***/
  cutoff = bw ;
  size = 2*bw;

  /* allocate memory */
  rdata = (double *) malloc(sizeof(double) * (size * size));
  idata = (double *) malloc(sizeof(double) * (size * size));
  rcoeffs = (double *) malloc(sizeof(double) * (bw * bw));
  icoeffs = (double *) malloc(sizeof(double) * (bw * bw));
  weights = (double *) malloc(sizeof(double) * 4 * bw);
  seminaive_naive_tablespace =
    (double *) malloc(sizeof(double) *
		      (Reduced_Naive_TableSize(bw,cutoff) +
		       Reduced_SpharmonicTableSize(bw,cutoff)));
  workspace = (double *) malloc(sizeof(double) * 
				((8 * (bw*bw)) + 
				 (7 * bw)));


  /****
       At this point, check to see if all the memory has been
       allocated. If it has not, there's no point in going further.
  ****/

  if ( (rdata == NULL) || (idata == NULL) ||
       (rcoeffs == NULL) || (icoeffs == NULL) ||
       (seminaive_naive_tablespace == NULL) ||
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

  /* construct fftw plans */

  /* make DCT plan -> note that I will be using the GURU
     interface to execute these plans within the routines*/

  /* forward DCT */
  dctPlan = fftw_plan_r2r_1d( 2*bw, weights, rdata,
			      FFTW_REDFT10, FFTW_ESTIMATE ) ;
      
  /*
    fftw "preamble" ;
    note that this plan places the output in a transposed array
  */
  rank = 1 ;
  dims[0].n = 2*bw ;
  dims[0].is = 1 ;
  dims[0].os = 2*bw ;
  howmany_rank = 1 ;
  howmany_dims[0].n = 2*bw ;
  howmany_dims[0].is = 2*bw ;
  howmany_dims[0].os = 1 ;
  
  /* forward fft */
  fftPlan = fftw_plan_guru_split_dft( rank, dims,
				      howmany_rank, howmany_dims,
				      rdata, idata,
				      workspace, workspace+(4*bw*bw),
				      FFTW_ESTIMATE );


  /* now make the weights */
  makeweights( bw, weights );


  /* now read in samples */
  fp = fopen(argv[1],"r");
  for(i = 0 ; i < size*size ; i++ )
    {
      /* first the real part of the sample */
      fscanf(fp, "%lf", rdata + i );
      /* now the imaginary part */
      fscanf(fp, "%lf", idata + i );
    }
  fclose( fp ) ;

  
  /* now do the forward spherical transform */
  tstart = csecond();
  FST_semi_memo(rdata, idata,
		rcoeffs, icoeffs,
		bw,
		seminaive_naive_table,
		workspace,
		0,
		cutoff,
		&dctPlan,
		&fftPlan,
		weights );
  tstop = csecond();
    
  fprintf(stdout,"forward time \t = %.4e\n", tstop - tstart);

  fprintf(stdout,"about to write out coefficients\n");

  /* now write out coefficients, but in what format ? */
  if ( argc == 5 )
    order = atoi(argv[4]);
  else
    order = 0 ;

  fp = fopen( argv[2], "w" );
  if ( order == 0 )     /* code format */
    for( i = 0 ; i < bw*bw ; i ++ )
      fprintf(fp, "%.15f\n%.15f\n", 
	      rcoeffs[i], icoeffs[i]);
  else                  /* human format */
    for ( l = 0 ; l < bw ; l++ )
      for ( m = -l ; m < l + 1 ; m++ )
	{
	  dummy = seanindex(m, l, bw);
	  fprintf(fp, "l = %d\t m = %d\t %.15f + %.15f I\n",
		  l, m, rcoeffs[dummy], icoeffs[dummy]);
	}
  fclose(fp);

  fprintf(stdout,"finished writing coefficients\n");


  /* clean up */
  fftw_destroy_plan( fftPlan );
  fftw_destroy_plan( dctPlan );

  free(workspace);
  free(seminaive_naive_table);
  free(seminaive_naive_tablespace);
  free(weights);
  free(icoeffs);
  free(rcoeffs);
  free(idata);
  free(rdata);

  return 0 ;
  
}

