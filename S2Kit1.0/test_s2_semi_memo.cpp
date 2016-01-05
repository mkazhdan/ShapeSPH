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

Source code to test full spherical harmonic transform
using the seminaive and naive algorithms coded up during October, 1995.

In its current state, assuming that will seminaive at ALL orders.
If you wish to change this, modify the CUTOFF variable in this file, i.e.
cutoff = at what order to switch from semi-naive to naive algorithm

WILL PRECOMPUTE IN MEMORY EVERYTHING BEFORE DOING TRANSFORMS!

Idea is to record the execution time for the spectral->grid->spectral
roundtrip, after precomputations have been done.

The strategy is to generate random coefficients representing
a complete spherical harmonic series of funcions Y(m,l).
The ordering of the these coefficients is assumed to be the
same as that of the output of the FST_seminaive() routine, namely
  

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

This means that there are (bw*bw) coefficients.

Once the coefficients are generated, the corresponding
function is synthesized using InvFST_semi_memo(), then
transformed (analyzed) using FST_semi_memo(). Timing data
is printed.

Sample call

% test_s2_semi_memo bw loops [error_file]


Appropriate timing data will be printed out.


NOTE: In this program, the coefficients generated are such that the grid points
(sample values) produced are REAL. The routines InvFST_semi_fly() and
FST_semi_fly() will take advantage of this. If you wish to change this,
change the 7th argument of each function further down from "1" to "0".
This is also documented in the file FST_semi_fly.c .

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
  FILE *errorsfp ;
  int bw, size, loops, cutoff;
  int i, j;
  int l, m, dummy;
  double *rcoeffs, *icoeffs, *rdata, *idata, *rresult, *iresult;
  double *seminaive_naive_tablespace, *trans_seminaive_naive_tablespace,
    *workspace;
  double **seminaive_naive_table, **trans_seminaive_naive_table;
  double dumx, dumy, fudge;
  double *relerror, *curmax, granderror, grandrelerror;
  double realtmp, imagtmp,origmag, tmpmag;
  double total_time, for_time, inv_time;
  double tstart, tstop;
  double ave_error, ave_relerror, stddev_error, stddev_relerror;
  time_t seed;
  fftw_plan dctPlan, idctPlan ;
  fftw_plan fftPlan, ifftPlan ;
  double *weights ;
  fftw_iodim dims[1], howmany_dims[1];
  int rank, howmany_rank ;

  granderror = 0.0;
  grandrelerror = 0.0; 

  if (argc < 3)
    {
      fprintf(stdout,"Usage: test_s2_semi_memo bw loops [error_file]\n");
      exit(0);
    }


  bw = atoi(argv[1]);

  /*** ASSUMING WILL SEMINAIVE ALL ORDERS ***/
  cutoff = bw ;
  
  size = 2*bw;
  loops = atoi(argv[2]);

  total_time = 0.0;
  for_time = 0.0;
  inv_time = 0.0;


  /* allocate lots of memory */
  rcoeffs = (double *) malloc(sizeof(double) * (bw * bw));
  icoeffs = (double *) malloc(sizeof(double) * (bw * bw));
  rdata = (double *) malloc(sizeof(double) * (size * size));
  idata = (double *) malloc(sizeof(double) * (size * size));
  rresult = (double *) malloc(sizeof(double) * (bw * bw));
  iresult = (double *) malloc(sizeof(double) * (bw * bw));

  seminaive_naive_tablespace =
    (double *) malloc(sizeof(double) *
		      (Reduced_Naive_TableSize(bw,cutoff) +
		       Reduced_SpharmonicTableSize(bw,cutoff)));

  trans_seminaive_naive_tablespace =
    (double *) malloc(sizeof(double) *
		      (Reduced_Naive_TableSize(bw,cutoff) +
		       Reduced_SpharmonicTableSize(bw,cutoff)));
  
  workspace = (double *) malloc(sizeof(double) * 
				((8 * (bw*bw)) + 
				 (16 * bw)));

  /** space for errors **/
  relerror = (double *) malloc(sizeof(double) * loops);
  curmax = (double *) malloc(sizeof(double) * loops);

  /****
       At this point, check to see if all the memory has been
       allocated. If it has not, there's no point in going further.
  ****/

  if ( (rdata == NULL) || (idata == NULL) ||
       (rresult == NULL) || (iresult == NULL) ||
       (rcoeffs == NULL) || (icoeffs == NULL) ||
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
  
  /*** generate a seed, needed to generate random data ***/
  time(&seed);
  srand48( seed );


  /* make array for weights, and construct fftw plans */
  weights = (double *) malloc(sizeof(double) * 4 * bw);

  /* make DCT plans -> note that I will be using the GURU
     interface to execute these plans within the routines*/

  /* forward DCT */
  dctPlan = fftw_plan_r2r_1d( 2*bw, weights, rdata,
			      FFTW_REDFT10, FFTW_ESTIMATE ) ;
      
  /* inverse DCT */
  idctPlan = fftw_plan_r2r_1d( 2*bw, weights, rdata,
			       FFTW_REDFT01, FFTW_ESTIMATE );

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


  fprintf(stdout,"about to enter loop\n\n");
  for(i=0; i<loops; i++){

    /**** loop to generate spherical harmonic coefficients *****/
    for(m=0;m<bw;m++)
      for(l=m;l<bw;l++){
	dumx = 2.0 * (drand48()-0.5) ;
	dumy = 2.0 * (drand48()-0.5) ;
	dummy = seanindex(m,l,bw);
	rcoeffs[dummy] = dumx;
	icoeffs[dummy] = dumy;
	dummy = seanindex(-m,l,bw);
	rcoeffs[dummy] = ((double) pow(-1.0, (double) m)) * dumx;
	icoeffs[dummy] = ((double) pow(-1.0, (double) (m + 1))) * dumy;
      }

    /* have to zero out the m=0 coefficients, since those are real */
    for(m=0;m<bw;m++)
      icoeffs[m] = 0.0;

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
    inv_time += (tstop - tstart);
    
    fprintf(stdout,"inv time \t = %.4e\n", tstop - tstart);

    /* now do the forward spherical transform */
    tstart = csecond();

    FST_semi_memo(rdata, idata,
		  rresult, iresult,
		  bw,
		  seminaive_naive_table,
		  workspace,
		  0,
		  cutoff,
		  &dctPlan,
		  &fftPlan,
		  weights );

    tstop = csecond();
    for_time += (tstop - tstart);
    
    fprintf(stdout,"forward time \t = %.4e\n", tstop - tstart);


    /* now to compute the error */
    relerror[i] = 0.0;
    curmax[i] = 0.0;
    for(j=0;j<(bw*bw);j++){
      realtmp = rresult[j]-rcoeffs[j];
      imagtmp = iresult[j]-icoeffs[j];
      origmag = sqrt((rcoeffs[j]*rcoeffs[j]) + (icoeffs[j]*icoeffs[j]));
      tmpmag  = sqrt((realtmp*realtmp) + (imagtmp*imagtmp));
      relerror[i] = max(relerror[i],tmpmag/(origmag + pow(10.0, -50.0)));
      curmax[i]  = max(curmax[i],tmpmag);
    }
    
    fprintf(stdout,"r-o error\t = %.12f\n", curmax[i]);
    fprintf(stdout,"(r-o)/o error\t = %.12f\n\n", relerror[i]);
    
    granderror += curmax[i];
    grandrelerror += relerror[i];
        
  }
  
  total_time = inv_time + for_time;

  ave_error = granderror / ( (double) loops );
  ave_relerror = grandrelerror / ( (double) loops );
  stddev_error = 0.0 ; stddev_relerror = 0.0;
  for( i = 0 ; i < loops ; i ++ )
    {
      stddev_error += pow( ave_error - curmax[i] , 2.0 );
      stddev_relerror += pow( ave_relerror - relerror[i] , 2.0 );
    }
  /*** this won't work if loops == 1 ***/
  if( loops != 1 )
    {
      stddev_error = sqrt(stddev_error / ( (double) (loops - 1) ) );
      stddev_relerror = sqrt(stddev_relerror / ( (double) (loops - 1) ) );
    }


  fprintf(stdout,"Program: test_s2_semi_memo\n");
  fprintf(stdout,"Bandwidth = %d\n", bw);

#ifndef WALLCLOCK
  fprintf(stdout,"Total elapsed cpu time :\t\t %.4e seconds.\n",
	  total_time);
  fprintf(stdout,"Average cpu forward per iteration:\t %.4e seconds.\n",
	  for_time/((double) loops));  
  fprintf(stdout,"Average cpu inverse per iteration:\t %.4e seconds.\n",
	  inv_time/((double) loops));
#else
  fprintf(stdout,"Total elapsed wall time :\t\t %.4e seconds.\n",
	  total_time);
  fprintf(stdout,"Average wall forward per iteration:\t %.4e seconds.\n",
	  for_time/((double) loops));  
  fprintf(stdout,"Average wall inverse per iteration:\t %.4e seconds.\n",
	  inv_time/((double) loops));
#endif

  fprintf(stdout,"Average r-o error:\t\t %.4e\t",
	  granderror/((double) loops));
  fprintf(stdout,"std dev: %.4e\n",stddev_error);
  fprintf(stdout,"Average (r-o)/o error:\t\t %.4e\t",
	  grandrelerror/((double) loops));
  fprintf(stdout,"std dev: %.4e\n\n",stddev_relerror);


  if (argc == 4)
    {
      fudge = -1.0;
      errorsfp = fopen(argv[3],"w");
      for(m = 0 ; m < bw ; m++ )
	{
	  fudge *= -1.0;
	  for(l = m ; l< bw ; l++ )
	    {
	      dummy = seanindex(m,l,bw);
	      fprintf(errorsfp,
		      "dummy = %d\t m = %d\tl = %d\t%.10f  %.10f\n",
		      dummy, m, l,
		      fabs(rcoeffs[dummy] - rresult[dummy]),
		      fabs(icoeffs[dummy] - iresult[dummy]));
	      
	      dummy = seanindex(-m,l,bw);	      
	      fprintf(errorsfp,
		      "dummy = %d\t m = %d\tl = %d\t%.10f  %.10f\n",
		      dummy, -m, l,
		      fabs(rcoeffs[dummy] - rresult[dummy]),
		      fabs(icoeffs[dummy] - iresult[dummy]));
	      
	    }
	}
      
      fclose(errorsfp);
      
    }

  fftw_destroy_plan( ifftPlan );
  fftw_destroy_plan( fftPlan );
  fftw_destroy_plan( idctPlan );
  fftw_destroy_plan( dctPlan );


  free( weights );

  free(trans_seminaive_naive_table); free(seminaive_naive_table);
  free(curmax); free(relerror);  free(workspace);
  free(trans_seminaive_naive_tablespace);
  free(seminaive_naive_tablespace);
  free(iresult); free(rresult);
  

  free(idata);

  
  free(rdata);
  free(icoeffs); free(rcoeffs);

  return 0 ;
  
}


