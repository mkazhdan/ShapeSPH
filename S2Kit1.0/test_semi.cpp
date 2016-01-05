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
  test_semi.c - top-level code for computing Legendre transform using the
  seminaive algorithm; will do forward and inverse DLTs, and return
  error and timing results

  m    - order of the problem

  bw   - bandwidth

  loops - number of loops thru timed portion of code.  Intended
          to reduce noise due to multiprocessing and 
          discretization errors
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>   /* to declare memcpy */
#include <time.h>

#include "fftw3.h"
#include "makeweights.h"
#include "cospmls.h"
#include "primitive.h"
#include "seminaive.h"
#include "csecond.h"

#define mymax(a, b) ((a) > (b) ? (a) : (b))

/**************************************************/
/**************************************************/

int main( int argc, char **argv)
{
  int i, j;
  int n, bw, m, k, loops;
  double *samples, *coeffs, *newcoeffs;
  double *workspace;
  double *cos_pml_table, *transpose_cos_pml_table ;
  double *eval_args;
  double *sin_values ;
  double *curmax;
  double tmp_error, ave_error, sum_error;
  double tmp_relerror, ave_relerror, sum_relerror;
  double stddev_error, stddev_relerror;
  double *error, *relerror;
  double tstartf, tstopf, totaltimef ;
  double tstarti, tstopi, totaltimei ;
  double *weights ;
  long int seed ;
  fftw_plan forwardDCT, inverseDCT ;

  if (argc < 4)
    {
      fprintf(stdout,"Usage: test_semi m bw loops\n");
      exit(0);
    }

  m = atoi(argv[1]);
  bw = atoi(argv[2]);
  loops = atoi(argv[3]);

  n = 2 * bw;

  /* space for samples and coefficients */
  samples = (double *) malloc(sizeof(double) * n);  
  coeffs = (double *) malloc(sizeof(double) * (bw - m));
  newcoeffs = (double *) malloc(sizeof(double) * (bw - m));

  /* space for precomputed Legendres */
  cos_pml_table =
    (double *) malloc(sizeof(double) * TableSize(m,bw)) ;
  transpose_cos_pml_table =
    (double *) malloc(sizeof(double) * TableSize(m,bw));

  /* need sin values if order of transform is odd */
  eval_args = (double *) malloc(sizeof(double) * n );
  sin_values = (double *) malloc(sizeof(double) * n );

  /* space for error stuff */
  error = (double *) malloc(sizeof(double) * loops);
  relerror = (double *) malloc(sizeof(double) * loops);
  curmax = (double *) malloc(sizeof(double) * loops);

  /* temp space */
  workspace = (double *) malloc(sizeof(double) * 9 * bw);

  /* for weights */
  weights = (double *) malloc(sizeof(double) * 4 * bw);

  /*
    generate sin values (need if order of transform is odd)
  */
  ArcCosEvalPts(n,eval_args);
  for (i=0; i<n; i++)
    sin_values[i] = sin(eval_args[i]);

  /* make DCT plans -> note the arrays!!! since I'll be
     using the guru fftw interface, I don't care what the
     arrays are - I just want to make sure they're big enough */

  /* forward DCT */
  forwardDCT = fftw_plan_r2r_1d( 2*bw, workspace, weights,
				 FFTW_REDFT10, FFTW_ESTIMATE ) ;
      
  /* inverse DCT */
  inverseDCT = fftw_plan_r2r_1d( 2*bw, workspace, weights,
				 FFTW_REDFT01, FFTW_ESTIMATE );

  /* make the weights */
  makeweights( bw, weights );

  /*
    precompute cosine series for Pml (Gml) functions
    necessary for the forward transform
  */
  CosPmlTableGen(bw,
		 m,
		 cos_pml_table,
		 workspace);
  /*
    take transpose of the precomputed cosine series,
    necessary for the inverse transform
  */
  Transpose_CosPmlTableGen(bw,
			   m,
			   cos_pml_table,
			   transpose_cos_pml_table );

  sum_error = 0.0 ;
  sum_relerror = 0.0 ;

  totaltimef = 0.0 ;
  totaltimei = 0.0 ;

  /* generate seed for random number generator */
  time ( &seed ) ;
  srand48( seed ) ;


  /* now the for-loop */
  for ( k = 0 ; k < loops ; k ++ )
    {
      /* generate random coefficients */
      for( i = 0 ; i < (bw - m) ; i ++ )
	coeffs[ i ] = 2.0 * ( drand48() - 0.5 ) ;

      tstarti = csecond() ;
      /* do inverse semi-naive transform */
      InvSemiNaiveReduced(coeffs,
			  bw,
			  m,
			  samples,
			  transpose_cos_pml_table,
			  sin_values,
			  workspace,
			  &inverseDCT );
      tstopi = csecond();
      totaltimei += (tstopi - tstarti) ;

      tstartf = csecond() ;
      /* now do forward semi-naive transform */
      SemiNaiveReduced(samples,
		       bw,
		       m,
		       newcoeffs,
		       workspace,
		       cos_pml_table,
		       weights,
		       &forwardDCT );
      tstopf = csecond();
      totaltimef += (tstopf - tstartf) ;

      /*
	now tally up the error between the original
	coefficients and the new ones
      */
	 
      relerror[ k ] = 0.0 ;
      curmax[ k ] = 0.0 ;
      /* now figure out errors */
      for( j = 0 ; j < bw - m ; j ++ )
	{
	  tmp_error = fabs( coeffs[j] - newcoeffs[j] );
	  tmp_relerror = tmp_error / ( fabs( coeffs[j] ) +
				       pow( 10.0, -50.0 ) );
	  curmax[ k ] = mymax( curmax[ k ], tmp_error );
	  relerror[ k ] = mymax( relerror[ k ], tmp_relerror );
	}
      sum_error += curmax[ k ] ;
      sum_relerror += relerror[ k ] ;

    }

  ave_error = sum_error / ( (double) loops );
  ave_relerror = sum_relerror / ( (double) loops );
  stddev_error = 0.0 ; stddev_relerror = 0.0;
  for( i = 0 ; i < loops ; i ++ )
    {
      stddev_error += pow( ave_error - curmax[ i ] , 2.0 );
      stddev_relerror += pow( ave_relerror - relerror[ i ] , 2.0 );
    }
  /*** this won't work if loops == 1 ***/
  if( loops != 1 )
    {
      stddev_error = sqrt(stddev_error / ( (double) (loops - 1) ) );
      stddev_relerror = sqrt(stddev_relerror / ( (double) (loops - 1) ) );
    }

  fprintf(stdout,"bw = %d\tm = %d\n",bw, m);
  fprintf(stdout,"loops = %d\n", loops );

  fprintf(stdout,"Average r-o error:\t\t %.4e\t",
	  sum_error/((double) loops));
  fprintf(stdout,"std dev: %.4e\n",stddev_error);
  fprintf(stdout,"Average (r-o)/o error:\t\t %.4e\t",
	  sum_relerror/((double) loops));
  fprintf(stdout,"std dev: %.4e\n\n",stddev_relerror);

  fprintf(stdout,"average forward time = %.4e\n",
	  totaltimef/((double) loops) );
  fprintf(stdout,"average inverse time = %.4e\n",
	  totaltimei/((double) loops) );


  fftw_destroy_plan( inverseDCT );
  fftw_destroy_plan( forwardDCT );


  free(weights);
  free(workspace) ;
  free(curmax) ;
  free(relerror) ;
  free(error) ;
  free(sin_values );
  free(eval_args) ;
  free(transpose_cos_pml_table) ;
  free(cos_pml_table) ;
  free(newcoeffs);
  free(coeffs);
  free(samples);

  return 0 ;

}

