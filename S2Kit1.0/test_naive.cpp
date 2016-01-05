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


/*************************************************************************/

/* test_naive.c - top-level code for computing forward and inversee
   Legendre transforms using the naive algorithm; will return timing
   and error information

   m    - order of the problem
   bw   - bandwidth
   loops - number of loops thru timed portion of code.  Intended
           to reduce noise due to multiprocessing and 
	   discretization errors
   
   Sample calls:

   test_naive m bw loops

   test_naive 0 32 10

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "pmls.h"
#include "makeweights.h"
#include "naive_synthesis.h"
#include "csecond.h"

#define mymax(a, b) ((a) > (b) ? (a) : (b))


int main(argc, argv)
     int argc;
     char **argv;
{
  int i, j, k;
  int bw, m, loops;
  double *samples, *coeffs, *newcoeffs ;
  double *plm, *weights, *workspace ;
  double sum_error, sum_relerror ;
  double tmp_error, ave_error ;
  double tmp_relerror, ave_relerror ;
  double stddev_error, stddev_relerror;
  double *error, *relerror, *curmax ;
  double tstartf, tstopf, totaltimef ;
  double tstarti, tstopi, totaltimei ;
  long int seed ;

  if (argc < 4)
    {
      fprintf(stdout,"Usage: test_naive m bw loops\n");
      return(0);
    }

  m = atoi(argv[1]);
  bw = atoi(argv[2]);
  loops = atoi(argv[3]);

  /* space for samples and coefficients */
  samples = (double *) malloc(sizeof(double) * 2 * bw);
  coeffs = (double *) malloc(sizeof(double) * (bw - m) );
  newcoeffs = (double *) malloc(sizeof(double) * (bw - m) );

  /* space for precomputed Plms */
  plm = (double *) malloc(sizeof(double) * 2 * bw * (bw - m) );

  /* for weights */
  weights = (double *) malloc(sizeof(double) * 4 * bw);

  /* workspace space */
  workspace = (double *) malloc(sizeof(double) * 18 * bw);

  /* space for error stuff */
  error = (double *) malloc(sizeof(double) * loops);
  relerror = (double *) malloc(sizeof(double) * loops);
  curmax = (double *) malloc(sizeof(double) * loops);
  
  sum_error = 0.0 ;
  sum_relerror = 0.0 ;

  /* precompute the Plms */
  PmlTableGen( bw, m, plm, workspace ) ;

  totaltimef = 0.0 ;
  totaltimei = 0.0 ;

  /* generate seed for random number generator */
  time ( &seed ) ;
  srand48( seed ) ;

  /* make the weights */
  makeweights( bw, weights );

  /* now the for-loop */
  for ( k = 0 ; k < loops ; k ++ )
    {
      /* generate random coefficients */
      for( i = 0 ; i < (bw - m) ; i ++ )
	coeffs[ i ] = 2.0 * ( drand48() - 0.5 ) ;
     
      tstarti = csecond() ;
      /* do inverse naive transform */
      Naive_SynthesizeX(coeffs,
			bw,
			m,
			samples,
			plm ) ;
      tstopi = csecond();
      totaltimei += (tstopi - tstarti) ;

      tstartf = csecond() ;
      /* now do forward naive transform */
      Naive_AnalysisX( samples,
		       bw,
		       m,
		       weights,
		       newcoeffs,
		       plm,
		       workspace );
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


  /* now free up memory */
  free( curmax );
  free( relerror );
  free( error );
  free( workspace );
  free( weights );
  free( newcoeffs );
  free( coeffs );
  free( samples );

  return 0 ;
}










