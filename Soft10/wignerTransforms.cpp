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
 
  functions involved in doing Wigner (little d) transforms:

  wigNaiveAnalysis(): forward wigner transform (spatial -> spectral)
  wigNaiveSynthesis(): inverse wigner transform (spectral -> spatial)

*/

#include "utils_so3.h"
#include "weights.h"


/****************************************************************/
/**
   wigNaiveAnalysis

   Given a bandwidth bw, and orders m1, m2, this function will analyse
   the signal of length n = 2*bw by projecting it onto the L2-normed
   little d wigners of degrees l = Max(|m1|, |m2|) ... bw - 1, i.e. I
   am doing the integrals

   <signal, d_{m1,m2}^l> for l = Max(|m1|, |m2|) ... bw - 1

   where the little d is normalized such that

   \int_0^pi (d_{m1,m2}^l) (d_{m1,m2}^{lx}) \sin\theta d\theta = \delta_{l,lx}

   NOTE: to use this routine as part of the full transform on SO(3), I still
   need to multiply these coefficients by (PI/bw) ... the PI is because
   my basis functions are the big D's now, the bw because of the adjustment
   in scaling doing the Fourier transform introduces. This multiplying HAS
   TO HAPPEN outside of this routine.


   NOTE that this function was written to be part of a full
   SO(3) harmonic transform, so a certain amount of precomputation
   has been assumed.

   arguments: m1, m2 = orders of the transform
              bw = bandwidth
	      signal = array of sample values, length n = 2*bw
	      wigners = array of length (bw - max(|m1|,|m2|))*( 2 * bw );
                        this array holds the wigner little d's that I
			need - PRECOMPUTED by genWig_L2() 
	      coeffs = array of length (bw - max(|m1|,|m2|)), to store the
	               final results
	      workspace = scratch area, of length n

**/

void wigNaiveAnalysis( int m1,
		       int m2,
		       int bw,
		       double *signal,
		       double *wigners,
		       double *coeffs,
		       double *workspace )
{
  int i, j, l, m, n ;
  double tmp, *wignersPtr ;
  double *weightedSignal;
  const double *weights ;

  /* l is the degree of the "first" wigner function at
     this order */
  l = MAX( ABS( m1 ) , ABS( m2 ) ) ;
  m = l ;
  n = 2 * bw ;


  weightedSignal = workspace ;
  wignersPtr = wigners ;

  /* weight the signal with the appropriate quadrature weights */
  weights = get_weights( bw ) ;
  vec_pt_mul_so3( (double *) weights, signal, weightedSignal, n ) ;
  /****
  for( i = 0 ; i < n ; i ++ )
    weightedSignal[ i ] = weights[ i ] * signal[ i ] ;
  ****/


  /* and now I start analysing, i.e. multiply the matrix "wigners" by
     the vector "weightedSignal" */

  for( i = 0 ; i < bw - m ; i ++ )
    {
      tmp = 0.0 ;
      for( j = 0 ; j < n ; j ++ )
	tmp += *wignersPtr++ * weightedSignal[ j ] ;
      coeffs[ i ] = tmp ;
    }

  /* and that should be all */

}

/****************************************************************/
/**
   wigNaiveSynthesis

   Given a bandwidth bw, and orders m1, m2, this function will synthesize
   the signal of length n = 2*bw by "summing up the coefficients". More
   plainly, this is the inverse transform of wigNaiveAnalysis.

   Let l = Max(|m1|, |m2|). In matrix-lingo, wigNaiveAnalysis may be
   written as:

   c = P W f

   where f is the data vector, W is the quadrature matrix (i.e. weights),
   P is the (bw-l) x n matrix of sample values of the L2-normed wigners
   d_{m1,m2}^l d_{m1,m2}^{l+1} ... d_{m1,m2}^{bw-1}, and c is the
   wigner series representation (i.e. coefficients) of f (c is
   a vector of length bw-l).

   So wigNaiveSynthesis can be written as

   f = Transpose(P) c

   No quadrature matrix is necessary.

   NOTE that this function was written to be part of a full
   SO(3) harmonic transform, so a certain amount of precomputation
   has been assumed.

   NOTE: to use this routine as part of the full transform on SO(3), I still
   need to multiply these coefficients by (PI/bw) ... the PI is because
   my basis functions are the big D's now, the bw because of the adjustment
   in scaling doing the Fourier transform introduces. This multiplying HAS
   TO HAPPEN outside of this routine.

   arguments: m1, m2 = orders of the transform
              bw = bandwidth
	      coeffs = array of coefficients, length (bw - max(|m1|,|m2|))
	      wignersTrans = array of length (bw - max(|m1|,|m2|))*( 2 * bw );
                             this array holds the wigner little d's that I
			     need at this order - PRECOMPUTED genWigTrans_L2()
	      signal = array of length n = 2*bw, to store the final results,
                       the reconstructed sample values
	      workspace = scratch area, of length 0 * n
                          (that's right, 0 * n ... I'm keeping this
			   argument here just so that this function
			   call looks identical to wigNaiveAnalysis)

**/

void wigNaiveSynthesis( int m1,
			int m2,
			int bw,
			double *coeffs,
			double *wignersTrans,
			double *signal,
			double *workspace )
{
  int i, j, m, n ;
  double tmp, *wignersTransPtr ;

  /* l is the degree of the "first" wigner function at
     this order */
  m = MAX( ABS( m1 ) , ABS( m2 ) ) ;
  n = 2 * bw ;

  wignersTransPtr = wignersTrans ;

  /* just do the sums */
  for ( i = 0 ; i < n ; i ++ )
    {
      tmp = 0.0 ;
      for ( j = 0 ; j < (bw - m) ; j ++ )
	tmp += *wignersTransPtr++ * coeffs[j] ;
      signal[ i ] = tmp ;
    }

  /* that's it */
}
