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

  functions that do the FORWARD and INVERSE transforms on the
  full group SO(3)

  - uses Wigner-d symmetries

  sample size = (2*bw)^3
  coefficient size = on the order of bw^3 (a little more, actually)

  functions in here:

  Forward_SO3_Naive_sym();
  Inverse_SO3_Naive_sym();

*/

/* the following is for memset */
#include <string.h>

#include <math.h>
#include "utils_so3.h"
#include "fft_grids_so3.h"
#include "makeWigner.h"
#include "wignerTransforms.h"
#include "wignerTransforms_sym.h"

/*
  ok, the forward transform

  Function arguments are as follows:

  bw = bandwidth of transform
  rdata, idata: real and imaginary parts of the input signal,
                EACH array of size (2*bw)^3
  rcoeffs, icoeffs: real and imaginary parts of the coefficients,
                    EACH array of size bw*(4*bw^2-1)/3
  workspace1: array for tmp storage, of size (gulp) 4 * (2*bw)^3
  workspace2: another array for tmp storage, of size 24*bw + 2*bw^2 ;

  flag: = 0 : data is COMPLEX
        = 1 : data is REAL

*/

void Forward_SO3_Naive_sym( int bw,
			    double *rdata, double *idata,
			    double *rcoeffs, double *icoeffs,
			    double *workspace1, double *workspace2,
			    int flag )
{
  int j, n ;
  int m1, m2;
  double fudge ;
  int sampHere, coefHere, coefHere2 ;
  int tCoeffs ;
  double *t1r, *t1i, *t2r, *t2i ;
  double *sinPts, *cosPts, *sinPts2, *cosPts2 ;
  double *wigners, *scratch ;
  double *rcoeffsPtr, *icoeffsPtr ;
  double *rdataPtr, *idataPtr ;
  double dn ;

  n = 2 * bw ;
  t1r = workspace1 ;
  t1i = workspace1 + (n * n * n) ;
  t2r = t1i + (n * n * n) ;
  t2i = t2r + (n * n * n) ;

  /* I'll need these for later */
  rcoeffsPtr = rcoeffs ;
  icoeffsPtr = icoeffs ;
  rdataPtr = rdata ;
  idataPtr = idata ;

  sinPts = workspace2 ;
  cosPts = sinPts + n ;
  sinPts2 = cosPts + n ;
  cosPts2 = sinPts2 + n ;
  wigners = cosPts2 + n ;
  scratch = wigners + ( bw * n ) ; /* wigners need at most bw*n space AT
				      ANY given orders m1, m2 */

  tCoeffs = totalCoeffs_so3( bw ) ;

  /*
    before going further, let's precompute all the sines
    and cosines I'll need. No matter what order transform
    I'm doing, these'll stay the same.
  */
  SinEvalPts( n, sinPts );
  CosEvalPts( n, cosPts );
  SinEvalPts2( n, sinPts2 );
  CosEvalPts2( n, cosPts2 );

    
  /*
    Stage 1: FFT the "rows". Instead of treating the signal as
    3-D object, I can also think of it as an array of size
    (n^2) x n. This means all I'm doing in the first stage
    is taking n^2-many FFTs, each of length n.

    NOTE: Since I'm reusing the FFT code from SpharmonicKit,
    even though I'm doing the FORWARD SO(3) transform
    here, I need to call grid_invfourier_so3  -> the signs
    on the complex exponentials are switched (detailed
    explanation to be put here eventually, but trust
    me)
  */
  
  grid_invfourier_so3( rdata, idata,
		       t1r, t1i,
		       n*n, n,
		       scratch ) ;

  /* normalize the Fourier coefficients (sorry, have to do it) */
  /*
    But wait! Instead of normalizing here, and right after the second
    fft, I'm going to wait till the end, and combine the 3 normalizations
    I have to do. Cut down on the for-loops ...
  */
  /*
    dn = 1. / sqrt( (double) n ) ;
    for ( j = 0 ; j < n*n*n; j++ )
    {
    t1r[ j ] *= dn ;
    t1i[ j ] *= dn ;
    }
  */

  /*
    Stage 2: transpose!
  */
  
  transpose_so3( t1r, t2r, n*n, n ) ;
  transpose_so3( t1i, t2i, n*n, n ) ;

  /*
    Stage 3: FFT again. Note that I'm using the tmp space
    of t1r, t1i again
  */

  grid_invfourier_so3( t2r, t2i,
		       t1r, t1i,
		       n*n, n,
		       scratch ) ;

  /* normalize the Fourier coefficients (sorry, have to do it) */
  /*
    for ( j = 0 ; j < n*n*n; j++ )
    {
    t1r[ j ] *= dn ;
    t1i[ j ] *= dn ;
    }
  */

  /*
    Stage 4: transpose again! And note I'm using the tmp space
    of t2r, t2i again.
  */
  
  transpose_so3( t1r, t2r, n*n, n ) ;
  transpose_so3( t1i, t2i, n*n, n ) ;
 
  /*
    Stage 5: Do the Wigner transforms. This is the tricky bit.

    Since I'm working with two order indeces, m1 and m2, the
    for-loops will be more intricate than in the case of the
    "ordinary" spherical transform on S^2.

    Also, I will be taking advantage of the symmetries of the
    Wigner-d functions. As long as I keep my signs and flips
    right, the Wigner-d's I precompute for an order (m1, m2)
    transform can generally  be used in seven more transforms:
    (m1,-m2), (m2,m1), (m2,-m1), (-m2,m1), (-m2,-m1), (-m1,m2)
    and (-m1,-m2).

    I say "general" because, of course, I'll be transforming
    at orders (m1,m1), (m1,0) and (0,m1), so I won't get such
    a huge reduction. Still, this should save time.

    If assumptions are made regarding the original input signal,
    e.g. it's strictly real, then one may take advantage of
    symmetries of the big D wigners (i.e. function of all 3
    parameters alpha, beta, gamma) and so simplify the for-loops
    some and hence increase the speed of the program. However,
    the for-loops to follow will make no such assumptions.
    Whether the signal is real or complex, these for-loops will
    handle it.

    The for-loops will be "designed" as follows. They will be
    divided into cases according to the orders:


    0) {f_{0,0}}

    1) for 0 <= m1 <= bw-1
    compute the coefficients
    i)   {f_{ m1, m1}}
    ii)  {f_{-m1,-m1}}
    iii) {f_{-m1, m1}}
    iv)  {f_{ m1,-m1}}

    2) for 1 <= m1 <= bw-1
    compute the coefficients
    i)   {f_{ m1,  0}}
    ii)  {f_{-m1,  0}}
    iii) {f_{  0, m1}}
    iv)  {f_{  0,-m1}}

    3) for 1 <= m1 <= bw-1
    for m1+1 <= m2 <= bw-1
    compute the coefficients
    i)    {f_{ m1, m2}}
    ii)   {f_{-m1,-m2}}
    iii)  {f_{ m1,-m2}}
    iv)   {f_{-m1, m2}}
    v)    {f_{ m2, m1}}
    vi)   {f_{-m2,-m1}}
    vii)  {f_{ m2,-m1}}
    viii) {f_{-m2, m1}}

    Fasten your seatbelt, folks. It's going to be a bumpy ride.

  */


  /***************************/
  /*                         */
  /* {f_{0,0}} coefficient   */
  /*                         */
  /***************************/


  /* compute the wigners I'll need */
  genWig_L2( 0, 0, bw,
	     sinPts, cosPts,
	     sinPts2, cosPts2,
	     wigners, scratch ) ;
  
  /* now, get the locations of where the
     samples I have to transform are, and
     where the coefficients have to go */
  
  sampHere = sampLoc_so3( 0, 0, bw ) ;
  coefHere = coefLoc_so3( 0, 0, bw ) ;

  /* ok, reset sample, coef ptrs */
  rcoeffsPtr = rcoeffs ;
  icoeffsPtr = icoeffs ;
  rdataPtr = t2r ;
  idataPtr = t2i ;
  
  /* now advance by the computed amounts */
  rdataPtr += sampHere ;
  idataPtr += sampHere ;
  rcoeffsPtr += coefHere ;
  icoeffsPtr += coefHere ;
  
  /* now transform the real and imaginary parts
     of the data */
  
  wigNaiveAnalysis( 0, 0, bw, rdataPtr,
		    wigners, rcoeffsPtr,
		    scratch ) ;
  
  wigNaiveAnalysis( 0, 0, bw, idataPtr,
		    wigners, icoeffsPtr,
		    scratch ) ;



  /*** 0 <= m1 <= bw-1 ***/
  for ( m1 = 1 ; m1 < bw ; m1 ++ )
    {

      /* compute the wigners I'll need */
      genWig_L2( m1, m1, bw,
		 sinPts, cosPts,
		 sinPts2, cosPts2,
		 wigners, scratch ) ;

      /***************************/
      /*                         */
      /* {f_{m1,m1}} coefficient */
      /*                         */
      /***************************/


      /* now, get the locations of where the
	 samples I have to transform are, and
	 where the coefficients have to go */

      sampHere = sampLoc_so3( m1, m1, bw ) ;
      coefHere = coefLoc_so3( m1, m1, bw ) ;

      /* ok, reset sample, coef ptrs */
      rcoeffsPtr = rcoeffs ;
      icoeffsPtr = icoeffs ;
      rdataPtr = t2r ;
      idataPtr = t2i ;
      
      /* now advance by the computed amounts */
      rdataPtr += sampHere ;
      idataPtr += sampHere ;
      rcoeffsPtr += coefHere ;
      icoeffsPtr += coefHere ;

      /* now transform the real and imaginary parts
	 of the data */

      wigNaiveAnalysis( m1, m1, bw, rdataPtr,
			wigners, rcoeffsPtr,
			scratch ) ;
      
      wigNaiveAnalysis( m1, m1, bw, idataPtr,
			wigners, icoeffsPtr,
			scratch ) ;

      /*****************************/
      /*                           */
      /* {f_{-m1,-m1}} coefficient */
      /*                           */
      /*****************************/

      if ( flag == 0 ) /* if data is complex */
	{

	  /* now, get the locations of where the
	     samples I have to transform are, and
	     where the coefficients have to go */

	  sampHere = sampLoc_so3( -m1, -m1, bw ) ;
	  coefHere = coefLoc_so3( -m1, -m1, bw ) ;

	  /* ok, reset sample, coef ptrs */
	  rcoeffsPtr = rcoeffs ;
	  icoeffsPtr = icoeffs ;
	  rdataPtr = t2r ;
	  idataPtr = t2i ;
      
	  /* now advance by the computed amounts */
	  rdataPtr += sampHere ;
	  idataPtr += sampHere ;
	  rcoeffsPtr += coefHere ;
	  icoeffsPtr += coefHere ;

	  /* now transform the real and imaginary parts
	     of the data */

	  wigNaiveAnalysis( -m1, -m1, bw, rdataPtr,
			    wigners, rcoeffsPtr,
			    scratch ) ;
      
	  wigNaiveAnalysis( -m1, -m1, bw, idataPtr,
			    wigners, icoeffsPtr,
			    scratch ) ;

	}
      else  /* data is real, so use symmetry */
	{
	  coefHere = coefLoc_so3( m1, m1, bw ) ;
	  coefHere2 = coefLoc_so3( -m1, -m1, bw ) ;

	  for ( j = 0 ; j < bw - m1 ; j ++ )
	    {
	      rcoeffs[coefHere2+j] = rcoeffs[coefHere+j];
	      icoeffs[coefHere2+j] = -icoeffs[coefHere+j];
	    }

	}


      /*****************************/
      /*                           */
      /* {f_{-m1,m1}} coefficient  */
      /*                           */
      /*****************************/


      /* now, get the locations of where the
	 samples I have to transform are, and
	 where the coefficients have to go */

      sampHere = sampLoc_so3( -m1, m1, bw ) ;
      coefHere = coefLoc_so3( -m1, m1, bw ) ;

      /* ok, reset sample, coef ptrs */
      rcoeffsPtr = rcoeffs ;
      icoeffsPtr = icoeffs ;
      rdataPtr = t2r ;
      idataPtr = t2i ;
      
      /* now advance by the computed amounts */
      rdataPtr += sampHere ;
      idataPtr += sampHere ;
      rcoeffsPtr += coefHere ;
      icoeffsPtr += coefHere ;

      /* now transform the real and imaginary parts
	 of the data */

      wigNaiveAnalysis_symY( -m1, m1, bw, rdataPtr,
			     wigners, rcoeffsPtr,
			     scratch ) ;
      
      wigNaiveAnalysis_symY( -m1, m1, bw, idataPtr,
			     wigners, icoeffsPtr,
			     scratch ) ;


      /*****************************/
      /*                           */
      /* {f_{m1,-m1}} coefficient  */
      /*                           */
      /*****************************/

      if ( flag == 0 ) /* data is complex */
	{

	  /* now, get the locations of where the
	     samples I have to transform are, and
	     where the coefficients have to go */

	  sampHere = sampLoc_so3( m1, -m1, bw ) ;
	  coefHere = coefLoc_so3( m1, -m1, bw ) ;

	  /* ok, reset sample, coef ptrs */
	  rcoeffsPtr = rcoeffs ;
	  icoeffsPtr = icoeffs ;
	  rdataPtr = t2r ;
	  idataPtr = t2i ;
      
	  /* now advance by the computed amounts */
	  rdataPtr += sampHere ;
	  idataPtr += sampHere ;
	  rcoeffsPtr += coefHere ;
	  icoeffsPtr += coefHere ;

	  /* now transform the real and imaginary parts
	     of the data */

	  wigNaiveAnalysis_symY( m1, -m1, bw, rdataPtr,
				 wigners, rcoeffsPtr,
				 scratch ) ;
      
	  wigNaiveAnalysis_symY( m1, -m1, bw, idataPtr,
				 wigners, icoeffsPtr,
				 scratch ) ;

	}
      else /* data is real, so use symmetry */
	{
	  coefHere = coefLoc_so3( -m1, m1, bw );
	  coefHere2 = coefLoc_so3( m1, -m1, bw );

	  for ( j = 0 ; j < bw - m1 ; j ++ )
	    {
	      rcoeffs[coefHere2+j] = rcoeffs[coefHere+j];
	      icoeffs[coefHere2+j] = -icoeffs[coefHere+j];
	    }
	}
    }

  /*** for 1 <= m1 <= bw-1 ***/
  for ( m1 = 1 ; m1 < bw ; m1 ++ )
    {
      
      /* compute the wigners I'll need */
      genWig_L2( m1, 0, bw,
		 sinPts, cosPts,
		 sinPts2, cosPts2,
		 wigners, scratch ) ;


      /***************************/
      /*                         */
      /* {f_{m1,0}} coefficient */
      /*                         */
      /***************************/


      /* get the locations of where the
	 samples I have to transform are, and
	 where the coefficients have to go */

      sampHere = sampLoc_so3( m1, 0, bw ) ;
      coefHere = coefLoc_so3( m1, 0, bw ) ;

      /* ok, reset sample, coef ptrs */
      rcoeffsPtr = rcoeffs ;
      icoeffsPtr = icoeffs ;
      rdataPtr = t2r ;
      idataPtr = t2i ;
      
      /* now advance by the computed amounts */
      rdataPtr += sampHere ;
      idataPtr += sampHere ;
      rcoeffsPtr += coefHere ;
      icoeffsPtr += coefHere ;

      /* now transform the real and imaginary parts
	 of the data */

      wigNaiveAnalysis( m1, 0, bw, rdataPtr,
			wigners, rcoeffsPtr,
			scratch ) ;
      
      wigNaiveAnalysis( m1, 0, bw, idataPtr,
			wigners, icoeffsPtr,
			scratch ) ;

      /***************************/
      /*                         */
      /* {f_{-m1,0}} coefficient */
      /*                         */
      /***************************/

      if ( flag == 0 ) /* data is complex */
	{
      
	  /* get the locations of where the
	     samples I have to transform are, and
	     where the coefficients have to go */

	  sampHere = sampLoc_so3( -m1, 0, bw ) ;
	  coefHere = coefLoc_so3( -m1, 0, bw ) ;

	  /* ok, reset sample, coef ptrs */
	  rcoeffsPtr = rcoeffs ;
	  icoeffsPtr = icoeffs ;
	  rdataPtr = t2r ;
	  idataPtr = t2i ;
      
	  /* now advance by the computed amounts */
	  rdataPtr += sampHere ;
	  idataPtr += sampHere ;
	  rcoeffsPtr += coefHere ;
	  icoeffsPtr += coefHere ;

	  /* now transform the real and imaginary parts
	     of the data */

	  wigNaiveAnalysis_symX( -m1, 0, bw, rdataPtr,
				 wigners, rcoeffsPtr,
				 scratch ) ;
      
	  wigNaiveAnalysis_symX( -m1, 0, bw, idataPtr,
				 wigners, icoeffsPtr,
				 scratch ) ;

	}
      else  /* data is real, so use symmetry */
	{
	  coefHere = coefLoc_so3( m1, 0, bw );
	  coefHere2 = coefLoc_so3( -m1, 0, bw );

	  if ( (m1 % 2) == 0 )
	    fudge = 1.0 ;
	  else
	    fudge = -1.0 ;

	  for ( j = 0 ; j < bw - m1 ; j ++ )
	    {
	      rcoeffs[coefHere2+j] = fudge * rcoeffs[coefHere+j];
	      icoeffs[coefHere2+j] = -fudge * icoeffs[coefHere+j];
	    }
	  
	}

      /***************************/
      /*                         */
      /* {f_{0,m1}} coefficient */
      /*                         */
      /***************************/

      /* get the locations of where the
	 samples I have to transform are, and
	 where the coefficients have to go */

      sampHere = sampLoc_so3( 0, m1, bw ) ;
      coefHere = coefLoc_so3( 0, m1, bw ) ;

      /* ok, reset sample, coef ptrs */
      rcoeffsPtr = rcoeffs ;
      icoeffsPtr = icoeffs ;
      rdataPtr = t2r ;
      idataPtr = t2i ;
      
      /* now advance by the computed amounts */
      rdataPtr += sampHere ;
      idataPtr += sampHere ;
      rcoeffsPtr += coefHere ;
      icoeffsPtr += coefHere ;

      /* now transform the real and imaginary parts
	 of the data */

      wigNaiveAnalysis_symX( 0, m1, bw, rdataPtr,
			     wigners, rcoeffsPtr,
			     scratch ) ;
      
      wigNaiveAnalysis_symX( 0, m1, bw, idataPtr,
			     wigners, icoeffsPtr,
			     scratch ) ;


      /***************************/
      /*                         */
      /* {f_{0,-m1}} coefficient */
      /*                         */
      /***************************/

      if ( flag == 0 ) /* data is complex */
	{
      
	  /* get the locations of where the
	     samples I have to transform are, and
	     where the coefficients have to go */

	  sampHere = sampLoc_so3( 0, -m1, bw ) ;
	  coefHere = coefLoc_so3( 0, -m1, bw ) ;

	  /* ok, reset sample, coef ptrs */
	  rcoeffsPtr = rcoeffs ;
	  icoeffsPtr = icoeffs ;
	  rdataPtr = t2r ;
	  idataPtr = t2i ;
      
	  /* now advance by the computed amounts */
	  rdataPtr += sampHere ;
	  idataPtr += sampHere ;
	  rcoeffsPtr += coefHere ;
	  icoeffsPtr += coefHere ;

	  /* now transform the real and imaginary parts
	     of the data */

	  wigNaiveAnalysis( 0, -m1, bw, rdataPtr,
			    wigners, rcoeffsPtr,
			    scratch ) ;
      
	  wigNaiveAnalysis( 0, -m1, bw, idataPtr,
			    wigners, icoeffsPtr,
			    scratch ) ;

	}
      else  /* data is real, so use symmetry */
	{
   	  coefHere = coefLoc_so3( 0, m1, bw );
	  coefHere2 = coefLoc_so3( 0, -m1, bw );

	  if ( (m1 % 2) == 0 )
	    fudge = 1.0 ;
	  else
	    fudge = -1.0 ;

	  for ( j = 0 ; j < bw - m1 ; j ++ )
	    {
	      rcoeffs[coefHere2+j] = fudge * rcoeffs[coefHere+j];
	      icoeffs[coefHere2+j] = -fudge * icoeffs[coefHere+j];
	    }
	  
	}

    }


  /***
      1 <= m1 <= bw-1
      m1+1 <= m2 <= bw-1
  ***/

  for ( m1 = 1 ; m1 < bw ; m1 ++ )
    {
      for ( m2 = m1 + 1 ; m2 < bw ; m2 ++ )
	{

	  
	  /* compute the wigners I'll need */
	  genWig_L2( m1, m2, bw,
		     sinPts, cosPts,
		     sinPts2, cosPts2,
		     wigners, scratch ) ;


	  /***************************/
	  /*                         */
	  /* {f_{m1,m2}} coefficient */
	  /*                         */
	  /***************************/

	  /* get the locations of where the
	     samples I have to transform are, and
	     where the coefficients have to go */
	  
	  sampHere = sampLoc_so3( m1, m2, bw ) ;
	  coefHere = coefLoc_so3( m1, m2, bw ) ;
	  
	  /* ok, reset sample, coef ptrs */
	  rcoeffsPtr = rcoeffs ;
	  icoeffsPtr = icoeffs ;
	  rdataPtr = t2r ;
	  idataPtr = t2i ;
	  
	  /* now advance by the computed amounts */
	  rdataPtr += sampHere ;
	  idataPtr += sampHere ;
	  rcoeffsPtr += coefHere ;
	  icoeffsPtr += coefHere ;
	  
	  /* now transform the real and imaginary parts
	     of the data */
	  
	  wigNaiveAnalysis( m1, m2, bw, rdataPtr,
			    wigners, rcoeffsPtr,
			    scratch ) ;
	  
	  wigNaiveAnalysis( m1, m2, bw, idataPtr,
			    wigners, icoeffsPtr,
			    scratch ) ;

	  /*****************************/
	  /*                           */
	  /* {f_{-m1,-m2}} coefficient */
	  /*                           */
	  /*****************************/

	  if ( flag == 0 ) /* data is complex */
	    {

	      /* get the locations of where the
		 samples I have to transform are, and
		 where the coefficients have to go */
	  
	      sampHere = sampLoc_so3( -m1, -m2, bw ) ;
	      coefHere = coefLoc_so3( -m1, -m2, bw ) ;
	  
	      /* ok, reset sample, coef ptrs */
	      rcoeffsPtr = rcoeffs ;
	      icoeffsPtr = icoeffs ;
	      rdataPtr = t2r ;
	      idataPtr = t2i ;


	      /* now advance by the computed amounts */
	      rdataPtr += sampHere ;
	      idataPtr += sampHere ;
	      rcoeffsPtr += coefHere ;
	      icoeffsPtr += coefHere ;
	  
	      /* now transform the real and imaginary parts
		 of the data */
	  
	      wigNaiveAnalysis_symX( -m1, -m2, bw, rdataPtr,
				     wigners, rcoeffsPtr,
				     scratch ) ;
	  
	      wigNaiveAnalysis_symX( -m1, -m2, bw, idataPtr,
				     wigners, icoeffsPtr,
				     scratch ) ;

	    }
	  else  /* data is real, so use symmetry */
	    {
	      coefHere = coefLoc_so3( m1, m2, bw );
	      coefHere2 = coefLoc_so3( -m1, -m2, bw );
	  
	      if ( ((m2-m1) % 2) == 0 )
		fudge = 1.0 ;
	      else
		fudge = -1.0 ;

	      for ( j = 0 ; j < bw - m2 ; j ++ )
		{
		  rcoeffs[coefHere2+j] = fudge * rcoeffs[coefHere+j];
		  icoeffs[coefHere2+j] = -fudge * icoeffs[coefHere+j];
		}
	    }


	  /****************************/
	  /*                          */
	  /* {f_{m1,-m2}} coefficient */
	  /*                          */
	  /****************************/
  

	  /* get the locations of where the
	     samples I have to transform are, and
	     where the coefficients have to go */
	  
	  sampHere = sampLoc_so3( m1, -m2, bw ) ;
	  coefHere = coefLoc_so3( m1, -m2, bw ) ;

	  /* ok, reset sample, coef ptrs */
	  rcoeffsPtr = rcoeffs ;
	  icoeffsPtr = icoeffs ;
	  rdataPtr = t2r ;
	  idataPtr = t2i ;
      
	  /* now advance by the computed amounts */
	  rdataPtr += sampHere ;
	  idataPtr += sampHere ;
	  rcoeffsPtr += coefHere ;
	  icoeffsPtr += coefHere ;
	  
	  /* now transform the real and imaginary parts
	     of the data */
	  
	  wigNaiveAnalysis_symY( m1, -m2, bw, rdataPtr,
				 wigners, rcoeffsPtr,
				 scratch ) ;
	  
	  wigNaiveAnalysis_symY( m1, -m2, bw, idataPtr,
				 wigners, icoeffsPtr,
				 scratch ) ;


	  /*****************************/
	  /*                           */
	  /* {f_{-m1,m2}} coefficient  */
	  /*                           */
	  /*****************************/


	  if ( flag == 0 ) /* data is complex */
	    {

	      /* get the locations of where the
		 samples I have to transform are, and
		 where the coefficients have to go */
	  
	      sampHere = sampLoc_so3( -m1, m2, bw ) ;
	      coefHere = coefLoc_so3( -m1, m2, bw ) ;


	      /* ok, reset sample, coef ptrs */
	      rcoeffsPtr = rcoeffs ;
	      icoeffsPtr = icoeffs ;
	      rdataPtr = t2r ;
	      idataPtr = t2i ;
      
	      /* now advance by the computed amounts */
	      rdataPtr += sampHere ;
	      idataPtr += sampHere ;
	      rcoeffsPtr += coefHere ;
	      icoeffsPtr += coefHere ;
	  
	      /* now transform the real and imaginary parts
		 of the data */
	  
	      wigNaiveAnalysis_symY( -m1, m2, bw, rdataPtr,
				     wigners, rcoeffsPtr,
				     scratch ) ;
	  
	      wigNaiveAnalysis_symY( -m1, m2, bw, idataPtr,
				     wigners, icoeffsPtr,
				     scratch ) ;
	  

	    }
	  else  /* data is real, so use symmetry */
	    {
	      coefHere = coefLoc_so3( m1, -m2, bw );
	      coefHere2 = coefLoc_so3( -m1, m2, bw );

	      if ( ((m2-m1) % 2) == 0 )
		fudge = 1.0 ;
	      else
		fudge = -1.0 ;
	      
	      for ( j = 0 ; j < bw - m2 ; j ++ )
		{
		  rcoeffs[coefHere2+j] = fudge * rcoeffs[coefHere+j];
		  icoeffs[coefHere2+j] = -fudge * icoeffs[coefHere+j];
		}
	    }


	  /***************************/
	  /*                         */
	  /* {f_{m2,m1}} coefficient */
	  /*                         */
	  /***************************/
	  
	  /* get the locations of where the
	     samples I have to transform are, and
	     where the coefficients have to go */
	  
	  sampHere = sampLoc_so3( m2, m1, bw ) ;
	  coefHere = coefLoc_so3( m2, m1, bw ) ;
	  
	  /* ok, reset sample, coef ptrs */
	  rcoeffsPtr = rcoeffs ;
	  icoeffsPtr = icoeffs ;
	  rdataPtr = t2r ;
	  idataPtr = t2i ;
	  
	  /* now advance by the computed amounts */
	  rdataPtr += sampHere ;
	  idataPtr += sampHere ;
	  rcoeffsPtr += coefHere ;
	  icoeffsPtr += coefHere ;
	  
	  /* now transform the real and imaginary parts
	     of the data */
	  
	  wigNaiveAnalysis_symX( m2, m1, bw, rdataPtr,
				 wigners, rcoeffsPtr,
				 scratch ) ;
	  
	  wigNaiveAnalysis_symX( m2, m1, bw, idataPtr,
				 wigners, icoeffsPtr,
				 scratch ) ;


	  /*****************************/
	  /*                           */
	  /* {f_{-m2,-m1}} coefficient */
	  /*                           */
	  /*****************************/

	  if ( flag == 0 ) /* data is complex */
	    {
	  
	      /* get the locations of where the
		 samples I have to transform are, and
		 where the coefficients have to go */
	  
	      sampHere = sampLoc_so3( -m2, -m1, bw ) ;
	      coefHere = coefLoc_so3( -m2, -m1, bw ) ;
	  
	      /* ok, reset sample, coef ptrs */
	      rcoeffsPtr = rcoeffs ;
	      icoeffsPtr = icoeffs ;
	      rdataPtr = t2r ;
	      idataPtr = t2i ;
	  
	      /* now advance by the computed amounts */
	      rdataPtr += sampHere ;
	      idataPtr += sampHere ;
	      rcoeffsPtr += coefHere ;
	      icoeffsPtr += coefHere ;
	  
	      /* now transform the real and imaginary parts
		 of the data */
	  
	      wigNaiveAnalysis( -m2, -m1, bw, rdataPtr,
				wigners, rcoeffsPtr,
				scratch ) ;
	  
	      wigNaiveAnalysis( -m2, -m1, bw, idataPtr,
				wigners, icoeffsPtr,
				scratch ) ;

	    }
	  else  /* data is real, so use symmetry */
	    {
	      coefHere = coefLoc_so3( m2, m1, bw );
	      coefHere2 = coefLoc_so3( -m2, -m1, bw );

	      if ( ((m2-m1) % 2) == 0 )
		fudge = 1.0 ;
	      else
		fudge = -1.0 ;

	      for ( j = 0 ; j < bw - m2 ; j ++ )
		{
		  rcoeffs[coefHere2+j] = fudge * rcoeffs[coefHere+j];
		  icoeffs[coefHere2+j] = -fudge * icoeffs[coefHere+j];
		}
	    }

	  /****************************/
	  /*                          */
	  /* {f_{m2,-m1}} coefficient */
	  /*                          */
	  /****************************/
  

	  /* get the locations of where the
	     samples I have to transform are, and
	     where the coefficients have to go */
	  
	  sampHere = sampLoc_so3( m2, -m1, bw ) ;
	  coefHere = coefLoc_so3( m2, -m1, bw ) ;
	  
	  /* ok, reset sample, coef ptrs */
	  rcoeffsPtr = rcoeffs ;
	  icoeffsPtr = icoeffs ;
	  rdataPtr = t2r ;
	  idataPtr = t2i ;
      
	  /* now advance by the computed amounts */
	  rdataPtr += sampHere ;
	  idataPtr += sampHere ;
	  rcoeffsPtr += coefHere ;
	  icoeffsPtr += coefHere ;
	  
	  /* now transform the real and imaginary parts
	     of the data */
	  
	  wigNaiveAnalysis_symY( m1, -m2, bw, rdataPtr,
				 wigners, rcoeffsPtr,
				 scratch ) ;
	  
	  wigNaiveAnalysis_symY( m1, -m2, bw, idataPtr,
				 wigners, icoeffsPtr,
				 scratch ) ;


	  /****************************/
	  /*                          */
	  /* {f_{-m2,m1}} coefficient */
	  /*                          */
	  /****************************/

	  if ( flag == 0 ) /* data is complex */
	    {  

	      /* get the locations of where the
		 samples I have to transform are, and
		 where the coefficients have to go */
	  
	      sampHere = sampLoc_so3( -m2, m1, bw ) ;
	      coefHere = coefLoc_so3( -m2, m1, bw ) ;
	  
	      /* ok, reset sample, coef ptrs */
	      rcoeffsPtr = rcoeffs ;
	      icoeffsPtr = icoeffs ;
	      rdataPtr = t2r ;
	      idataPtr = t2i ;
      
	      /* now advance by the computed amounts */
	      rdataPtr += sampHere ;
	      idataPtr += sampHere ;
	      rcoeffsPtr += coefHere ;
	      icoeffsPtr += coefHere ;
	  
	      /* now transform the real and imaginary parts
		 of the data */
	  
	      wigNaiveAnalysis_symY( -m1, m2, bw, rdataPtr,
				     wigners, rcoeffsPtr,
				     scratch ) ;
	  
	      wigNaiveAnalysis_symY( -m1, m2, bw, idataPtr,
				     wigners, icoeffsPtr,
				     scratch ) ;


	    }
	  else  /* data is real, so use symmetry */
	    {
	      coefHere = coefLoc_so3( m2, -m1, bw );
	      coefHere2 = coefLoc_so3( -m2, m1, bw );
	      
	      if ( ((m2-m1) % 2) == 0 )
		fudge = 1.0 ;
	      else
		fudge = -1.0 ;
	      
	      for ( j = 0 ; j < bw - m2 ; j ++ )
		{
		  rcoeffs[coefHere2+j] = fudge * rcoeffs[coefHere+j];
		  icoeffs[coefHere2+j] = -fudge * icoeffs[coefHere+j];
		}

	    }

	}
    }

  	  
  /* reset coef ptrs */
  rcoeffsPtr = rcoeffs ;
  icoeffsPtr = icoeffs ;

  /* need to normalize, one last time */
  /* This loop combines the two normalizations *not* done
     after the ffts near the start of this function */
  dn = (M_PI / ((double) (bw * n)));
  for ( j = 0 ; j < tCoeffs ; j ++ )
    {
      rcoeffsPtr[ j ] *= dn ;
      icoeffsPtr[ j ] *= dn ;
    }
  
  /*** and we're done ! ***/
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/

/*
  ok, the inverse transform

  Function arguments are as follows:

  bw = bandwidth of transform
  rcoeffs, icoeffs: real and imaginary parts of the coefficients,
                    EACH array of size (1/3 * bw * ( 4 bw^2 - 1 ))
  rdata, idata: real and imaginary parts of the synthesized (output)
                signal, EACH array of size (2*bw)^3
  workspace1: array for tmp storage, of size (gulp) 2 * (2*bw)^3

       (I don't need as much space as Forward_SO3_Naive because I'll
        be able to use the output rdata, idata arrays as tmp storage
	within this function. I couldn't use the rcoeffs, icoeffs arrays
	in Forward_SO3_Naive() because they weren't large enough ... I
	needed (2 bw)^3 space and they're only (1/3 * bw * ( 4 bw^2 - 1 ))
  workspace2: another array for tmp storage, of size 24*bw + 2*bw^2 ;

  flag: = 0 : data is COMPLEX
        = 1 : data is REAL

*/

void Inverse_SO3_Naive_sym( int bw,
			    double *rcoeffs, double *icoeffs,
			    double *rdata, double *idata,
			    double *workspace1, double *workspace2 ,
			    int flag )
{
  int j, n ;
  int m1, m2 ;
  int sampHere, sampHere2, coefHere ;
  double *t1r, *t1i ;
  double *sinPts, *cosPts, *sinPts2, *cosPts2 ;
  double *wignersTrans, *scratch ;
  double *rdataPtr, *idataPtr ;
  double *rcoeffsPtr, *icoeffsPtr ;
  double dn ;

  n = 2 * bw ;
  t1r = workspace1 ;
  t1i = workspace1 + (n * n * n) ;

  sinPts = workspace2 ;
  cosPts = sinPts + n ;
  sinPts2 = cosPts + n ;
  cosPts2 = sinPts2 + n ;
  wignersTrans = cosPts2 + n ;
  scratch = wignersTrans + ( bw * n ) ; /* wignersTrans need at most bw*n
					   space AT ANY given orders m1, m2 */

  /*
    before going further, let's precompute all the sines
    and cosines I'll need. No matter what order transform
    I'm doing, these'll stay the same.
  */
  SinEvalPts( n, sinPts );
  CosEvalPts( n, cosPts );
  SinEvalPts2( n, sinPts2 );
  CosEvalPts2( n, cosPts2 );


  /* Stage 0.5: Need to normalize the numbers before
     doing the IDWT
  */
  /* Nope! I'll do at the end of the function call,
     when I have to normalize for the two ffts anyway.
     This should cut down on the for-loops. */
  /*
    dn = ( ((double) bw) / M_PI ) ;
    dn = 1. ;
    for ( j = 0 ; j < totalCoeffs_so3( bw ) ; j++ )
    {
    t1r[ j ] = rcoeffs[ j ] * dn ;
    t1i[ j ] = icoeffs[ j ] * dn ;
    }
  */

  /*
    Stage 1: Do the Inverse Wigner transform. The rcoeffs, icoeffs
    arrays are assumed to be in the same "arrangement" as that produced
    by Forward_SO3_Naive().

    Since I'm working with two order indeces, m1 and m2, the
    for-loops will be more intricate than in the case of the
    "ordinary" spherical transform on S^2.

    Also, I will be taking advantage of the symmetries of the
    Wigner-d functions. As long as I keep my signs and flips
    right, the Wigner-d's I precompute for an order (m1, m2)
    transform can generally  be used in seven more transforms:
    (m1,-m2), (m2,m1), (m2,-m1), (-m2,m1), (-m2,-m1), (-m1,m2)
    and (-m1,-m2).


    The for-loops will be "designed" as follows. They will be
    divided into cases according to the orders:


    0) {f_{0,0}} inverse transform

    1) for 0 <= m1 <= bw-1
    compute inverse transform of
    i)   {f_{ m1, m1}}
    ii)  {f_{-m1,-m1}}
    iii) {f_{-m1, m1}}
    iv)  {f_{ m1,-m1}}

    2) for 1 <= m1 <= bw-1
    compute inverse transform of
    i)   {f_{ m1,  0}}
    ii)  {f_{-m1,  0}}
    iii) {f_{  0, m1}}
    iv)  {f_{  0,-m1}}

    3) for 1 <= m1 <= bw-1
    for m1+1 <= m2 <= bw-1
    compute inverse transform 
    i)    {f_{ m1, m2}}
    ii)   {f_{-m1,-m2}}
    iii)  {f_{ m1,-m2}}
    iv)   {f_{-m1, m2}}
    v)    {f_{ m2, m1}}
    vi)   {f_{-m2,-m1}}
    vii)  {f_{ m2,-m1}}
    viii) {f_{-m2, m1}}

    If assumptions are made regarding the original input signal,
    e.g. it's strictly real, then one may take advantage of
    symmetries of the big D wigners (i.e. function of all 3
    parameters alpha, beta, gamma) and so simplify the for-loops
    some and hence increase the speed of the program. However,
    the for-loops to follow will make no such assumptions.
    Whether the signal is real or complex, these for-loops will
    handle it.

    Fasten your seatbelt, folks. It's going to be a bumpy ride.

  */


  /* NOTE that I'm using the rdata, idata arrays as tmp space
     in the early going of the function */


  /***************************/
  /*                         */
  /* {f_{0,0}} coefficient   */
  /*                         */
  /***************************/
 
     
  /* compute the wigners I'll need */
  genWigTrans_L2( 0, 0, bw,
		  sinPts, cosPts,
		  sinPts2, cosPts2,
		  wignersTrans, scratch ) ;
  
  /* now, get the locations of where the
     samples I have to transform are, and
     where the coefficients have to go */
  
  sampHere = sampLoc_so3( 0, 0, bw ) ;
  coefHere = coefLoc_so3( 0, 0, bw ) ;
  
  /* ok, reset sample, coef ptrs */
  rcoeffsPtr = rcoeffs ;
  icoeffsPtr = icoeffs ;
  rdataPtr = rdata ;
  idataPtr = idata ;
  
  /* now advance by the computed amounts */
  rdataPtr += sampHere ;
  idataPtr += sampHere ;
  rcoeffsPtr += coefHere ;
  icoeffsPtr += coefHere ;
  
  /* now transform the real and imaginary parts
     of the data */
  
  
  wigNaiveSynthesis( 0, 0, bw, rcoeffsPtr,
		     wignersTrans, rdataPtr,
		     scratch ) ;
  
  wigNaiveSynthesis( 0, 0, bw, icoeffsPtr,
		     wignersTrans, idataPtr,
		     scratch ) ;


  /*** 0 <= m1 <= bw-1 ***/
  for ( m1 = 1 ; m1 < bw ; m1 ++ )
    {
      
      /* compute the wigners I'll need */
      genWigTrans_L2( m1, m1, bw,
		      sinPts, cosPts,
		      sinPts2, cosPts2,
		      wignersTrans, scratch ) ;
      
      /***************************/
      /*                         */
      /* {f_{m1,m1}} coefficient */
      /*                         */
      /***************************/


      /* now, get the locations of where the
	 samples I have to transform are, and
	 where the coefficients have to go */

      sampHere = sampLoc_so3( m1, m1, bw ) ;
      coefHere = coefLoc_so3( m1, m1, bw ) ;

      /* ok, reset sample, coef ptrs */
      rcoeffsPtr = rcoeffs ;
      icoeffsPtr = icoeffs ;
      rdataPtr = rdata ;
      idataPtr = idata ;
      
      /* now advance by the computed amounts */
      rdataPtr += sampHere ;
      idataPtr += sampHere ;
      rcoeffsPtr += coefHere ;
      icoeffsPtr += coefHere ;

      /* now transform the real and imaginary parts
	 of the data */


      wigNaiveSynthesis( m1, m1, bw, rcoeffsPtr,
			 wignersTrans, rdataPtr,
			 scratch ) ;
      
      wigNaiveSynthesis( m1, m1, bw, icoeffsPtr,
			 wignersTrans, idataPtr,
			 scratch ) ;
     
      /*****************************/
      /*                           */
      /* {f_{-m1,-m1}} coefficient */
      /*                           */
      /*****************************/

      if ( flag == 0 ) /* if data is complex */
	{

	  /* now, get the locations of where the
	     samples I have to transform are, and
	     where the coefficients have to go */

	  sampHere = sampLoc_so3( -m1, -m1, bw ) ;
	  coefHere = coefLoc_so3( -m1, -m1, bw ) ;

	  /* ok, reset sample, coef ptrs */
	  rcoeffsPtr = rcoeffs ;
	  icoeffsPtr = icoeffs ;
	  rdataPtr = rdata ;
	  idataPtr = idata ;
      
	  /* now advance by the computed amounts */
	  rdataPtr += sampHere ;
	  idataPtr += sampHere ;
	  rcoeffsPtr += coefHere ;
	  icoeffsPtr += coefHere ;

	  /* now transform the real and imaginary parts
	     of the data */


	  wigNaiveSynthesis( -m1, -m1, bw, rcoeffsPtr,
			     wignersTrans, rdataPtr,
			     scratch ) ;
      
	  wigNaiveSynthesis( -m1, -m1, bw, icoeffsPtr,
			     wignersTrans, idataPtr,
			     scratch ) ;

	}
      else  /* otherwise, use symmetry */
	{
	  sampHere = sampLoc_so3( m1, m1, bw );
	  sampHere2 = sampLoc_so3( -m1, -m1, bw );
	  for ( j = 0 ; j < 2*bw ; j ++ )
	    {
	      rdata[sampHere2+j] = rdata[sampHere+j];
	      idata[sampHere2+j] = -idata[sampHere+j];
	    }

	}

      

      /*****************************/
      /*                           */
      /* {f_{-m1,m1}} coefficient  */
      /*                           */
      /*****************************/


      /* now, get the locations of where the
	 samples I have to transform are, and
	 where the coefficients have to go */

      sampHere = sampLoc_so3( -m1, m1, bw ) ;
      coefHere = coefLoc_so3( -m1, m1, bw ) ;

      /* ok, reset sample, coef ptrs */
      rcoeffsPtr = rcoeffs ;
      icoeffsPtr = icoeffs ;
      rdataPtr = rdata ;
      idataPtr = idata ;
      
      /* now advance by the computed amounts */
      rdataPtr += sampHere ;
      idataPtr += sampHere ;
      rcoeffsPtr += coefHere ;
      icoeffsPtr += coefHere ;

      /* now transform the real and imaginary parts
	 of the data */


      wigNaiveSynthesis_symY( -m1, m1, bw, rcoeffsPtr,
			      wignersTrans, rdataPtr,
			      scratch ) ;
      
      wigNaiveSynthesis_symY( -m1, m1, bw, icoeffsPtr,
			      wignersTrans, idataPtr,
			      scratch ) ;
      


      /*****************************/
      /*                           */
      /* {f_{m1,-m1}} coefficient  */
      /*                           */
      /*****************************/

      if ( flag == 0 )  /* if data is complex */
	{

	  /* now, get the locations of where the
	     samples I have to transform are, and
	     where the coefficients have to go */

	  sampHere = sampLoc_so3( m1, -m1, bw ) ;
	  coefHere = coefLoc_so3( m1, -m1, bw ) ;

	  /* ok, reset sample, coef ptrs */
	  rcoeffsPtr = rcoeffs ;
	  icoeffsPtr = icoeffs ;
	  rdataPtr = rdata ;
	  idataPtr = idata ;
      
	  /* now advance by the computed amounts */
	  rdataPtr += sampHere ;
	  idataPtr += sampHere ;
	  rcoeffsPtr += coefHere ;
	  icoeffsPtr += coefHere ;

	  /* now transform the real and imaginary parts
	     of the data */


	  wigNaiveSynthesis_symY( m1, -m1, bw, rcoeffsPtr,
				  wignersTrans, rdataPtr,
				  scratch ) ;
      
	  wigNaiveSynthesis_symY( m1, -m1, bw, icoeffsPtr,
				  wignersTrans, idataPtr,
				  scratch ) ;

	}
      else  /* otherwise, use symmetry */
	{
	  sampHere = sampLoc_so3( -m1, m1, bw );
	  sampHere2 = sampLoc_so3( m1, -m1, bw );
	  for ( j = 0 ; j < 2*bw ; j ++ )
	    {
	      rdata[sampHere2+j] = rdata[sampHere+j];
	      idata[sampHere2+j] = -idata[sampHere+j];
	    }
	  
	}
 
    }

  /*** for 1 <= m1 <= bw-1 ***/
  for ( m1 = 1 ; m1 < bw ; m1 ++ )
    {
      
      /* compute the wigners I'll need */
      genWigTrans_L2( m1, 0, bw,
		      sinPts, cosPts,
		      sinPts2, cosPts2,
		      wignersTrans, scratch ) ;
      

      /***************************/
      /*                         */
      /* {f_{m1,0}} coefficient */
      /*                         */
      /***************************/


      /* get the locations of where the
	 samples I have to transform are, and
	 where the coefficients have to go */

      sampHere = sampLoc_so3( m1, 0, bw ) ;
      coefHere = coefLoc_so3( m1, 0, bw ) ;

      /* ok, reset sample, coef ptrs */
      rcoeffsPtr = rcoeffs ;
      icoeffsPtr = icoeffs ;
      rdataPtr = rdata ;
      idataPtr = idata ;
      
      /* now advance by the computed amounts */
      rdataPtr += sampHere ;
      idataPtr += sampHere ;
      rcoeffsPtr += coefHere ;
      icoeffsPtr += coefHere ;

      /* now transform the real and imaginary parts
	 of the data */


      wigNaiveSynthesis( m1, 0, bw, rcoeffsPtr,
			 wignersTrans, rdataPtr,
			 scratch ) ;
      
      wigNaiveSynthesis( m1, 0, bw, icoeffsPtr,
			 wignersTrans, idataPtr,
			 scratch ) ;



      /***************************/
      /*                         */
      /* {f_{-m1,0}} coefficient */
      /*                         */
      /***************************/


      if ( flag == 0 ) /* if data is complex */
	{

      
	  /* get the locations of where the
	     samples I have to transform are, and
	     where the coefficients have to go */

	  sampHere = sampLoc_so3( -m1, 0, bw ) ;
	  coefHere = coefLoc_so3( -m1, 0, bw ) ;

	  /* ok, reset sample, coef ptrs */
	  rcoeffsPtr = rcoeffs ;
	  icoeffsPtr = icoeffs ;
	  rdataPtr = rdata ;
	  idataPtr = idata ;
      
	  /* now advance by the computed amounts */
	  rdataPtr += sampHere ;
	  idataPtr += sampHere ;
	  rcoeffsPtr += coefHere ;
	  icoeffsPtr += coefHere ;

	  /* now transform the real and imaginary parts
	     of the data */

	  wigNaiveSynthesis_symX( -m1, 0, bw, rcoeffsPtr,
				  wignersTrans, rdataPtr,
				  scratch ) ;
      
	  wigNaiveSynthesis_symX( -m1, 0, bw, icoeffsPtr,
				  wignersTrans, idataPtr,
				  scratch ) ;


	}
      else  /* otherwise, use symmetry */
	{
	  sampHere = sampLoc_so3( m1, 0, bw );
	  sampHere2 = sampLoc_so3( -m1, 0, bw );

	  for ( j = 0 ; j < 2*bw ; j ++ )
	    {
	      rdata[sampHere2+j] = rdata[sampHere+j];
	      idata[sampHere2+j] = -idata[sampHere+j];
	    }
	  
	}

      /***************************/
      /*                         */
      /* {f_{0,m1}} coefficient  */
      /*                         */
      /***************************/

      
      /* get the locations of where the
	 samples I have to transform are, and
	 where the coefficients have to go */

      sampHere = sampLoc_so3( 0, m1, bw ) ;
      coefHere = coefLoc_so3( 0, m1, bw ) ;

      /* ok, reset sample, coef ptrs */
      rcoeffsPtr = rcoeffs ;
      icoeffsPtr = icoeffs ;
      rdataPtr = rdata ;
      idataPtr = idata ;
      
      /* now advance by the computed amounts */
      rdataPtr += sampHere ;
      idataPtr += sampHere ;
      rcoeffsPtr += coefHere ;
      icoeffsPtr += coefHere ;

      /* now transform the real and imaginary parts
	 of the data */


      wigNaiveSynthesis_symX( 0, m1, bw, rcoeffsPtr,
			      wignersTrans, rdataPtr,
			      scratch ) ;
      
      wigNaiveSynthesis_symX( 0, m1, bw, icoeffsPtr,
			      wignersTrans, idataPtr,
			      scratch ) ;
      
      /***************************/
      /*                         */
      /* {f_{0,-m1}} coefficient */
      /*                         */
      /***************************/
      
      if ( flag == 0 ) /* if data is complex */
	{
	  

	  /* get the locations of where the
	     samples I have to transform are, and
	     where the coefficients have to go */

	  sampHere = sampLoc_so3( 0, -m1, bw ) ;
	  coefHere = coefLoc_so3( 0, -m1, bw ) ;

	  /* ok, reset sample, coef ptrs */
	  rcoeffsPtr = rcoeffs ;
	  icoeffsPtr = icoeffs ;
	  rdataPtr = rdata ;
	  idataPtr = idata ;
      
	  /* now advance by the computed amounts */
	  rdataPtr += sampHere ;
	  idataPtr += sampHere ;
	  rcoeffsPtr += coefHere ;
	  icoeffsPtr += coefHere ;

	  /* now transform the real and imaginary parts
	     of the data */


	  wigNaiveSynthesis( 0, -m1, bw, rcoeffsPtr,
			     wignersTrans, rdataPtr,
			     scratch ) ;
      
	  wigNaiveSynthesis( 0, -m1, bw, icoeffsPtr,
			     wignersTrans, idataPtr,
			     scratch ) ;

	}
      else  /* otherwise, use symmetry */
	{
	  sampHere = sampLoc_so3( 0, m1, bw );
	  sampHere2 = sampLoc_so3( 0, -m1, bw );

	  for ( j = 0 ; j < 2*bw ; j ++ )
	    {
	      rdata[sampHere2+j] = rdata[sampHere+j];
	      idata[sampHere2+j] = -idata[sampHere+j];
	    }
	  
	}

    }


  /***
      1 <= m1 <= bw-1
      m1+1 <= m2 <= bw-1
  ***/

  for ( m1 = 1 ; m1 < bw ; m1 ++ )
    {
      for ( m2 = m1 + 1 ; m2 < bw ; m2 ++ )
	{

	  
	  /* compute the wigners I'll need */
	  genWigTrans_L2( m1, m2, bw,
			  sinPts, cosPts,
			  sinPts2, cosPts2,
			  wignersTrans, scratch ) ;


	  /***************************/
	  /*                         */
	  /* {f_{m1,m2}} coefficient */
	  /*                         */
	  /***************************/

	  /* get the locations of where the
	     samples I have to transform are, and
	     where the coefficients have to go */
	  
	  sampHere = sampLoc_so3( m1, m2, bw ) ;
	  coefHere = coefLoc_so3( m1, m2, bw ) ;
	  
	  /* ok, reset sample, coef ptrs */
	  rcoeffsPtr = rcoeffs ;
	  icoeffsPtr = icoeffs ;
	  rdataPtr = rdata ;
	  idataPtr = idata ;
	  
	  /* now advance by the computed amounts */
	  rdataPtr += sampHere ;
	  idataPtr += sampHere ;
	  rcoeffsPtr += coefHere ;
	  icoeffsPtr += coefHere ;
	  
	  /* now transform the real and imaginary parts
	     of the data */
	  
	  wigNaiveSynthesis( m1, m2, bw, rcoeffsPtr,
			     wignersTrans, rdataPtr,
			     scratch ) ;
	  
	  wigNaiveSynthesis( m1, m2, bw, icoeffsPtr,
			     wignersTrans, idataPtr,
			     scratch ) ;
	  

	  /*****************************/
	  /*                           */
	  /* {f_{-m1,-m2}} coefficient */
	  /*                           */
	  /*****************************/

	  if ( flag == 0 ) /* if data is complex */
	    {


	      /* get the locations of where the
		 samples I have to transform are, and
		 where the coefficients have to go */
	  
	      sampHere = sampLoc_so3( -m1, -m2, bw ) ;
	      coefHere = coefLoc_so3( -m1, -m2, bw ) ;
	  
	      /* ok, reset sample, coef ptrs */
	      rcoeffsPtr = rcoeffs ;
	      icoeffsPtr = icoeffs ;
	      rdataPtr = rdata ;
	      idataPtr = idata ;


	      /* now advance by the computed amounts */
	      rdataPtr += sampHere ;
	      idataPtr += sampHere ;
	      rcoeffsPtr += coefHere ;
	      icoeffsPtr += coefHere ;
	  
	      /* now transform the real and imaginary parts
		 of the data */

	      wigNaiveSynthesis_symX( -m1, -m2, bw, rcoeffsPtr,
				      wignersTrans, rdataPtr,
				      scratch ) ;
	  
	      wigNaiveSynthesis_symX( -m1, -m2, bw, icoeffsPtr,
				      wignersTrans, idataPtr,
				      scratch ) ;


	    }
	  else  /* otherwise, use symmetry */
	    {
	      sampHere = sampLoc_so3( m1, m2, bw );
	      sampHere2 = sampLoc_so3( -m1, -m2, bw );

	      for ( j = 0 ; j < 2*bw ; j ++ )
		{
		  rdata[sampHere2+j] = rdata[sampHere+j];
		  idata[sampHere2+j] = -idata[sampHere+j];
		}
	    }

	  /****************************/
	  /*                          */
	  /* {f_{m1,-m2}} coefficient */
	  /*                          */
	  /****************************/
  

	  /* get the locations of where the
	     samples I have to transform are, and
	     where the coefficients have to go */
	  
	  sampHere = sampLoc_so3( m1, -m2, bw ) ;
	  coefHere = coefLoc_so3( m1, -m2, bw ) ;


	  /* ok, reset sample, coef ptrs */
	  rcoeffsPtr = rcoeffs ;
	  icoeffsPtr = icoeffs ;
	  rdataPtr = rdata ;
	  idataPtr = idata ;
      
	  /* now advance by the computed amounts */
	  rdataPtr += sampHere ;
	  idataPtr += sampHere ;
	  rcoeffsPtr += coefHere ;
	  icoeffsPtr += coefHere ;
	  
	  /* now transform the real and imaginary parts
	     of the data */

	  wigNaiveSynthesis_symY( m1, -m2, bw, rcoeffsPtr,
				  wignersTrans, rdataPtr,
				  scratch ) ;
	  
	  wigNaiveSynthesis_symY( m1, -m2, bw, icoeffsPtr,
				  wignersTrans, idataPtr,
				  scratch ) ;

	  /*****************************/
	  /*                           */
	  /* {f_{-m1,m2}} coefficient  */
	  /*                           */
	  /*****************************/

	  if ( flag == 0 ) /* if data is complex */
	    {

	      /* get the locations of where the
		 samples I have to transform are, and
		 where the coefficients have to go */
	  
	      sampHere = sampLoc_so3( -m1, m2, bw ) ;
	      coefHere = coefLoc_so3( -m1, m2, bw ) ;


	      /* ok, reset sample, coef ptrs */
	      rcoeffsPtr = rcoeffs ;
	      icoeffsPtr = icoeffs ;
	      rdataPtr = rdata ;
	      idataPtr = idata ;
      
	      /* now advance by the computed amounts */
	      rdataPtr += sampHere ;
	      idataPtr += sampHere ;
	      rcoeffsPtr += coefHere ;
	      icoeffsPtr += coefHere ;
	  
	      /* now transform the real and imaginary parts
		 of the data */

	      wigNaiveSynthesis_symY( -m1, m2, bw, rcoeffsPtr,
				      wignersTrans, rdataPtr,
				      scratch ) ;
	  
	      wigNaiveSynthesis_symY( -m1, m2, bw, icoeffsPtr,
				      wignersTrans, idataPtr,
				      scratch ) ;


	    }
	  else  /* otherwise, use symmetry */
	    {
	      sampHere = sampLoc_so3( m1, -m2, bw );
	      sampHere2 = sampLoc_so3( -m1, m2, bw );

	      for ( j = 0 ; j < 2*bw ; j ++ )
		{
		  rdata[sampHere2+j] = rdata[sampHere+j];
		  idata[sampHere2+j] = -idata[sampHere+j];
		}

	    }
	  

	  /***************************/
	  /*                         */
	  /* {f_{m2,m1}} coefficient */
	  /*                         */
	  /***************************/
	  
	  /* get the locations of where the
	     samples I have to transform are, and
	     where the coefficients have to go */
	  
	  sampHere = sampLoc_so3( m2, m1, bw ) ;
	  coefHere = coefLoc_so3( m2, m1, bw ) ;
	  
	  /* ok, reset sample, coef ptrs */
	  rcoeffsPtr = rcoeffs ;
	  icoeffsPtr = icoeffs ;
	  rdataPtr = rdata ;
	  idataPtr = idata ;
	  
	  /* now advance by the computed amounts */
	  rdataPtr += sampHere ;
	  idataPtr += sampHere ;
	  rcoeffsPtr += coefHere ;
	  icoeffsPtr += coefHere ;
	  
	  /* now transform the real and imaginary parts
	     of the data */

	  wigNaiveSynthesis_symX( m2, m1, bw, rcoeffsPtr,
				  wignersTrans, rdataPtr,
				  scratch ) ;
	  
	  wigNaiveSynthesis_symX( m2, m1, bw, icoeffsPtr,
				  wignersTrans, idataPtr,
				  scratch ) ;



	  /*****************************/
	  /*                           */
	  /* {f_{-m2,-m1}} coefficient */
	  /*                           */
	  /*****************************/
	  
	  if ( flag == 0 ) /* if data is complex */
	    {

	      /* get the locations of where the
		 samples I have to transform are, and
		 where the coefficients have to go */
	  
	      sampHere = sampLoc_so3( -m2, -m1, bw ) ;
	      coefHere = coefLoc_so3( -m2, -m1, bw ) ;
	  
	      /* ok, reset sample, coef ptrs */
	      rcoeffsPtr = rcoeffs ;
	      icoeffsPtr = icoeffs ;
	      rdataPtr = rdata ;
	      idataPtr = idata ;
	  
	      /* now advance by the computed amounts */
	      rdataPtr += sampHere ;
	      idataPtr += sampHere ;
	      rcoeffsPtr += coefHere ;
	      icoeffsPtr += coefHere ;
	  
	      /* now transform the real and imaginary parts
		 of the data */
	  
	      wigNaiveSynthesis( -m2, -m1, bw, rcoeffsPtr,
				 wignersTrans, rdataPtr,
				 scratch ) ;
	  
	      wigNaiveSynthesis( -m2, -m1, bw, icoeffsPtr,
				 wignersTrans, idataPtr,
				 scratch ) ;


	    }
	  else  /* otherwise, use symmetry */
	    {
	      sampHere = sampLoc_so3( m2, m1, bw );
	      sampHere2 = sampLoc_so3( -m2, -m1, bw );

	      for ( j = 0 ; j < 2*bw ; j ++ )
		{
		  rdata[sampHere2+j] = rdata[sampHere+j];
		  idata[sampHere2+j] = -idata[sampHere+j];
		}
	    }
	  
	  
	  /****************************/
	  /*                          */
	  /* {f_{m2,-m1}} coefficient */
	  /*                          */
	  /****************************/
  

	  /* get the locations of where the
	     samples I have to transform are, and
	     where the coefficients have to go */
	  
	  sampHere = sampLoc_so3( m2, -m1, bw ) ;
	  coefHere = coefLoc_so3( m2, -m1, bw ) ;
	  
	  /* ok, reset sample, coef ptrs */
	  rcoeffsPtr = rcoeffs ;
	  icoeffsPtr = icoeffs ;
	  rdataPtr = rdata ;
	  idataPtr = idata ;
      
	  /* now advance by the computed amounts */
	  rdataPtr += sampHere ;
	  idataPtr += sampHere ;
	  rcoeffsPtr += coefHere ;
	  icoeffsPtr += coefHere ;
	  
	  /* now transform the real and imaginary parts
	     of the data */

	  wigNaiveSynthesis_symY( m1, -m2, bw, rcoeffsPtr,
				  wignersTrans, rdataPtr,
				  scratch ) ;
	  
	  wigNaiveSynthesis_symY( m1, -m2, bw, icoeffsPtr,
				  wignersTrans, idataPtr,
				  scratch ) ;
	  

	  /****************************/
	  /*                          */
	  /* {f_{-m2,m1}} coefficient */
	  /*                          */
	  /****************************/
  

	  if ( flag == 0 ) /* if data is complex */
	    {


	      /* get the locations of where the
		 samples I have to transform are, and
		 where the coefficients have to go */
	  
	      sampHere = sampLoc_so3( -m2, m1, bw ) ;
	      coefHere = coefLoc_so3( -m2, m1, bw ) ;
	  
	      /* ok, reset sample, coef ptrs */
	      rcoeffsPtr = rcoeffs ;
	      icoeffsPtr = icoeffs ;
	      rdataPtr = rdata ;
	      idataPtr = idata ;
      
	      /* now advance by the computed amounts */
	      rdataPtr += sampHere ;
	      idataPtr += sampHere ;
	      rcoeffsPtr += coefHere ;
	      icoeffsPtr += coefHere ;
	  
	      /* now transform the real and imaginary parts
		 of the data */

	      wigNaiveSynthesis_symY( -m1, m2, bw, rcoeffsPtr,
				      wignersTrans, rdataPtr,
				      scratch ) ;
	  
	      wigNaiveSynthesis_symY( -m1, m2, bw, icoeffsPtr,
				      wignersTrans, idataPtr,
				      scratch ) ;

	    }
	  else  /* otherwise, use symmetry */
	    {
	      sampHere = sampLoc_so3( m2, -m1, bw );
	      sampHere2 = sampLoc_so3( -m2, m1, bw );

	      for ( j = 0 ; j < 2*bw ; j ++ )
		{
		  rdata[sampHere2+j] = rdata[sampHere+j];
		  idata[sampHere2+j] = -idata[sampHere+j];
		}
	    }
		    
	}
    }

  /* I need to set some zeros in the rdata, idata  arrays,
     so that I can take the fft correctly */

  /* reset ptrs to correct starting positions */
  rdataPtr = rdata + (n)*(bw) ;
  idataPtr = idata + (n)*(bw) ;

  for ( m1 = 0 ; m1 < bw  ; m1 ++ )
    {
      memset( rdataPtr, 0, sizeof(double) * n );
      memset( idataPtr, 0, sizeof(double) * n );
      
      rdataPtr += (2*n)*(bw) ;
      idataPtr += (2*n)*(bw) ;
    }
  
  rdataPtr = rdata + bw*n*(n);
  idataPtr = idata + bw*n*(n);
  
  memset( rdataPtr, 0, sizeof(double) * n * n );
  memset( idataPtr, 0, sizeof(double) * n * n );
  
  rdataPtr += n * n + n*bw;
  idataPtr += n * n + n*bw;
  
  for ( m1 = 1 ; m1 < bw  ; m1 ++ )
    {
      memset( rdataPtr, 0, sizeof(double) * n );
      memset( idataPtr, 0, sizeof(double) * n );
      
      rdataPtr += (2*n)*(bw) ;
      idataPtr += (2*n)*(bw) ;
    }

  /*
    Stage 2: transpose! Note I'm using the rdata, idata arrays
    as tmp space
  */
  
  transpose_so3( rdata, t1r, n, n*n ) ;
  transpose_so3( idata, t1i, n, n*n ) ;


  /*
    Stage 3: FFT the "rows". Instead of treating the signal as
    3-D object, I can also think of it as an array of size
    (n^2) x n. This means all I'm doing in the first stage
    is taking n^2-many FFTs, each of length n.

    NOTE: Since I'm reusing the FFT code from SpharmonicKit,
    even though I'm doing the INVERSE SO(3) transform
    here, I need to call grid_fourier_so3  -> the signs
    on the complex exponentials are switched (detailed
    explanation to be put here eventually, but trust
    me)
  */
  
  grid_fourier_so3( t1r, t1i,
		    rdata, idata,
		    n*n, n,
		    scratch ) ;

  /* normalize the Fourier coefficients (sorry, have to do it) */
  /*
    dn = sqrt( (double) n );
    for ( j = 0 ; j < n*n*n; j++ )
    {
    rdata[ j ] *= dn ;
    idata[ j ] *= dn ;
    }
  */

  /*
    Stage 4: transpose! Note I'm using the rdata, idata arrays
    as tmp space
  */
  
  transpose_so3( rdata, t1r, n, n*n ) ;
  transpose_so3( idata, t1i, n, n*n ) ;


  /*
    Stage 5: FFT again. Note that THIS TIME, the rdata, idata
    arrays will hold the final answers I want
  */

  grid_fourier_so3( t1r, t1i,
		    rdata, idata,
		    n*n, n,
		    scratch ) ;

  /* normalize the Fourier coefficients (sorry, have to do it) */
  dn = (double) n ;
  dn *= ( ((double) bw) / M_PI ) ;
  for ( j = 0 ; j < n*n*n; j++ )
    {
      rdata[ j ] *= dn ;
      idata[ j ] *= dn ;
    }

  /* and that's all, folks */


}
