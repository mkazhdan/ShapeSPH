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
  functions having to do with rotating functions on the
  sphere by massaging their Fourier coefficients with
  the Wigner D functions

  NOTE: In the (vain?) attempt to conserve memory, the function

  rotateCoefAll_mem

  defined in this file will write over the input coefficients.
  If this bothers you, you can always use

  rotate_so3.{h,c}



*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "complex.h"
#include "primitive_FST.h"
#include "rotate_so3_mem.h"
#include "FST_semi_memo.h"

/*****************************************************************/
/*****************************************************************/
/*
  genExp_mem: given an angle phi, and an bandlimit bw, this routine
          will generate the values

	  Exp[-I k phi], where k = -(bw-1) ... -1 0 1 ... (bw-1)

	  This Exp represents the one found in the definition of the
	  Wigner-D's:

          D_{mm'}^l (alpha, beta, gamma) =
               exp(-I alpha m) d_{mm'}^l(beta) exp(-I gamma m')

	  It will be used to massage the spherical coefficients
	  of a function f in L^2(S^2) which is rotated, e.g.
	  in rotateCoefDegree. Note that it will generate all
	  the Exp values I'll possibly need to massage all the
	  coefficients ... no need to regenerate values for each
	  degree l (provided I know where in the array the k = 0
	  term is )

  arguments:

         INPUT:

	 phi: the angle

	 bw : the bandlimit

         OUTPUT:

	 expValsR, expValsI:
            Two REAL arrays, each of length 2*bw-1, will contain
	    the real and imaginary parts of the values in this order:

	    Exp[-I*phi*(-(bw-1))] ... Exp[0] ... Exp[-I*phi*(bw-1)]

*/

void genExp_mem( int bw ,
		 REAL phi ,
		 REAL *expValsR,
		 REAL *expValsI )
{
  int i, expZero ;

  expZero = bw - 1 ;

  for ( i = -expZero ; i < expZero + 1 ; i++ )
    {
      expValsR[expZero+i] = (REAL) cos(((REAL) i)*phi);
      expValsI[expZero+i] = (REAL) -sin(((REAL) i)*phi);
    }
}

/*****************************************************************

 wignerdmat_mem: given a wigner little-d matrix of degree L-1 evaluated
             at some angle beta, construct the wigner little-d
             matrix of degree L, EVALUATED AT THAT SAME BETA.

 arguments:

         INPUT:

         L: degree of wigner little-d's that you want

         matIn: array of REALs of length (2*L-1)*(2*L-1), containing
                the values of d_{m1,m2}^(L-1)(beta) for
                -(L-1) <= m1, m2 <= L-1

                So matIn contains the Wigner litte-d's for the
                PREVIOUS degree, L - 1.

         trigs: array of REALS of length 2, containing, IN THIS ORDER,
                cos(beta/2) and sin(beta/2)

         sqrts: array of REALS of length >= L+1. It contains the square roots
                of 0 through L+1

         workspace: scratch array of REALS of size (2*L)^2


         OUTPUT:

         matOut: array of REALS of length (2*L+1)*(2*L+1), will contain
                 the values of d_{m1,m2}^(L)(beta) for
                 -(L+1) <= m1, m2 <= (L+1)


         The 1-d matIn/matOut arrays are in ROW-MAJOR format, with m1 indexing the
         rows, and m2 indexing the columns of the 2-D matrix. The indexing
         goes from -L to L top to bottom/left to right. E.g. For the
         case of L = 1, the 2-d array is

         ( d_{-1,-1}^1   d_{-1,0}^1   d_{-1,1}^1 )
         (                                       )
         ( d_{ 0,-1}^1   d_{ 0,0}^1   d_{ 0,1}^1 )
         (                                       )
         ( d_{ 1,-1}^1   d_{ 1,0}^1   d_{ 1,1}^1 )

         and it's stored in the 1-d array in this manner:

         ( d_{-1,-1}^1   d_{-1,0}^1 ... d_{ 1,0}^1  d_{ 1,1}^1 )


 NOTE: 1) These wigner little-ds are NOT normalized! They will be used only
          for rotations.

       2) If L = 0 or 1, then matIn can be any ol' array (since it isn't used
          in this case). The array matOut will be what you want. The L = 0 and 1
          arrays are constructed by hand.

       3) This algorithm is a slight variation of the one given by T. Risbo in

          "Fourier transform summation of Legendre series and D-Functions"
          T. Risbo
          Journal of Geodesy
          1996
          volume 70: p. 383 - 396.

*****************************************************************/

void wignerdmat_mem( int L,
		     REAL *matIn,
		     REAL *matOut,
		     REAL *trigs,
		     REAL *sqrts,
		     REAL *workspace )
{

  int i, j ;
  int deg, tmpdim ;
  REAL cosVal, sinVal ;
  REAL *ptrIn, *ptrOut, *ptrTmp ;
  REAL rdeg ;

  cosVal = trigs[0] ;
  sinVal = trigs[1] ;

  ptrIn = matIn ;
  ptrOut = matOut ;
  ptrTmp = workspace ;

  if ( L > 1 )
    {
      /* zero out tmp and final arrays */
      memset( ptrTmp, 0, sizeof(REAL) * (2*L)*(2*L) );
      memset( matOut, 0, sizeof(REAL) * (2*L+1)*(2*L+1) );

      memcpy( ptrTmp, matIn, sizeof(REAL) * (2*L-1)*(2*L-1) ); 

      for ( deg = 2*L-1 ; deg < 2*L+1 ; deg++ )
	{
	  rdeg = 1./(( REAL ) deg) ;

	  tmpdim = deg + 1 ;
	  memset( matOut, 0, sizeof(REAL) * (tmpdim)*(tmpdim) );

	  for ( i = 0 ; i < deg ; i++ )
	    for ( j = 0 ; j < deg ; j++ )
	      {
		matOut[(i*tmpdim)+j] +=
		  rdeg*sqrts[deg-i]*sqrts[deg-j]*
		  ptrTmp[(i*deg)+j]*cosVal ;

		matOut[((i+1)*tmpdim)+j] -=
		  rdeg*sqrts[i+1]*sqrts[deg-j]*
		  ptrTmp[(i*deg)+j]*sinVal ;

		matOut[(i*tmpdim)+(j+1)] +=
		  rdeg*sqrts[deg-i]*sqrts[j+1]*
		  ptrTmp[(i*deg)+j]*sinVal ;

		matOut[((i+1)*tmpdim)+(j+1)] +=
		  rdeg*sqrts[i+1]*sqrts[j+1]*
		  ptrTmp[(i*deg)+j]*cosVal ;

	      }
	  if ( deg == 2*L - 1 )
	    memcpy( ptrTmp, matOut, sizeof(REAL) * (tmpdim)*(tmpdim) ); 

	}
    }
  else if ( L == 0 )
    {
      matOut[0] = 1 ;
    }
  else   /* L == 1 */
    {
      matOut[0] = cosVal*cosVal ;
      matOut[1] = sqrts[2]*cosVal*sinVal ;
      matOut[2] = sinVal*sinVal ;
      
      matOut[3] = -matOut[1] ;
      matOut[4] = matOut[0]-matOut[2];
      matOut[5] = matOut[1] ;

      matOut[6] = matOut[2] ;
      matOut[7] = -matOut[1] ;
      matOut[8] = matOut[0] ;
    }

}

/*****************************************************************/
/*****************************************************************/
/*
  rotateCoefDegree_mem: given the degree L spherical coefficients {f_L^m},
                    where -L <= m <= L of some function on the sphere,
		    and three Euler angles alpha, beta, gamma, multiply
		    those coefficients by the wigner D-mat(alpha,beta,gamma).
		    In other words, I'm rotating f by the Euler angles,
		    and I'm massaging the coefficients as representation
		    theory says I should.
		    
  Arguments:

	L: degree of coefficients (i.e. the representation)

	cfR, cfI: two 1-D REAL arrays, each of length 2L+1,
	      containing the real and imaginary parts of
	      the coefficients of the signal to be rotated.
	      Each is in the following order:

	      [ f_L^(-L) f_L^(-L+1) ... f_L^(-1) 0 f_L^(1) ... f_L^(L) ]

	expAR, expAI: two 1-D REAL arrays, each of length at least
	        2L+1, containing the real and imaginary parts of the
		following array in the following order:

		[exp(-I*alpha*(-M)) ... exp(-I*alpha*M)] 

		where M >= L.

		I saw "at least" because, since all the coefficients
		will be massaged, of all degrees, you might as well
		precompute all the numbers you'll possibly need:

		[exp(-I*alpha*(-(bw-1))) ... exp(-I*alpha*(bw-1))]
		
		at the start, and be done with it ... no need to compute
		the same stuff again and again, for each L

	expGR, expGI: two 1-D REAL arrays, each of length at least
	        2L+1, containing the real and imaginary parts of the
		following array in the following order:

		[exp(-I*gamma*(-M)) ... exp(-I*gamma*M)] 

		where M >= L.

		See the explanation for expA. It's the same case here.

        expZero -> WHERE in expA, expG does the 0 term, i.e. exp(-I*0),
                occur? I need to know where to start, especially since
		expA, expG might be arrays longer then I need for the
		particular degree L this routine is given.


	dmat -> REAL array of length (2*L+1)^2: it contains the
	        "beta" array generated by wignerdmat, i.e. the array
		you get when plugging beta into that routine.

	OUTPUT:

	cfrR, cfrI: two 1-D REAL arrays, each of length 2*L+1,
	       containing the real and imaginary parts of the
	       degree L coefficients of the rotated function



  Basically, this routine does a matrix-vector multiply, where the
  matrix is [Wigner-D_{M,M'}(alpha,beta,gamma)] of dimension (2L+1)x(2L+1),
  and the vector is the vector of spherical coefficients f_l^m. To
  be a little more precise ... the rows are indexed by the order M,
  and the columns by M', where -L <= M, M' <= L. E.g.

  [ D_{-1,-1}  D_{-1,0}  D_{-1,1} ]   [ f_1^{-1} ]
  [                               ]   [          ]
  [ D_{ 0,-1}  D_{ 0,0}  D_{ 0,1} ] . [ f_1^{ 0} ]
  [                               ]   [          ]
  [ D_{ 1,-1}  D_{ 1,0}  D_{ 1,1} ]   [ f_1^{ 1} ]

  NOTE that these Wigner-D's are not normalized like they are when doing
  the forward/inverse SO(3) transform!


*/

void rotateCoefDegree_mem( int L,
			   REAL *cfR, REAL *cfI,
			   REAL *cfrR, REAL *cfrI,
			   REAL *expAR, REAL *expAI,
			   REAL *expGR, REAL *expGI,
			   int expZero,
			   REAL *dmat )
{
  int i, j, dim ;
  REAL *expAPtrR, *expAPtrI ;
  REAL *expGPtrR, *expGPtrI ;
  REAL tmpAR, tmpAI, tmpR, tmpI ;
  REAL tmp2R, tmp2I, tmp3R, tmp3I ;


  expAPtrR = &expAR[ expZero - L ];
  expAPtrI = &expAI[ expZero - L ];

  expGPtrR = &expGR[ expZero - L ];
  expGPtrI = &expGI[ expZero - L ];

  dim = 2*L + 1 ;

  for ( i = 0 ; i < dim ; i ++ )
    {
      tmpAR = expAPtrR[ i ] ;
      tmpAI = expAPtrI[ i ] ;

      tmpR = (REAL) 0 ;
      tmpI = (REAL) 0 ;

      for ( j = 0 ; j < dim ; j ++ )
	{
	  tmp2R = dmat[i*dim+j]*expGPtrR[j] ;
	  tmp2I = dmat[i*dim+j]*expGPtrI[j] ;

	  tmp3R = tmp2R*tmpAR - tmp2I*tmpAI ;
	  tmp3I = tmp2R*tmpAI + tmp2I*tmpAR ;

	  tmp2R = tmp3R*cfR[j] - tmp3I*cfI[j];
	  tmp2I = tmp3R*cfI[j] + tmp3I*cfR[j];

	  tmpR += tmp2R ;
	  tmpI += tmp2I ;
	}

      cfrR[i] = tmpR ;
      cfrI[i] = tmpI ;
    }

}

/*****************************************************************/
/*****************************************************************/
/*
  rotateCoefAll_mem: given the spherical harmonic coefficients of
                 a bandlimited function f defined on the sphere,
		 massage *all* the coefficients in such a way as
		 to get the coefficients of a *rotated* version
		 of that function.

		 More plainly, let g be an element of SO(3). We can
		 write it in terms of the Euler angles alpha, beta,
		 and gamma, where alpha and gamma are the angles of
		 rotation about the z-axis, and beta is rotation
		 about the y-aixs. So, given a bandlimited function
		 f defined on the sphere, we can rotate it:

		 f(omega) -> \Lambda(g)f(omega) = f(g^(-1) omega).

		 What rotateCoefAll does is massage all the coefficients
		 of f in such a way as to get coefficients of
		 \Lambda(g)f(omega). Representation theory says
		 we can, and so we will.

		 NOTE: Here are order of rotation events:
		       1) rotate by gamma about the z-axis
		       2) rotate by beta about the y-axis
		       3) rotate by alpha about the z-axis.
		 

		 NOTE: This routine will WRITE OVER the input
		       coefficients.


  arguments:

           bw: bandlimit of the function being rotated, i.e.
	       what's the bandlimit of the function you want
	       to rotate?

	   degOut: Now, instead of rotating the coefficients of
	           all the degrees 0 <= L < bw, for whatever
		   reason you may want to go only "so far" and
		   rotate only through L = degOut. The rotated
		   coefficients for degrees degOut+1 <= L < bwIn
		   will be set to 0.

	   alpha, beta, gamma: the Euler angles

	   cfR, cfI: two REAL arrays containing the real and imaginary
	          parts of the spherical coefficients of f BEFORE AND
		  AFTER MASSAGING, each of size (bw)^2

	   scratch: REAL scratch space of size
	            2*bw + 2*(2*bw-1)^2 + (2*degOut)^2 +
		    2*4*(2*bw-1)


  NOTE:	The coefficients in cf will be in the same order
        as that produced by a forward spherical transform
	routine found in SpharmonicKit. This will require
	the use of the indexing functions
	
	seanindex
	seanindex2
	
	defined in primitive_FST.[c,h].

  NOTE: Since the routines in SpharmonicKit produce two separate double
        arrays, for the real and imaginary parts of the coefficients,
	they will need to be interleaved prior to entering rotateCoefAll.
	Sorry.

*/

void rotateCoefAll_mem( int bw,
			int degOut,
			REAL alpha,
			REAL beta,
			REAL gamma,
			REAL *cfR, REAL *cfI,
			REAL *scratch )
{
  int i, m ;
  int deg, dim, dummy ;
  REAL *sqrts, trigs[2] ;
  REAL *matIn, *matOut, *workspace ;
  REAL *cfTmpR, *cfTmpI ;
  REAL *cfrTmpR, *cfrTmpI ;
  REAL *expAR, *expAI, *expGR, *expGI ;


  sqrts = scratch ;
  matIn = sqrts + 2*bw ;
  matOut = matIn + (2*bw-1)*(2*bw-1);

  cfTmpR = matOut + (2*bw-1)*(2*bw-1);
  cfTmpI = cfTmpR + (2*bw - 1);
  cfrTmpR = cfTmpI + (2*bw - 1);
  cfrTmpI = cfrTmpR + (2*bw - 1);
  expAR = cfrTmpI + (2*bw - 1);
  expAI = expAR + (2*bw - 1);
  expGR = expAI + (2*bw - 1);
  expGI = expGR + (2*bw - 1);
  workspace = expGI + (2*bw - 1) ;
  /* workspace still needs (2*degOut)^2 space */

  /*
    We'll first precompute stuff needed for the
    wignerdmat routine: the square roots, and the trig values
  */

  for ( i = 0 ; i < 2*bw ; i ++ )
    sqrts[i] = (REAL) sqrt( (REAL) i ) ;

  trigs[0] = (REAL) cos(0.5*beta);
  trigs[1] = (REAL) sin(0.5*beta);
 
  /*
    For the rotateCoefDegree routine, we'll need the values
    of the complex exponentials
  */

  genExp_mem( bw, alpha, expAR, expAI ) ;
  genExp_mem( bw, gamma, expGR, expGI ) ;


  /*
    Now, to massage the coefficients, one degree at a time

    Note that this will be done for degrees L = 0 through
    L = degOut (provided that also L < bwOut). The remaining
    degrees will be set to 0.
  */

  for ( deg = 0 ; deg <= degOut ; deg ++ )
    {
      dim = 2*deg - 1 ;

      if ( deg > 0 )
	memcpy( matIn, matOut, sizeof(REAL)*(dim*dim) );

      /* grab the spherical coefficients at this degree */
      for ( m = -deg ; m < deg+1 ; m++ )
	{
	  dummy = seanindex(m, deg, bw);
	  cfTmpR[m+deg] = cfR[dummy];
	  cfTmpI[m+deg] = cfI[dummy];
	}

      /* generate the required wigner-d mat */
      wignerdmat_mem( deg, matIn, matOut, trigs, sqrts, workspace );

      /* massage the coefficients at this degree */
      rotateCoefDegree_mem( deg,
			    cfTmpR, cfTmpI,
			    cfrTmpR, cfrTmpI,
			    expAR, expAI,
			    expGR, expGI,
			    bw - 1,
			    matOut ) ;

      /* save the rotated coefficients */
      for ( m = -deg ; m < deg+1 ; m++ )
	{
	  dummy = seanindex(m, deg, bw);
	  cfR[dummy] = cfrTmpR[m+deg];
	  cfI[dummy] = cfrTmpI[m+deg];
	}

    }

  /* zero out remaining coefficients */
  for ( deg = degOut + 1 ; deg < bw ; deg ++ )
    for ( m = -deg ; m < deg + 1 ; m ++ )
      {
	dummy = seanindex(m, deg, bw);
	cfR[dummy] = (REAL) 0;
	cfI[dummy] = (REAL) 0;
      }
}

/*****************************************************************/
/*****************************************************************/
/*
  rotateFct_mem: given a function defined on the sphere,
             this routine will rotate it. The rotation
	     will be described of the Euler angles
	     alpha, beta, gamma, with
	     
	     0 <= alpha, gamma < 2*pi
	     0 <= beta <= pi
	     
	     Here are order of rotation events:
		   1) rotate by gamma about the z-axis
		   2) rotate by beta about the y-axis
		   3) rotate by alpha about the z-axis.

  arguments:

           bw: bandlimit of function


	   degOut: Now, instead of rotating the coefficients of
	           all the degrees 0 <= L < bwIn, for whatever
		   reason you may want to go only "so far" and
		   rotate only through L = degOut. The rotated
		   coefficients for degrees degOut+1 <= L < bwIn
		   will be considered to be 0.

           sigR: REAL array of the *real* parts of the signal sampled
	           on the usual 2bwIn * 2bwIn grid on the sphere;
		 WILL BE WRITTEN OVER WITH THE REAL PARTS OF THE SAMPLE
		 POINTS OF THE ROTATED SIGNAL

	   sigI: REAL array of the *imaginary* parts of the signal
	           sampled on the usual 2bwIn * 2bwIn grid on the sphere;
		 WILL BE WRITTEN OVER WITH THE IMAGINARY PARTS OF THE
		 SAMPLE POINTS OF THE ROTATED SIGNAL

           sigOutR: REAL array of the *real* parts of the rotated signal
	            sampled on the usual 2bwOut * 2bwOut grid on the sphere

	   sigOutI: REAL array of the *imaginary* parts of the rotated
	            signal sampled on the usual 2bwOut * 2bwOut grid on
		    the sphere

	   alpha, beta, gamma: the Euler angles defining the rotation

	   scratch: REAL array of size (gulp)
		    10*bwIn^2 + 48*bwIn - 8

	   spharmonic_pml_table: 
                  should be a (double **) pointer to
		  the result of a call to Spharmonic_Pml_Table. Because this
		  table is re-used in the inverse transform, and because for
		  timing purposes the computation of the table is not included,
		  it is passed in as an argument.
		  
		  spharmonic_pml_table will be an array of (double *) pointers
		  the array being of length TableSize(m,bwIn)

	   transpose_spharmonic_pml_table:
	          should be the (double **) result of a call
		  to Transpose_Spharmonic_Pml_Table()


  NOTE: Since this routine requires huge amounts of memory, it might
        be more convenient if, instead of using rotateFct, you instead
	use rotateCoefAll on its own, e.g.

	1) use your favourite S^2 program to get the coefficients
	   of the function, and write them to disk

	2) in a separate program, read in the coefficients and then
	   call rotateCoefAll, then write out the massaged coefficients
	   to disk

	3) use your favourite S^2 program to take the inverse transform
	   of the massaged coefficients and hence obtain the rotated
	   function.

	Most of the memory in this function is devoted to doing the forward
	and inverse S^2 transforms, e.g. arrays to contain the precomputed
	DCTs of the Legendre functions.



*/

void rotateFct_mem( int bw, int degOut,
		REAL *sigR, REAL *sigI,
		REAL alpha, REAL beta, REAL gamma,
		REAL *scratch,
		REAL **spharmonic_pml_table,
		REAL **transpose_spharmonic_pml_table )
{
  REAL *rcoeffs, *icoeffs, *workspace ;

  rcoeffs = scratch ;
  icoeffs = rcoeffs + (bw*bw) ;
  workspace = icoeffs + (bw*bw) ;

  /* compute spherical coefficients of input signal */
  FST_semi_memo( sigR, sigI,
		 rcoeffs, icoeffs,
		 2*bw, spharmonic_pml_table,
		 workspace, 1, bw ) ;

  /* now massage the coefficients */
  /* note that I'm using the workspace array again */
  rotateCoefAll_mem( bw, degOut,
		 alpha, beta, gamma,
		 rcoeffs, icoeffs,
		 workspace ) ;

  /* take inverse spherical transform */
  InvFST_semi_memo( rcoeffs, icoeffs,
		    sigR, sigI,
		    2*bw,
		    transpose_spharmonic_pml_table,
		    workspace, 1, bw ) ;

  /* and that should be that ... */

}


