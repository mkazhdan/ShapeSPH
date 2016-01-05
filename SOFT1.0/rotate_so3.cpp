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




*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "complex.h"
#include "primitive_FST.h"
#include "rotate_so3.h"
#include "FST_semi_memo.h"

/*****************************************************************/
/*****************************************************************/
/*
  genExp: given an angle phi, and an bandlimit bw, this routine
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

void genExp( int bw ,
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

 wignerdmat: given a wigner little-d matrix of degree L-1 evaluated
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

void wignerdmat( int L,
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
  rotateCoefDegree: given the degree L spherical coefficients {f_L^m},
                    where -L <= m <= L of some function on the sphere,
		    and three Euler angles alpha, beta, gamma, multiply
		    those coefficients by the wigner D-mat(alpha,beta,gamma).
		    In other words, I'm rotating f by the Euler angles,
		    and I'm massaging the coefficients as representation
		    theory says I should.
		    
  Arguments:

        INPUT:
	
	L: degree of coefficients (i.e. the representation)

	cfR, cfI: two 1-D REAL arrays, each of length 2L+1,
	      containing the real and imaginary parts of
	      the coefficients of the signal to be rotated.
	      Each is in the following order:

	      [ f_L^(-L) f_L^(-L+1) ... f_L^(-1) 0 f_L^(1) ... f_L^(L) ]

        alpha, gamma -> the Euler angles


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

void rotateCoefDegree( int L,
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
  rotateCoefAll: given the spherical harmonic coefficients of
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
		 

  arguments:

           bwIn: bandlimit of the function being rotated, i.e.
	         what's the bandlimit of the function you want
		 to rotate?

	   degOut: Now, instead of rotating the coefficients of
	           all the degrees 0 <= L < bwIn, for whatever
		   reason you may want to go only "so far" and
		   rotate only through L = degOut. The rotated
		   coefficients for degrees degOut+1 <= L < bwIn
		   will be set to 0.

	   bwOut: Now, for whatever reason, you may want to "dumb down"
	          the rotated result. That is, you might want the
		  rotated function to be of bandlimit bwOut <= bwIn.
		  Now you can.

	   alpha, beta, gamma: the Euler angles

	   cfR, cfI: two REAL arrays containing the real and imaginary
	          parts of the spherical coefficients of f, each of
		  size (bwIn)^2

	   cfrR, cfrI: two REAL arrays containing the real and imaginary
	          parts of the spherical coefficients of *rotated* f,
		  each of size (bwIn)^2

	   scratch: REAL scratch space of size
	            2*bwIn + 2*(2*bwIn-1)^2 + (2*degOut)^2 +
		    2*4*(2*bwIn-1)


  NOTE: Like bwIn, bwOut must also be a power of 2 !

  NOTE:	The coefficients in cf will be in the same order
        as that produced by a forward spherical transform
	routine found in SpharmonicKit. This will require
	the use of the indexing functions
	
	seanindex
	seanindex2
	
	defined in primitive_FST.[c,h].


  NOTE: This function is an especially egregious memory hog. If
        you can guarantee that bwIn will always equal bwOut, and if
	you don't mind writing over the spherical coefficients of
	the input signal, you can lessen the memory requirements
	alot. By how much? I'll have to figure that out later, but
	I know it is alot.

  NOTE: Since the routines in SpharmonicKit produce two separate double
        arrays, for the real and imaginary parts of the coefficients,
	they will need to be interleaved prior to entering rotateCoefAll.
	Sorry.

*/

void rotateCoefAll( int bwIn,
		    int bwOut,
		    int degOut,
		    REAL alpha,
		    REAL beta,
		    REAL gamma,
		    REAL *cfR, REAL *cfI,
		    REAL *cfrR, REAL *cfrI,
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
  matIn = sqrts + 2*bwIn ;
  matOut = matIn + (2*bwIn-1)*(2*bwIn-1);

  cfTmpR = matOut + (2*bwIn-1)*(2*bwIn-1);
  cfTmpI = cfTmpR + (2*bwIn - 1);
  cfrTmpR = cfTmpI + (2*bwIn - 1);
  cfrTmpI = cfrTmpR + (2*bwIn - 1);
  expAR = cfrTmpI + (2*bwIn - 1);
  expAI = expAR + (2*bwIn - 1);
  expGR = expAI + (2*bwIn - 1);
  expGI = expGR + (2*bwIn - 1);
  workspace = expGI + (2*bwIn - 1) ;
  /* workspace still needs (2*degOut)^2 space */

  /*
    We'll first precompute stuff needed for the
    wignerdmat routine: the square roots, and the trig values
  */

  for ( i = 0 ; i < 2*bwIn ; i ++ )
    sqrts[i] = (REAL) sqrt( (REAL) i ) ;

  trigs[0] = (REAL) cos(0.5*beta);
  trigs[1] = (REAL) sin(0.5*beta);
 
  /*
    For the rotateCoefDegree routine, we'll need the values
    of the complex exponentials
  */

  genExp( bwIn, alpha, expAR, expAI ) ;
  genExp( bwIn, gamma, expGR, expGI ) ;


  /*
    Now, to massage the coefficients, one degree at a time

    Note that this will be done for degrees L = 0 through
    L = degOut (provided that also L < bwOut). The remaining
    degrees will be set to 0.
  */

  /*  for ( deg = 0 ; (deg <= degOut) && (deg < bwOut) ; deg ++ ) */

  for ( deg = 0 ; (deg <= degOut) && (deg < bwIn) ; deg ++ ) 
    {
      dim = 2*deg - 1 ;
      
      if ( deg > 0 )
	memcpy( matIn, matOut, sizeof(REAL)*(dim*dim) );

      /* grab the spherical coefficients at this degree */
      for ( m = -deg ; m < deg+1 ; m++ )
	{
	  dummy = seanindex(m, deg, bwIn);
	  cfTmpR[m+deg] = cfR[dummy];
	  cfTmpI[m+deg] = cfI[dummy];
	}

      /* generate the required wigner-d mat */
      wignerdmat( deg, matIn, matOut, trigs, sqrts, workspace );

      /* massage the coefficients at this degree */
      rotateCoefDegree( deg,
			cfTmpR, cfTmpI,
			cfrTmpR, cfrTmpI,
			expAR, expAI,
			expGR, expGI,
			bwIn - 1,
			matOut ) ;

      /* save the rotated coefficients */
      for ( m = -deg ; m < deg+1 ; m++ )
	{
	  dummy = seanindex(m, deg, bwOut);
	  cfrR[dummy] = cfrTmpR[m+deg];
	  cfrI[dummy] = cfrTmpI[m+deg];
	}
    }

  /* zero out remaining coefficients */
  /*  for ( deg = degOut + 1 ; deg < bwOut ; deg ++ ) */
  for ( ; deg < bwOut ; deg ++ )
    for ( m = -deg ; m < deg + 1 ; m ++ )
      {
	dummy = seanindex(m, deg, bwOut);
	cfrR[dummy] = (REAL) 0;
	cfrI[dummy] = (REAL) 0;
      }
}

/*****************************************************************/
/*****************************************************************/
/*
  rotateFct: given a function defined on the sphere,
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

           bwIn: bandlimit of the inputed function

	   bwOut: bandlimit of the outputed, rotated function

	   degOut: Now, instead of rotating the coefficients of
	           all the degrees 0 <= L < bwIn, for whatever
		   reason you may want to go only "so far" and
		   rotate only through L = degOut. The rotated
		   coefficients for degrees degOut+1 <= L < bwIn
		   will be considered to be 0.

           sigInR: REAL array of the *real* parts of the signal sampled
	           on the usual 2bwIn * 2bwIn grid on the sphere

	   sigInI: REAL array of the *imaginary* parts of the signal
	           sampled on the usual 2bwIn * 2bwIn grid on the sphere

           sigOutR: REAL array of the *real* parts of the rotated signal
	            sampled on the usual 2bwOut * 2bwOut grid on the sphere

	   sigOutI: REAL array of the *imaginary* parts of the rotated
	            signal sampled on the usual 2bwOut * 2bwOut grid on
		    the sphere

	   alpha, beta, gamma: the Euler angles defining the rotation

	   scratch: REAL array of size (gulp)
		    14*bwIn^2 + 48*bwIn - 8

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

void rotateFct( int bwIn, int bwOut, int degOut,
		REAL *sigInR, REAL *sigInI,
		REAL *sigOutR, REAL *sigOutI,
		REAL alpha, REAL beta, REAL gamma,
		REAL *scratch,
		REAL **spharmonic_pml_table,
		REAL **transpose_spharmonic_pml_table )
{
  REAL *rcoeffs, *icoeffs, *workspace ;
  REAL *coefOutR, *coefOutI ;
  int i ;

  rcoeffs = scratch ;
  icoeffs = rcoeffs + (bwIn*bwIn) ;
  coefOutR = icoeffs + (bwIn*bwIn) ;
  coefOutI = coefOutR + (bwOut*bwOut) ;
  workspace = coefOutI + (bwOut*bwOut) ;


  for(i = 0 ; i < (bwOut*bwOut) ; i ++ )
    {
      coefOutR[i] = (REAL) 0 ;
      coefOutI[i] = (REAL) 0 ;
    }

  /* compute spherical coefficients of input signal */
  FST_semi_memo( sigInR, sigInI,
		 rcoeffs, icoeffs,
		 2*bwIn, spharmonic_pml_table,
		 workspace, 1, bwIn ) ;

  /* now massage the coefficients */
  /* note that I'm using the workspace array again */
  rotateCoefAll( bwIn, bwOut, degOut,
		 alpha, beta, gamma,
		 rcoeffs, icoeffs,
		 coefOutR, coefOutI,
		 workspace ) ;


  /* take inverse spherical transform */
  InvFST_semi_memo( coefOutR, coefOutI,
		    sigOutR, sigOutI,
		    2*bwOut,
		    transpose_spharmonic_pml_table,
		    workspace, 1, bwOut ) ;

  /* and that should be that ... */

}


