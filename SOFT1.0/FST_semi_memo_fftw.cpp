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

  FST_semi_memo.c - routines to perform convolutions on the
  2-sphere using a combination of semi-naive and naive algorithms.

  ASSUMES THAT ALL PRECOMPUTED-DATA IS IN MEMORY, AND NOT TO BE
  READ FROM THE DISK.

  The primary functions in this package are

  1) FST_semi_memo() - computes the spherical harmonic expansion.
  2) InvFST_semi_memo() - computes the inverse spherical harmonic transform.
  3) FZT_semi_memo() - computes the zonal harmonic transform.
  4) TransMult() - Multiplies harmonic coefficients using Driscoll-Healy
                    result.  Dual of convolution in "time" domain.
  5) Conv2Sphere_semi_memo() - Convolves two functins defined on the 2-sphere,
                          using seminaive transform

  For descriptions on calling these functions, see the documentation
  preceding each function.

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cospmls.h"
#include "fft_grids.h"
#include "naive_synthesis.h"
#include "primitive.h"
#include "primitive_FST.h"
#include "seminaive.h"
#include "seminaive_fftw.h"
#include <fftw3.h>

#define compmult(a,b,c,d,e,f) (e) = ((a)*(c))-((b)*(d)); (f) = ((a)*(d))+((b)*(c))

/************************************************************************/
/* performs a spherical harmonic transform using the semi-naive
   and naive algorithms */
/* size is the dimension of the input array (size x size) and it is
   expected that size=2*bw.  The inputs rdata and idata are expected
   to be pointers to size x size arrays.  rcoeffs and icoeffs are expected
   to be pointers to bw x bw arrays, and will contain the harmonic
   coefficients in a "linearized" form.


   spharmonic_pml_table should be a (double **) pointer to
   the result of a call to Spharmonic_Pml_Table.  Because this
   table is re-used in the inverse transform, and because for
   timing purposes the computation of the table is not included,
   it is passed in as an argument.  Also, at some point this
   code may be used as par of a series of convolutions, so
   reducing repetitive computation is prioritized.

   spharmonic_pml_table will be an array of (double *) pointers
   the array being of length TableSize(m,bw)

   workspace needs to be a double pointer to an array of size
   (8 * bw^2) + (29 * bw).

   cutoff -> what order to switch from semi-naive to naive
             algorithm.

*/

/* 
   Output Ordering of coeffs f(m,l) is
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
    
   This only requires an array of size (bw*bw).  If zero-padding
   is used to make the indexing nice, then you need a an
   (2bw-1) * bw array - but that is not done here.
   Because of the amount of space necessary for doing
   large transforms, it is important not to use any
   more than necessary.
*/

/*      dataformat =0 -> samples are complex, =1 -> samples real */

void FST_semi_memo_fftw(double *rdata,
						double *rcoeffs, double *icoeffs, 
						int size, double **seminaive_naive_table,
						double *workspace)
{
	int bw, m;
	double *rres, *ires;
	double *rdataptr, *idataptr;
	double *fltres, *scratchpad;
	double *eval_pts;
	double *cos_even;
	int l, dummy ;
	double tmpA, tmpB ;

	bw = size/2;

	/* assign space */
	rres = workspace;                    /* needs (2 * bw * (bw+1)) */
	ires = rres + (size * (size/2+1));   /* needs (2 * bw * (bw+1)) */ 
	fltres = ires + (size * (size/2+1)); /* needs (bw)  */
	eval_pts = fltres + bw;              /* needs (2*bw)  */
	scratchpad = eval_pts + (2*bw);      /* needs (24 * bw)  */
	cos_even = scratchpad + (24*bw);     /* needs (bw) */


	/* total workspace is (4 * bw^2) + (32 * bw) */

	/* do the FFTs along phi */
	fftw_iodim dims1,dims2;
	dims1.n=dims2.n=size;
	dims1.is=dims2.os=1;
	dims1.os=dims2.is=size;
	fftw_plan plan=fftw_plan_guru_split_dft_r2c(1,&dims1,1,&dims2,rdata,rres,ires,FFTW_ESTIMATE);
	fftw_execute(plan);
	fftw_destroy_plan(plan);

	/* point to start of output data buffers */
	rdataptr = rcoeffs;
	idataptr = icoeffs;

	for (m=0; m<bw; m++) {
		/* do the real part */
		SemiNaiveReduced(rres+(m*size), 
			bw, 
			m, 
			fltres, 
			seminaive_naive_table[m],
			scratchpad,
			cos_even);

		/* now load real part of coefficients into output space */  
		memcpy(rdataptr, fltres, sizeof(double) * (bw - m));

		rdataptr += bw-m;

		/* do imaginary part */
		SemiNaiveReduced(ires+(m*size),
			bw, 
			m, 
			fltres, 
			seminaive_naive_table[m],
			scratchpad,
			cos_even);

		/* now load imaginary part of coefficients into output space */  
		memcpy(idataptr, fltres, sizeof(double) * (bw - m));

		idataptr += bw-m;
	}

	/*****************

	New in 2.6: Need to massage the coefficients one more time,
	to get that right factor in there (how embarrassing)

	******************/

	tmpA = 2. * sqrt( PI );
	tmpB = sqrt( 2. * PI );
	tmpA/=size;
	tmpB/=size;

	for(m=0;m<bw;m++)
		for(l=m;l<bw;l++)
		{
			dummy = seanindex(m,l,bw);
			if ( m == 0 )
			{
				rcoeffs[dummy] *= tmpA ;
				icoeffs[dummy] *= tmpA ;
			}
			else
			{
				rcoeffs[dummy] *= tmpB ;
				icoeffs[dummy] *= tmpB ;
			}
		}
}
void FST_semi_memo_fftw(double *data, fftw_complex* coeffs,
						int size, double **seminaive_naive_table,
						double *workspace)
{
	int bw, m;
	double *res;
	double *dataptr;
	double *fltres, *scratchpad;
	double *eval_pts;
	double *cos_even;
	int l, dummy ;
	double tmpA, tmpB ;

	bw = size/2;

	/* assign space */
	res = workspace;                    /* needs (4 * bw * (bw+1)) */
	fltres = res + ( 4 * bw * (bw+1));  /* needs (2 * bw)  */
	eval_pts = fltres + (2*bw);         /* needs (2*bw)  */
	scratchpad = eval_pts + (2*bw);     /* needs (24 * bw)  */
	cos_even = scratchpad + (24*bw);    /* needs (bw) */


	/* total workspace is (4 * bw^2) + (32 * bw) */
	fftw_iodim dims1,dims2;
	fftw_plan plan;
	dims1.n=dims2.n=size;
	dims1.is=dims2.os=1;
	dims1.os=dims2.is=size;
	plan=fftw_plan_guru_dft_r2c(1,&dims1,1,&dims2,data,(fftw_complex*)res,FFTW_ESTIMATE);
	fftw_execute(plan);
	fftw_destroy_plan(plan);

	/* point to start of output data buffers */
	dataptr = (double*)coeffs;
	for (m=0; m<bw; m++) {
		/* do the real part */
		SemiNaiveReduced_fftw(res+(2*m*size), 
			bw, 
			m, 
			fltres, 
			seminaive_naive_table[m],
			scratchpad,
			cos_even);
		/* do imaginary part */
		SemiNaiveReduced_fftw(res+1+(2*m*size),
			bw, 
			m, 
			fltres+1, 
			seminaive_naive_table[m],
			scratchpad,
			cos_even);
		/* now load real and imaginary part of coefficients into output space */  
		memcpy(dataptr, fltres, sizeof(double) * (bw - m) * 2);

		dataptr += (bw-m)*2;
	}

	/*****************

	New in 2.6: Need to massage the coefficients one more time,
	to get that right factor in there (how embarrassing)

	******************/

	tmpA = 2. * sqrt( PI );
	tmpB = sqrt( 2. * PI );
	tmpA/=size;
	tmpB/=size;

	for(m=0;m<bw;m++)
		for(l=m;l<bw;l++)
		{
			dummy = seanindex(m,l,bw);
			if ( m == 0 )
			{
				coeffs[dummy][0] *= tmpA ;
				coeffs[dummy][1] *= tmpA ;
			}
			else
			{
				coeffs[dummy][0] *= tmpB ;
				coeffs[dummy][1] *= tmpB ;
			}
		}
}
void FST_semi_memo_fftw(float *data, fftwf_complex* coeffs,
						int size, float **seminaive_naive_table,
						float *workspace)
{
	int bw, m;
	float *res;
	float *dataptr;
	float *fltres, *scratchpad;
	float *eval_pts;
	float *cos_even;
	int l, dummy ;
	float tmpA, tmpB ;

	bw = size/2;

	/* assign space */
	res = workspace;                    /* needs (4 * bw * (bw+1)) */
	fltres = res + ( 4 * bw * (bw+1));  /* needs (2 * bw)  */
	eval_pts = fltres + (2*bw);         /* needs (2*bw)  */
	scratchpad = eval_pts + (2*bw);     /* needs (24 * bw)  */
	cos_even = scratchpad + (24*bw);    /* needs (bw) */


	/* total workspace is (4 * bw^2) + (32 * bw) */
	fftw_iodim dims1,dims2;
	fftwf_plan plan;
	dims1.n=dims2.n=size;
	dims1.is=dims2.os=1;
	dims1.os=dims2.is=size;
	plan=fftwf_plan_guru_dft_r2c(1,&dims1,1,&dims2,data,(fftwf_complex*)res,FFTW_ESTIMATE);
	fftwf_execute(plan);
	fftwf_destroy_plan(plan);

	/* point to start of output data buffers */
	dataptr = (float*)coeffs;
	for (m=0; m<bw; m++) {
		/* do the real part */
		SemiNaiveReduced_fftw(res+(2*m*size), 
			bw, 
			m, 
			fltres, 
			seminaive_naive_table[m],
			scratchpad,
			cos_even);
		/* do imaginary part */
		SemiNaiveReduced_fftw(res+1+(2*m*size),
			bw, 
			m, 
			fltres+1, 
			seminaive_naive_table[m],
			scratchpad,
			cos_even);
		/* now load real and imaginary part of coefficients into output space */  
		memcpy(dataptr, fltres, sizeof(float) * (bw - m) * 2);

		dataptr += (bw-m)*2;
	}

	/*****************

	New in 2.6: Need to massage the coefficients one more time,
	to get that right factor in there (how embarrassing)

	******************/

	tmpA = float( 2. * sqrt( PI ) );
	tmpB = float( sqrt( 2. * PI ) );
	tmpA/=size;
	tmpB/=size;

	for(m=0;m<bw;m++)
		for(l=m;l<bw;l++)
		{
			dummy = seanindex(m,l,bw);
			if ( m == 0 )
			{
				coeffs[dummy][0] *= tmpA ;
				coeffs[dummy][1] *= tmpA ;
			}
			else
			{
				coeffs[dummy][0] *= tmpB ;
				coeffs[dummy][1] *= tmpB ;
			}
		}
}

void FST_semi_memo_fftw(float *rdata,
		   float *rcoeffs, float *icoeffs, 
		   int size, float **seminaive_naive_table,
		   float *workspace)
{
	int bw, m;
	float *rres, *ires;
	float *rdataptr, *idataptr;
	float *fltres, *scratchpad;
	float *eval_pts;
	float *cos_even;
	int l, dummy ;
	float tmpA, tmpB ;

	bw = size/2;

	/* assign space */
	rres = workspace;                    /* needs (2 * bw * (bw+1)) */
	ires = rres + (size * (size/2+1));   /* needs (2 * bw * (bw+1)) */ 
	fltres = ires + (size * (size/2+1)); /* needs (bw)  */
	eval_pts = fltres + bw;              /* needs (2*bw)  */
	scratchpad = eval_pts + (2*bw);      /* needs (24 * bw)  */
	cos_even = scratchpad + (24*bw);     /* needs (bw) */


	/* total workspace is (4 * bw^2) + (32 * bw) */

	/* do the FFTs along phi */
	fftw_iodim dims1,dims2;
	dims1.n=dims2.n=size;
	dims1.is=dims2.os=1;
	dims1.os=dims2.is=size;
	fftwf_plan plan=fftwf_plan_guru_split_dft_r2c(1,&dims1,1,&dims2,rdata,rres,ires,FFTW_ESTIMATE);
	fftwf_execute(plan);
	fftwf_destroy_plan(plan);

	/* point to start of output data buffers */
	rdataptr = rcoeffs;
	idataptr = icoeffs;

	for (m=0; m<bw; m++) {
		/* do the real part */
		SemiNaiveReduced(rres+(m*size), 
			bw, 
			m, 
			fltres, 
			seminaive_naive_table[m],
			scratchpad,
			cos_even);

		/* now load real part of coefficients into output space */  
		memcpy(rdataptr, fltres, sizeof(float) * (bw - m));

		rdataptr += bw-m;

		/* do imaginary part */
		SemiNaiveReduced(ires+(m*size),
			bw, 
			m, 
			fltres, 
			seminaive_naive_table[m],
			scratchpad,
			cos_even);

		/* now load imaginary part of coefficients into output space */  
		memcpy(idataptr, fltres, sizeof(float) * (bw - m));

		idataptr += bw-m;
	}

	/*****************

	New in 2.6: Need to massage the coefficients one more time,
	to get that right factor in there (how embarrassing)

	******************/

	tmpA = float(2. * sqrt( PI ));
	tmpB = float(sqrt( 2. * PI ));
	tmpA/=size;
	tmpB/=size;

	for(m=0;m<bw;m++)
		for(l=m;l<bw;l++)
		{
			dummy = seanindex(m,l,bw);
			if ( m == 0 )
			{
				rcoeffs[dummy] *= tmpA ;
				icoeffs[dummy] *= tmpA ;
			}
			else
			{
				rcoeffs[dummy] *= tmpB ;
				icoeffs[dummy] *= tmpB ;
			}
		}
}

void InvFST_semi_memo_fftw(fftw_complex *coeffs,double *data, 
					  int size, 
					  double **transpose_seminaive_naive_table,
					  double *workspace)
{
	int bw, m, i, n;
	double *dataptr;
	double *fourdata;
	double *invfltres, *scratchpad;
	double *sin_values, *eval_pts;
	int l, dummy ;
	double tmpA, tmpB ;

	bw = size/2;

	/* allocate space */

	fourdata = workspace;                      /* needs (4 * bw * (bw+1)) */
	invfltres = fourdata + (4 * bw * (bw+1));  /* needs (4 * bw) */
	sin_values = invfltres + (4 * bw);         /* needs (2 * bw) */
	eval_pts = sin_values + (2 * bw);          /* needs (2 * bw) */
	scratchpad = eval_pts + (2 * bw);          /* needs (24 * bw) */

	/* total workspace = (4 * bw^2) + (36 * bw) */

	/* load up the sin_values array */
	n = 2*bw;

	ArcCosEvalPts(n, eval_pts);
	for (i=0; i<n; i++)
		sin_values[i] = sin(eval_pts[i]);

	/**********************

	New in 2.6: Need to massage the coefficients, to get that
	right factor in there (how embarrassing)

	***********************/

	tmpA = 1./(2. * sqrt(PI) );
	tmpB = 1./sqrt(2. * PI ) ;

	for(m=0;m<bw;m++)
		for(l=m;l<bw;l++)
		{
			dummy = seanindex(m,l,bw);
			if ( m == 0 )
			{
				coeffs[dummy][0] *= tmpA ;
				coeffs[dummy][1] *= tmpA ;
			}
			else
			{
				coeffs[dummy][0] *= tmpB ;
				coeffs[dummy][1] *= tmpB ;
			}
		}

		/* Now do all of the inverse Legendre transforms */
		dataptr = (double*)coeffs;

		for (m=0; m<bw; m++)
		{
			/* do real part first */ 
			InvSemiNaiveReduced_fftw(dataptr,
				bw,
				m,
				invfltres,
				transpose_seminaive_naive_table[m],
				sin_values,
				scratchpad);

			/* now do imaginary part */

			InvSemiNaiveReduced_fftw(dataptr+1,
				bw,
				m,
				invfltres+2*bw,
				transpose_seminaive_naive_table[m],
				sin_values,
				scratchpad);

			/* will store normal, then tranpose before doing inverse fft */
			double *temp=fourdata+(2*m*size);
			for(int i=0;i<size;i++){
				temp[2*i  ]=invfltres[i];
				temp[2*i+1]=invfltres[size+i];
			}

			/* move to next set of coeffs */

			dataptr += (bw-m)*2;

		} /* closes m loop */


		/* now fill in zero values where m = bw (from problem definition) */
		memset(fourdata + (2*bw*size),0,sizeof(double) * size *2);

		/* now do inverse fourier grid computation */
		fftw_iodim dims1,dims2;
		dims1.n=dims2.n=size;
		dims1.is=dims2.os=size;
		dims1.os=dims2.is=1;
		// Generate the FFTW plan
		fftw_plan plan=fftw_plan_guru_dft_c2r(1,&dims1,1,&dims2, (fftw_complex*)fourdata, data,FFTW_ESTIMATE);
		// do the FFTs along phi
		fftw_execute(plan);
		fftw_destroy_plan(plan);
		/* amscray */
}

void InvFST_semi_memo_fftw(fftwf_complex *coeffs,float *data, 
					  int size, 
					  float **transpose_seminaive_naive_table,
					  float *workspace)
{
	int bw, m, i, n;
	float *dataptr;
	float *fourdata;
	float *invfltres, *scratchpad;
	float *sin_values, *eval_pts;
	int l, dummy ;
	float tmpA, tmpB ;

	bw = size/2;

	/* allocate space */

	fourdata = workspace;                      /* needs (4 * bw * (bw+1)) */
	invfltres = fourdata + (4 * bw * (bw+1));  /* needs (4 * bw) */
	sin_values = invfltres + (4 * bw);         /* needs (2 * bw) */
	eval_pts = sin_values + (2 * bw);          /* needs (2 * bw) */
	scratchpad = eval_pts + (2 * bw);          /* needs (24 * bw) */

	/* total workspace = (4 * bw^2) + (36 * bw) */

	/* load up the sin_values array */
	n = 2*bw;

	ArcCosEvalPts(n, eval_pts);
	for (i=0; i<n; i++)
		sin_values[i] = sin(eval_pts[i]);

	/**********************

	New in 2.6: Need to massage the coefficients, to get that
	right factor in there (how embarrassing)

	***********************/

	tmpA = float( 1./(2. * sqrt(PI) ) );
	tmpB = float( 1./sqrt(2. * PI )   );

	for(m=0;m<bw;m++)
		for(l=m;l<bw;l++)
		{
			dummy = seanindex(m,l,bw);
			if ( m == 0 )
			{
				coeffs[dummy][0] *= tmpA ;
				coeffs[dummy][1] *= tmpA ;
			}
			else
			{
				coeffs[dummy][0] *= tmpB ;
				coeffs[dummy][1] *= tmpB ;
			}
		}

		/* Now do all of the inverse Legendre transforms */
		dataptr = (float*)coeffs;

		for (m=0; m<bw; m++)
		{
			/* do real part first */ 
			InvSemiNaiveReduced_fftw(dataptr,
				bw,
				m,
				invfltres,
				transpose_seminaive_naive_table[m],
				sin_values,
				scratchpad);

			/* now do imaginary part */

			InvSemiNaiveReduced_fftw(dataptr+1,
				bw,
				m,
				invfltres+2*bw,
				transpose_seminaive_naive_table[m],
				sin_values,
				scratchpad);

			/* will store normal, then tranpose before doing inverse fft */
			float *temp=fourdata+(2*m*size);
			for(int i=0;i<size;i++){
				temp[2*i  ]=invfltres[i];
				temp[2*i+1]=invfltres[size+i];
			}

			/* move to next set of coeffs */

			dataptr += (bw-m)*2;

		} /* closes m loop */


		/* now fill in zero values where m = bw (from problem definition) */
		memset(fourdata + (2*bw*size),0,sizeof(float) * size *2);

		/* now do inverse fourier grid computation */
		fftw_iodim dims1,dims2;
		dims1.n=dims2.n=size;
		dims1.is=dims2.os=size;
		dims1.os=dims2.is=1;
		// Generate the FFTW plan
		fftwf_plan plan=fftwf_plan_guru_dft_c2r(1,&dims1,1,&dims2, (fftwf_complex*)fourdata, data,FFTW_ESTIMATE);
		// do the FFTs along phi
		fftwf_execute(plan);
		fftwf_destroy_plan(plan);
		/* amscray */

}

void InvFST_semi_memo_fftw(double *rcoeffs, double *icoeffs, 
					  double *rdata, 
					  int size, 
					  double **transpose_seminaive_naive_table,
					  double *workspace)
{
	int bw, m, i, n;
	double *rdataptr, *idataptr;
	double *rfourdata, *ifourdata;
	double *rinvfltres, *iminvfltres, *scratchpad;
	double *sin_values, *eval_pts;
	int l, dummy ;
	double tmpA, tmpB ;

	bw = size/2;

	/* allocate space */

	rfourdata = workspace;                      /* needs (size * (bw+1)) */
	ifourdata = rfourdata + (2 * bw * (bw+1));  /* needs (size * (bw+1)) */
	rinvfltres = ifourdata + (2 * bw * (bw+1)); /* needs (2 * bw) */
	iminvfltres = rinvfltres + (2 * bw);        /* needs (2 * bw) */
	sin_values = iminvfltres + (2 * bw);        /* needs (2 * bw) */
	eval_pts = sin_values + (2 * bw);           /* needs (2 * bw) */
	scratchpad = eval_pts + (2 * bw);           /* needs (24 * bw) */

	/* total workspace = (4 * bw^2) + (36 * bw) */

	/* load up the sin_values array */
	n = 2*bw;

	ArcCosEvalPts(n, eval_pts);
	for (i=0; i<n; i++)
		sin_values[i] = sin(eval_pts[i]);

	/**********************

	New in 2.6: Need to massage the coefficients, to get that
	right factor in there (how embarrassing)

	***********************/

	tmpA = 1./(2. * sqrt(PI) );
	tmpB = 1./sqrt(2. * PI ) ;

	for(m=0;m<bw;m++)
		for(l=m;l<bw;l++)
		{
			dummy = seanindex(m,l,bw);
			if ( m == 0 )
			{
				rcoeffs[dummy] *= tmpA ;
				icoeffs[dummy] *= tmpA ;
			}
			else
			{
				rcoeffs[dummy] *= tmpB ;
				icoeffs[dummy] *= tmpB ;
			}
		}

		/* Now do all of the inverse Legendre transforms */
		rdataptr = rcoeffs;
		idataptr = icoeffs;
		for (m=0; m<bw; m++)
		{
			/* do real part first */ 
			InvSemiNaiveReduced(rdataptr,
				bw,
				m,
				rinvfltres,
				transpose_seminaive_naive_table[m],
				sin_values,
				scratchpad);

			/* now do imaginary part */

			InvSemiNaiveReduced(idataptr,
				bw,
				m,
				iminvfltres,
				transpose_seminaive_naive_table[m],
				sin_values,
				scratchpad);

			/* will store normal, then tranpose before doing inverse fft */
			memcpy(rfourdata+(m*size), rinvfltres, sizeof(double) * size);
			memcpy(ifourdata+(m*size), iminvfltres, sizeof(double) * size);

			/* move to next set of coeffs */

			rdataptr += bw-m;
			idataptr += bw-m;

		} /* closes m loop */


		/* now fill in zero values where m = bw (from problem definition) */
		memset(rfourdata + (bw * size), 0, sizeof(double) * size);
		memset(ifourdata + (bw * size), 0, sizeof(double) * size);

		/* now do inverse fourier grid computation */
		fftw_iodim dims1,dims2;
		dims1.n=dims2.n=size;
		dims1.is=dims2.os=size;
		dims1.os=dims2.is=1;
		// Generate the FFTW plan
		fftw_plan plan=fftw_plan_guru_split_dft_c2r(1,&dims1,1,&dims2, rfourdata, ifourdata, rdata,FFTW_ESTIMATE);
		// do the FFTs along phi
		fftw_execute(plan);
		fftw_destroy_plan(plan);
		/* amscray */

}
void InvFST_semi_memo_fftw(float *rcoeffs, float *icoeffs, 
					  float *rdata, 
					  int size, 
					  float **transpose_seminaive_naive_table,
					  float *workspace)
{
	int bw, m, i, n;
	float *rdataptr, *idataptr;
	float *rfourdata, *ifourdata;
	float *rinvfltres, *iminvfltres, *scratchpad;
	float *sin_values, *eval_pts;
	int l, dummy ;
	float tmpA, tmpB ;

	bw = size/2;

	/* allocate space */

	rfourdata = workspace;                      /* needs (size * (bw+1)) */
	ifourdata = rfourdata + (2 * bw * (bw+1));  /* needs (size * (bw+1)) */
	rinvfltres = ifourdata + (2 * bw * (bw+1)); /* needs (2 * bw) */
	iminvfltres = rinvfltres + (2 * bw);        /* needs (2 * bw) */
	sin_values = iminvfltres + (2 * bw);        /* needs (2 * bw) */
	eval_pts = sin_values + (2 * bw);           /* needs (2 * bw) */
	scratchpad = eval_pts + (2 * bw);           /* needs (24 * bw) */

	/* total workspace = (4 * bw^2) + (36 * bw) */

	/* load up the sin_values array */
	n = 2*bw;

	ArcCosEvalPts(n, eval_pts);
	for (i=0; i<n; i++)
		sin_values[i] = sin(eval_pts[i]);

	/**********************

	New in 2.6: Need to massage the coefficients, to get that
	right factor in there (how embarrassing)

	***********************/

	tmpA = float(1./(2. * sqrt(PI)));
	tmpB = float(1./sqrt(2. * PI )) ;

	for(m=0;m<bw;m++)
		for(l=m;l<bw;l++)
		{
			dummy = seanindex(m,l,bw);
			if ( m == 0 )
			{
				rcoeffs[dummy] *= tmpA ;
				icoeffs[dummy] *= tmpA ;
			}
			else
			{
				rcoeffs[dummy] *= tmpB ;
				icoeffs[dummy] *= tmpB ;
			}
		}

		/* Now do all of the inverse Legendre transforms */
		rdataptr = rcoeffs;
		idataptr = icoeffs;
		for (m=0; m<bw; m++)
		{
			/* do real part first */ 
			InvSemiNaiveReduced(rdataptr,
				bw,
				m,
				rinvfltres,
				transpose_seminaive_naive_table[m],
				sin_values,
				scratchpad);

			/* now do imaginary part */

			InvSemiNaiveReduced(idataptr,
				bw,
				m,
				iminvfltres,
				transpose_seminaive_naive_table[m],
				sin_values,
				scratchpad);

			/* will store normal, then tranpose before doing inverse fft */
			memcpy(rfourdata+(m*size), rinvfltres, sizeof(float) * size);
			memcpy(ifourdata+(m*size), iminvfltres, sizeof(float) * size);

			/* move to next set of coeffs */

			rdataptr += bw-m;
			idataptr += bw-m;

		} /* closes m loop */


		/* now fill in zero values where m = bw (from problem definition) */
		memset(rfourdata + (bw * size), 0, sizeof(float) * size);
		memset(ifourdata + (bw * size), 0, sizeof(float) * size);

		/* now do inverse fourier grid computation */
		fftw_iodim dims1,dims2;
		dims1.n=dims2.n=size;
		dims1.is=dims2.os=size;
		dims1.os=dims2.is=1;
		// Generate the FFTW plan
		fftwf_plan plan=fftwf_plan_guru_split_dft_c2r(1,&dims1,1,&dims2, rfourdata, ifourdata, rdata,FFTW_ESTIMATE);
		// do the FFTs along phi
		fftwf_execute(plan);
		fftwf_destroy_plan(plan);
		/* amscray */

}
