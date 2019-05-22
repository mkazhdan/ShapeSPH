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


/* Source code for computing the Legendre transform where
   projections are carried out in cosine space, i.e., the
   "seminaive" algorithm.

   For a description, see the related paper or Sean's thesis.

*/


/************************************************************************

  NOTE THAT the inverse discrete Legendre transform routine

                      InvSemiNaiveReduced

  does *not* use *any* FFTPACK routines. This is because I need
  to take N coefficients and return 2*N samples. This requires
  an inverse dct routine that can return an output array whose
  length is greater than the input array. I don't know how to
  do that with FFTPACK's inverse dct routine.
  
  If I have to do an inverse dct whose input and output are the
  same lengths, *then* I can use FFTPACK's routine.

  ********************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>   /** for memcpy **/

#include "cospmls.h"
#include "csecond.h"

/* even if I use fftpack, I still need newFCT for
   its inverse dct (which allows for output longer
   than input); it's used in the inverse seminaive
   algorithm */

#include "newFCT.h"
#include "oddweights.h"
#include "weights.h"


/********************************************************************

  First define those routines which *do not* use any functions
  defined in FFTPACK.

  *******************************************************************/

/***************************************************
  Just like SemiNaiveReduced EXCEPT I assume that
  the input data has ALREADY BEEN WEIGHTED AND
  COSINE-TRANSFORMED. So I do not have to use any
  functions defined in FFTPACK. This routine is
  expected to be used in the hybrid fast legendre
  transform.
  
  Also, note that there's no "workspace" argument:
  I don't need it in the function.

  k = how many coefficient do you want to compute ?

  **********************/

/************************************************************************/
/* InvSemiNaiveReduced computes the inverse Legendre transform
   using the transposed seminaive algorithm.  Note that because
   the Legendre transform is orthogonal, the inverse can be
   computed by transposing the matrix formulation of the
   problem.

   The forward transform looks like

   l = PCWf

   where f is the data vector, W is a quadrature matrix,
   C is a cosine transform matrix, P is a matrix
   full of coefficients of the cosine series representation
   of each Pml function P(m,m) P(m,m+1) ... P(m,bw-1),
   and l is the (associated) Legendre series representation
   of f.

   So to do the inverse, you do

   f = trans(C) trans(P) l

   so you need to transpose the matrix P from the forward transform
   and then do a cosine series evaluation.  No quadrature matrix
   is necessary.  If order m is odd, then there is also a sin
   factor that needs to be accounted for.

   Note that this function was written to be part of a full
   spherical harmonic transform, so a lot of precomputation
   has been assumed.

   Input argument description

   assoc_legendre_series - a double pointer to an array of length
                           (bw-m) containing associated
			   Legendre series coefficients.  Assumed
			   that first entry contains the P(m,m)
			   coefficient.

   bw - problem bandwidth, assumed to be a power of two.
   m - order of the associated Legendre functions
   result - a double pointer to an array of (2*bw) samples
            representing the evaluation of the Legendre
	    series at (2*bw) Chebyshev nodes.
   trans_cos_pml_table - double pointer to array representing
                         the linearized form of trans(P) above.
			 See cospmls.{h,c} for a description
			 of the function Transpose_CosPmlTableGen()
			 which generates this array.
   sin_values - when m is odd, need to factor in the sin(x) that
                is factored out of the generation of the values
		in trans(P).
   workspace - a double array of size 5*bw
   

*/
void SemiNaiveReduced_fftw(double *data, 
						   int bw, 
						   int m, 
						   double *result, 
						   double *cos_pml_table, 
						   double *workspace,
						   double *cos_even)
{
	int i, j, n;
	const double *weights;
	double *weighted_data;
	double *cos_data;
	double *scratchpad;

	double *cos_odd;

	double eresult0, eresult1, eresult2, eresult3;
	double oresult0, oresult1;  
	double *pml_ptr_even, *pml_ptr_odd;

	int ctr;
	double d_bw;

	n = 2*bw;
	d_bw = (double) bw;

	cos_odd = cos_even + (bw/2);

	/* assign workspace */
	weighted_data = workspace;
	cos_data = workspace + (2*bw);
	scratchpad = cos_data + (2*bw); /* scratchpad needs (4*bw) */
	/* total workspace = (8 * bw) */

	/*
	need to apply quadrature weights to the data and compute
	the cosine transform: if the order is even, get the regular
	weights; if the order is odd, get the oddweights (those
	where I've already multiplied the sin factor in
	*/

	if( (m % 2) == 0)
		weights = get_weights(bw);
	else
		weights = get_oddweights(bw);


	/*
	smooth the weighted signal
	*/
	for ( i = 0; i < n    ; ++i )
		weighted_data[i] = data[ i*2 ] * weights[ i ];

	kFCT(weighted_data, cos_data, scratchpad, n, bw, 1);

	/* normalize the data for this problem */
	if (m == 0)
	{
		for (i=0; i<bw; i++)
			cos_data[i] *= (d_bw * 0.5);
	}
	else
	{
		for (i=0; i<bw; i++)
			cos_data[i] *= d_bw ;
	}

	cos_data[0] *= 2.0;


	/*** separate cos_data into even and odd indexed arrays ***/
	ctr = 0;
	for(i = 0; i < bw; i += 2)
	{
		cos_even[ctr] = cos_data[i];
		cos_odd[ctr]  = cos_data[i+1];
		ctr++;
	}

	/*
	do the projections; Note that the cos_pml_table has
	had all the zeroes stripped out so the indexing is
	complicated somewhat
	*/
	/*
	this loop is unrolled to compute two coefficients
	at a time (for efficiency)
	*/

	for(i = 0; i < (bw - m) - ( (bw - m) % 2 ); i += 2)
	{
		/* get ptrs to precomputed data */
		pml_ptr_even = cos_pml_table + NewTableOffset(m, m + i);
		pml_ptr_odd  = cos_pml_table + NewTableOffset(m, m + i + 1);

		eresult0 = 0.0; eresult1 = 0.0;
		oresult0 = 0.0; oresult1 = 0.0;

		for(j = 0; j < (((m + i)/2) + 1) % 2; ++j)
		{
			eresult0 += cos_even[j] * pml_ptr_even[j];
			oresult0 += cos_odd[j] * pml_ptr_odd[j];
		}
		for( ;  j < ((m + i)/2) + 1; j += 2)
		{
			eresult0 += cos_even[j] * pml_ptr_even[j];
			eresult1 += cos_even[j + 1] * pml_ptr_even[j + 1];
			oresult0 += cos_odd[j] * pml_ptr_odd[j];
			oresult1 += cos_odd[j + 1] * pml_ptr_odd[j + 1];
		}

		/* save the result */
		result[2*i] = eresult0 + eresult1;
		result[2*(i + 1)] = oresult0 + oresult1;
	}

	/* if there are an odd number of coefficients to compute
	at this order, get that last coefficient! */
	if( ( (bw - m) % 2 ) == 1)
	{
		pml_ptr_even = cos_pml_table + NewTableOffset(m, bw - 1);

		eresult0 = 0.0; eresult1 = 0.0;
		eresult2 = 0.0; eresult3 = 0.0;

		for(j = 0; j < (((bw - 1)/2) + 1) % 4; ++j)
		{
			eresult0 += cos_even[j] * pml_ptr_even[j];
		}

		for( ;  j < ((bw - 1)/2) + 1; j += 4)
		{
			eresult0 += cos_even[j] * pml_ptr_even[j];
			eresult1 += cos_even[j + 1] * pml_ptr_even[j + 1];
			eresult2 += cos_even[j + 2] * pml_ptr_even[j + 2];
			eresult3 += cos_even[j + 3] * pml_ptr_even[j + 3];
		}
		result[(bw - m - 1)*2] = eresult0 + eresult1 + eresult2 + eresult3;
	}
}



void SemiNaiveReduced_fftw(float *data, 
						   int bw, 
						   int m, 
						   float *result, 
						   float *cos_pml_table, 
						   float *workspace,
						   float *cos_even)
{
	int i, j, n;
	const double *weights;
	float *weighted_data;
	float *cos_data;
	float *scratchpad;

	float *cos_odd;

	float eresult0, eresult1, eresult2, eresult3;
	float oresult0, oresult1;  
	float *pml_ptr_even, *pml_ptr_odd;

	int ctr;
	float d_bw;

	n = 2*bw;
	d_bw = (float) bw;

	cos_odd = cos_even + (bw/2);

	/* assign workspace */
	weighted_data = workspace;
	cos_data = workspace + (2*bw);
	scratchpad = cos_data + (2*bw); /* scratchpad needs (4*bw) */
	/* total workspace = (8 * bw) */

	/*
	need to apply quadrature weights to the data and compute
	the cosine transform: if the order is even, get the regular
	weights; if the order is odd, get the oddweights (those
	where I've already multiplied the sin factor in
	*/

	if( (m % 2) == 0)
		weights = get_weights(bw);
	else
		weights = get_oddweights(bw);


	/*
	smooth the weighted signal
	*/
	for ( i = 0; i < n    ; ++i )
		weighted_data[i] = float(data[ i*2 ] * weights[ i ]);

	kFCT(weighted_data, cos_data, scratchpad, n, bw, 1);

	/* normalize the data for this problem */
	if (m == 0)
	{
		for (i=0; i<bw; i++)
			cos_data[i] *= float(d_bw * 0.5);
	}
	else
	{
		for (i=0; i<bw; i++)
			cos_data[i] *= d_bw ;
	}

	cos_data[0] *= 2.0;


	/*** separate cos_data into even and odd indexed arrays ***/
	ctr = 0;
	for(i = 0; i < bw; i += 2)
	{
		cos_even[ctr] = cos_data[i];
		cos_odd[ctr]  = cos_data[i+1];
		ctr++;
	}

	/*
	do the projections; Note that the cos_pml_table has
	had all the zeroes stripped out so the indexing is
	complicated somewhat
	*/
	/*
	this loop is unrolled to compute two coefficients
	at a time (for efficiency)
	*/

	for(i = 0; i < (bw - m) - ( (bw - m) % 2 ); i += 2)
	{
		/* get ptrs to precomputed data */
		pml_ptr_even = cos_pml_table + NewTableOffset(m, m + i);
		pml_ptr_odd  = cos_pml_table + NewTableOffset(m, m + i + 1);

		eresult0 = 0.0; eresult1 = 0.0;
		oresult0 = 0.0; oresult1 = 0.0;

		for(j = 0; j < (((m + i)/2) + 1) % 2; ++j)
		{
			eresult0 += cos_even[j] * pml_ptr_even[j];
			oresult0 += cos_odd[j] * pml_ptr_odd[j];
		}
		for( ;  j < ((m + i)/2) + 1; j += 2)
		{
			eresult0 += cos_even[j] * pml_ptr_even[j];
			eresult1 += cos_even[j + 1] * pml_ptr_even[j + 1];
			oresult0 += cos_odd[j] * pml_ptr_odd[j];
			oresult1 += cos_odd[j + 1] * pml_ptr_odd[j + 1];
		}

		/* save the result */
		result[2*i] = eresult0 + eresult1;
		result[2*(i + 1)] = oresult0 + oresult1;
	}

	/* if there are an odd number of coefficients to compute
	at this order, get that last coefficient! */
	if( ( (bw - m) % 2 ) == 1)
	{
		pml_ptr_even = cos_pml_table + NewTableOffset(m, bw - 1);

		eresult0 = 0.0; eresult1 = 0.0;
		eresult2 = 0.0; eresult3 = 0.0;

		for(j = 0; j < (((bw - 1)/2) + 1) % 4; ++j)
		{
			eresult0 += cos_even[j] * pml_ptr_even[j];
		}

		for( ;  j < ((bw - 1)/2) + 1; j += 4)
		{
			eresult0 += cos_even[j] * pml_ptr_even[j];
			eresult1 += cos_even[j + 1] * pml_ptr_even[j + 1];
			eresult2 += cos_even[j + 2] * pml_ptr_even[j + 2];
			eresult3 += cos_even[j + 3] * pml_ptr_even[j + 3];
		}
		result[(bw - m - 1)*2] = eresult0 + eresult1 + eresult2 + eresult3;
	}
}
void InvSemiNaiveReduced_fftw(double *assoc_legendre_series, 
							  int bw, 
							  int m, 
							  double *result, 
							  double *trans_cos_pml_table, 
							  double *sin_values,
							  double *workspace)
{
	double *fcos; /* buffer for cosine series rep of result */
	double *trans_tableptr;
	double *scratchpad;
	double *assoc_offset;
	int i, j, n, rowsize;

	double *p;
	double fcos0, fcos1, fcos2, fcos3;

	fcos = workspace; /* needs bw */
	scratchpad = fcos + bw; /* needs (4*bw) */

	/* total workspace is 5*bw */

	trans_tableptr = trans_cos_pml_table;
	p = trans_cos_pml_table;

	/* main loop - compute each value of fcos

	Note that all zeroes have been stripped out of the
	trans_cos_pml_table, so indexing is somewhat complicated.
	*/

	for (i=0; i<=m; i++)
	{
		rowsize = Transpose_RowSize(i, m, bw);

		/* need to point to correct starting value of Legendre series */
		assoc_offset = assoc_legendre_series + (i % 2)*2;

		fcos0 = 0.0 ; fcos1 = 0.0; fcos2 = 0.0; fcos3 = 0.0;

		for (j = 0; j < rowsize % 4; ++j)
			fcos0 += assoc_offset[2*2*j] * trans_tableptr[j];

		for ( ; j < rowsize; j += 4){
			fcos0 += assoc_offset[2*2*j] * trans_tableptr[j];
			fcos1 += assoc_offset[2*2*(j+1)] * trans_tableptr[j+1];
			fcos2 += assoc_offset[2*2*(j+2)] * trans_tableptr[j+2];
			fcos3 += assoc_offset[2*2*(j+3)] * trans_tableptr[j+3];
		}

		fcos[i] = fcos0 + fcos1 + fcos2 + fcos3;

		trans_tableptr += rowsize;      
	}

	for (i=m+1; i<bw-1; i++)
	{

		rowsize = Transpose_RowSize(i, m, bw);

		/* need to point to correct starting value of Legendre series */
		assoc_offset = assoc_legendre_series + ((i - m) + (m % 2))*2;

		fcos0 = 0.0 ; fcos1 = 0.0; fcos2 = 0.0; fcos3 = 0.0;

		for (j = 0; j < rowsize % 4; ++j)
			fcos0 += assoc_offset[2*2*j] * trans_tableptr[j];

		for ( ; j < rowsize; j += 4){
			fcos0 += assoc_offset[2*2*j] * trans_tableptr[j];
			fcos1 += assoc_offset[2*2*(j+1)] * trans_tableptr[j+1];
			fcos2 += assoc_offset[2*2*(j+2)] * trans_tableptr[j+2];
			fcos3 += assoc_offset[2*2*(j+3)] * trans_tableptr[j+3];
		}
		fcos[i] = fcos0 + fcos1 + fcos2 + fcos3;

		trans_tableptr += rowsize;
	}

	if((m % 2) == 1)
		/* if m odd, no need to do last row - all zeroes */
		fcos[bw-1] = 0.0;
	else
	{
		rowsize = Transpose_RowSize(bw-1, m, bw);

		/* need to point to correct starting value of Legendre series */
		assoc_offset = assoc_legendre_series + ((bw - 1 - m) + (m % 2))*2;

		fcos0 = 0.0; fcos1 = 0.0; fcos2 = 0.0; fcos3 = 0.0;

		for(j = 0; j < rowsize % 4; ++j)
			fcos0 += assoc_offset[2*2*j] * trans_tableptr[j];

		for ( ; j < rowsize; j += 4){
			fcos0 += assoc_offset[2*2*j] * trans_tableptr[j];
			fcos1 += assoc_offset[2*2*(j+1)] * trans_tableptr[j+1];
			fcos2 += assoc_offset[2*2*(j+2)] * trans_tableptr[j+2];
			fcos3 += assoc_offset[2*2*(j+3)] * trans_tableptr[j+3];
		}

		fcos[bw - 1] = fcos0 + fcos1 + fcos2 + fcos3;

		trans_tableptr += rowsize;
	}


	/* now we have the cosine series for the result,
	so now evaluate the cosine series at 2*bw Chebyshev nodes 
	*/

	ExpIFCT(fcos, result, scratchpad, 2*bw, bw, 1);

	/* if m is odd, then need to multiply by sin(x) at Chebyshev nodes */

	if ((m % 2) == 1)
	{
		n = 2*bw;

		for (j=0; j<n; j++)
			result[j] *= sin_values[j];
	}

	trans_tableptr = p;

	/* amscray */

}

void InvSemiNaiveReduced_fftw(float *assoc_legendre_series, 
							  int bw, 
							  int m, 
							  float *result, 
							  float *trans_cos_pml_table, 
							  float *sin_values,
							  float *workspace)
{
	float *fcos; /* buffer for cosine series rep of result */
	float *trans_tableptr;
	float *scratchpad;
	float *assoc_offset;
	int i, j, n, rowsize;

	float *p;
	float fcos0, fcos1, fcos2, fcos3;

	fcos = workspace; /* needs bw */
	scratchpad = fcos + bw; /* needs (4*bw) */

	/* total workspace is 5*bw */

	trans_tableptr = trans_cos_pml_table;
	p = trans_cos_pml_table;

	/* main loop - compute each value of fcos

	Note that all zeroes have been stripped out of the
	trans_cos_pml_table, so indexing is somewhat complicated.
	*/

	for (i=0; i<=m; i++)
	{
		rowsize = Transpose_RowSize(i, m, bw);

		/* need to point to correct starting value of Legendre series */
		assoc_offset = assoc_legendre_series + (i % 2)*2;

		fcos0 = 0.0 ; fcos1 = 0.0; fcos2 = 0.0; fcos3 = 0.0;

		for (j = 0; j < rowsize % 4; ++j)
			fcos0 += assoc_offset[2*2*j] * trans_tableptr[j];

		for ( ; j < rowsize; j += 4){
			fcos0 += assoc_offset[2*2*j] * trans_tableptr[j];
			fcos1 += assoc_offset[2*2*(j+1)] * trans_tableptr[j+1];
			fcos2 += assoc_offset[2*2*(j+2)] * trans_tableptr[j+2];
			fcos3 += assoc_offset[2*2*(j+3)] * trans_tableptr[j+3];
		}

		fcos[i] = fcos0 + fcos1 + fcos2 + fcos3;

		trans_tableptr += rowsize;      
	}

	for (i=m+1; i<bw-1; i++)
	{

		rowsize = Transpose_RowSize(i, m, bw);

		/* need to point to correct starting value of Legendre series */
		assoc_offset = assoc_legendre_series + ((i - m) + (m % 2))*2;

		fcos0 = 0.0 ; fcos1 = 0.0; fcos2 = 0.0; fcos3 = 0.0;

		for (j = 0; j < rowsize % 4; ++j)
			fcos0 += assoc_offset[2*2*j] * trans_tableptr[j];

		for ( ; j < rowsize; j += 4){
			fcos0 += assoc_offset[2*2*j] * trans_tableptr[j];
			fcos1 += assoc_offset[2*2*(j+1)] * trans_tableptr[j+1];
			fcos2 += assoc_offset[2*2*(j+2)] * trans_tableptr[j+2];
			fcos3 += assoc_offset[2*2*(j+3)] * trans_tableptr[j+3];
		}
		fcos[i] = fcos0 + fcos1 + fcos2 + fcos3;

		trans_tableptr += rowsize;
	}

	if((m % 2) == 1)
		/* if m odd, no need to do last row - all zeroes */
		fcos[bw-1] = 0.0;
	else
	{
		rowsize = Transpose_RowSize(bw-1, m, bw);

		/* need to point to correct starting value of Legendre series */
		assoc_offset = assoc_legendre_series + ((bw - 1 - m) + (m % 2))*2;

		fcos0 = 0.0; fcos1 = 0.0; fcos2 = 0.0; fcos3 = 0.0;

		for(j = 0; j < rowsize % 4; ++j)
			fcos0 += assoc_offset[2*2*j] * trans_tableptr[j];

		for ( ; j < rowsize; j += 4){
			fcos0 += assoc_offset[2*2*j] * trans_tableptr[j];
			fcos1 += assoc_offset[2*2*(j+1)] * trans_tableptr[j+1];
			fcos2 += assoc_offset[2*2*(j+2)] * trans_tableptr[j+2];
			fcos3 += assoc_offset[2*2*(j+3)] * trans_tableptr[j+3];
		}

		fcos[bw - 1] = fcos0 + fcos1 + fcos2 + fcos3;

		trans_tableptr += rowsize;
	}


	/* now we have the cosine series for the result,
	so now evaluate the cosine series at 2*bw Chebyshev nodes 
	*/

	ExpIFCT(fcos, result, scratchpad, 2*bw, bw, 1);

	/* if m is odd, then need to multiply by sin(x) at Chebyshev nodes */

	if ((m % 2) == 1)
	{
		n = 2*bw;

		for (j=0; j<n; j++)
			result[j] *= sin_values[j];
	}

	trans_tableptr = p;

	/* amscray */

}
