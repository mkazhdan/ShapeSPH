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



/* source code for generating cosine transforms 
   of Pml and Gml functions */

#include <math.h>
#include <string.h>   /* to declare memcpy */
#include <stdlib.h>
#include <stdio.h>

#include "primitive.h"

/************************************************************************/
/* generate all of the Pmls for a specified value of m.  

   storeplm points to a double array of size 2 * bw * (bw - m);

   Workspace needs to be
   16 * bw 

   P(m,l,j) respresents the associated Legendre function P_l^m
   evaluated at the j-th Chebyshev point (for the bandwidth bw)
   Cos((2 * j + 1) * PI / (2 * bw)).

   The array is placed in storeplm as follows:

   P(m,m,0)    P(m,m,1)  ... P(m,m,2*bw-1)
   P(m,m+1,0)  P(m,m+1,1)... P(m,m+1,2*bw-1)
   P(m,m+2,0)  P(m,m+2,1)... P(m,m+2,2*bw-1)
   ...
   P(m,bw-1,0)   P(m,bw-1,1) ... P(m,bw-1,2*bw-1)

   This array will eventually be used by the naive transform algorithm.
   This function will precompute the arrays necessary for the algorithm.
*/

void PmlTableGen(int bw,
		 int m,
		 double *storeplm,
		 double *workspace)
{
  double *prev, *prevprev;
  double *temp1, *temp2, *temp3, *temp4;
  double *x_i, *eval_args;
  int i;
  
  prevprev = workspace;
  prev = prevprev + (2*bw);
  temp1 = prev + (2*bw);
  temp2 = temp1 + (2*bw);
  temp3 = temp2 + (2*bw);
  temp4 = temp3 + (2*bw);
  x_i = temp4 + (2*bw);
  eval_args = x_i + (2*bw);
  

  /* get the evaluation nodes */
  EvalPts(2*bw,x_i);
  ArcCosEvalPts(2*bw,eval_args);
  
  /* set initial values of first two Pmls */
  for (i=0; i<2*bw; i++) 
    prevprev[i] = 0.0;
  if (m == 0)
    for (i=0; i<2*bw; i++)
      prev[i] = 0.707106781186547 ;
  else 
    Pmm_L2(m, eval_args, 2*bw, prev);

  memcpy(storeplm, prev, sizeof(double) * 2 * bw);

  for(i = 0; i < bw - m - 1; i++)
    {
      vec_mul(L2_cn(m,m+i),prevprev,temp1,2*bw);
      vec_pt_mul(prev, x_i, temp2, 2*bw);
      vec_mul(L2_an(m,m+i), temp2, temp3, 2*bw);
      vec_add(temp3, temp1, temp4, 2*bw); /* temp4 now contains P(m,m+i+1) */
      
      storeplm += (2 * bw);
      memcpy(storeplm, temp4, sizeof(double) * 2 * bw);
      memcpy(prevprev, prev, sizeof(double) * 2 * bw);
      memcpy(prev, temp4, sizeof(double) * 2 * bw);
    }
}


