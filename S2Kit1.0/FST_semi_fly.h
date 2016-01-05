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



/* external interface for FST_semi_fly.c */

#ifndef _FSTSEMI_FLY_H
#define _FSTSEMI_FLY_H

#define compmult(a,b,c,d,e,f) (e) = ((a)*(c))-((b)*(d)); (f) = ((a)*(d))+((b)*(c))

extern int seanindex( int ,
		      int ,
		      int );

extern void TransMult( double *, double *,
		       double *, double *,
		       double *, double *,
		       int );

extern void FST_semi_fly( double *, double *,
			  double *, double *,
			  int ,
			  double *,
			  int ,
			  int ,
			  fftw_plan *,
			  fftw_plan *,
			  double * );

extern void InvFST_semi_fly( double *, double *, 
			     double *, double *,
			     int , 
			     double *,
			     int ,
			     int ,
			     fftw_plan *,
			     fftw_plan * );

extern void FZT_semi_fly( double *, double *,
			  double *, double *,
			  int ,
			  double *,
			  int ,
			  fftw_plan *,
			  double * );

extern void Conv2Sphere_semi_fly( double *, double *,
				  double *, double *,
				  double *, double *,
				  int ,
				  double *);

#endif /* _FSTSEMI_FLY_H */
