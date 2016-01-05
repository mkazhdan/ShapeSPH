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



/* external interface for FST_semi_memo.c */

#ifndef _FSTSEMI_MEMO_FFTW_H
#define _FSTSEMI_MEMO_FFTW_H

// Strictly real transforms (Misha Added)
extern void FST_semi_memo_fftw		(double *, double *, double *, int, double **, double *);
extern void FST_semi_memo_fftw		(double *, fftw_complex *, int, double **, double *);
extern void FST_semi_memo_fftw		(float *, float *, float *, int, float **, float *);
extern void FST_semi_memo_fftw		(float *, fftwf_complex *, int, float **, float *);
extern void InvFST_semi_memo_fftw	(double *, double *, double *, int, double **, double *);
extern void InvFST_semi_memo_fftw	(fftw_complex *, double *, int, double **, double *);
extern void InvFST_semi_memo_fftw	(float *, float *, float *, int, float **, float *);
extern void InvFST_semi_memo_fftw	(fftwf_complex *, float *, int, float **, float *);
// Done misha added
#endif /* _FSTSEMI_MEMO_FFTW_H */
