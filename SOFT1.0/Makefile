##                 SOFT: SO(3) Fourier transform code
##   
##                 Version 1.0
##
##
##    Peter Kostelec, Dan Rockmore
##    {geelong,rockmore}@cs.dartmouth.edu
##   
##    Contact: Peter Kostelec
##             geelong@cs.dartmouth.edu
##   
##   
##    Copyright 2003 Peter Kostelec, Dan Rockmore
##   
##   
##      This program is free software; you can redistribute it and/or modify
##      it under the terms of the GNU General Public License as published by
##      the Free Software Foundation; either version 2 of the License, or
##      (at your option) any later version.
##   
##      This program is distributed in the hope that it will be useful,
##      but WITHOUT ANY WARRANTY; without even the implied warranty of
##      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##      GNU General Public License for more details.
##   
##      You should have received a copy of the GNU General Public License
##      along with this program; if not, write to the Free Software
##      Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
##   
##   
##    Commercial use is absolutely prohibited.
##   
##    See the accompanying LICENSE file for details.
  

CC = cc

# define for fftw, otherwise leave blank
## FFTWINC = -I/net/misc/geelong/local/linux/include
## FFTWLIB = -L/net/misc/geelong/local/linux/lib -lfftw3

## FFTWINC = -I/net/misc/geelong/local/athlon/include
## FFTWLIB = -L/net/misc/geelong/local/athlon/lib -lfftw3

FFTWINC =
FFTWLIB =


## define WALLCLOCK if want walltime and not cputime
## CFLAGS = -Wall -g -DWALLCLOCK
## CFLAGS = -Wall -g
## CFLAGS = -O3 -pg
## CFLAGS = -O3 -pg ${FFTWINC}
CFLAGS = -O3 ${FFTWINC}
## CFLAGS = -Wall -g ${FFTWINC}

LDFLAGS = -lm


# generating Wigner little-d functions
GENWIGSRC = utils_so3.c makeWigner.c
GENWIGOBJ = utils_so3.o makeWigner.o

# naive Wigner-d transforms
WIGSRC = $(GENWIGSRC) weights.c wignerTransforms.c
WIGOBJ = $(GENWIGOBJ) weights.o wignerTransforms.o

# naive Wigner-d transforms, but using symmetries
WIGSYMSRC = $(WIGSRC) \
	wignerTransforms_sym.c
WIGSYMOBJ = $(WIGOBJ) \
	wignerTransforms_sym.o

# naive Wigner-d transforms, but using fftw (and symmetries)
WIGFFTWSRC = $(GENWIGSRC) weights.c \
	utils_vec_cx.c wignerTransforms_fftw.c
WIGFFTWOBJ = $(GENWIGSRC) weights.o \
	utils_vec_cx.o wignerTransforms_fftw.o

#
# SO(3) "soft" transforms
#

# soft transform
SOFTSRC =  indextables.c permroots.c FFTcode.c fft_grids_so3.c \
	$(WIGSRC) soft.c
SOFTOBJ =  indextables.o permroots.o FFTcode.o fft_grids_so3.o \
	$(WIGOBJ) soft.o

# soft transform using symmetries of Wigner little-d's
SOFTSYMSRC =  indextables.c permroots.c FFTcode.c fft_grids_so3.c \
	$(WIGSYMSRC) soft_sym.c
SOFTSYMOBJ =  indextables.o permroots.o FFTcode.o fft_grids_so3.o \
	$(WIGSYMOBJ) soft_sym.o

# soft transform using fftw (and symmetries)
SOFTFFTWSRC = $(WIGFFTWSRC) soft_fftw.c
SOFTFFTWOBJ = $(WIGFFTWOBJ) soft_fftw.o

# soft transform using fftw (and symmetries, and FFT - FOR CORRELATION)
SOFTFFTW2SRC = $(WIGFFTWSRC) soft_fftw.c indextables.c \
	permroots.c FFTcode.c fft_grids_so3.c
SOFTFFTW2OBJ = $(WIGFFTWOBJ) soft_fftw.o indextables.o \
	permroots.o FFTcode.o fft_grids_so3.o \

# soft transform using fftw (and symmetries) and PRECOMPUTED wigners
SOFTFFTWPCSRC = $(WIGFFTWSRC) soft_fftw_pc.c
SOFTFFTWPCOBJ = $(WIGFFTWOBJ) soft_fftw_pc.o

# soft transform using fftw (and symmetries) - WILL WRITE
# OVER FUNCTION INPUTS!!!
SOFTFFTWWOSRC = $(WIGFFTWSRC) soft_fftw_wo.c
SOFTFFTWWOOBJ = $(WIGFFTWOBJ) soft_fftw_wo.o



#
# correlation routines (needs files from SpharmonicKit which are
# included within this distribution)
#

# correlation using Wigner-d symmetries - assumes you will already have
# the S^2 coefficients using other means
CORSRC = $(SOFTSYMSRC) primitive_FST.c so3_correlate_sym.c
COROBJ = $(SOFTSYMOBJ) primitive_FST.o so3_correlate_sym.o

# correlation using Wigner-d symmetries - start from scratch
COR2SRC = $(SOFTSYMSRC) primitive_FST.c so3_correlate_sym.c \
	OURperms.c OURmods.c newFCT.c cospmls.c seminaive.c oddweights.c \
	primitive.c FST_semi_memo.c fft_grids.c naive_synthesis.c
COR2OBJ = $(SOFTSYMOBJ) primitive_FST.o so3_correlate_sym.o \
	OURperms.o OURmods.o newFCT.o cospmls.o seminaive.o oddweights.o \
	primitive.o FST_semi_memo.o fft_grids.o naive_synthesis.o

# correlation using Wigner-d symmetries and a bit of fftw - start from scratch
COR2FFTWSRC = $(SOFTFFTW2SRC) primitive_FST.c so3_correlate_fftw.c \
	OURperms.c OURmods.c newFCT.c cospmls.c seminaive.c oddweights.c \
	primitive.c FST_semi_memo.c fft_grids.c naive_synthesis.c
COR2FFTWOBJ = $(SOFTFFTW2OBJ) primitive_FST.o so3_correlate_fftw.o \
	OURperms.o OURmods.o newFCT.o cospmls.o seminaive.o oddweights.o \
	primitive.o FST_semi_memo.o fft_grids.o naive_synthesis.o


#
# rotating S^2 functions via massaging their spherical coefficients
# using the Wigner-d *and* Wigner-D functions; depends on code
# from SpharmonicKit which is included within this distribution
#

# first, to construct the S^2 coefficients using semi-naive transform
SEMISRC = OURperms.c OURmods.c newFCT.c weights.c \
	cospmls.c seminaive.c oddweights.c \
	primitive.c

SEMIOBJ = OURperms.o OURmods.o newFCT.o weights.o \
	cospmls.o seminaive.o oddweights.o \
	primitive.o

# seminaive spherical transform and convolution
FSTSEMISRC = $(SEMISRC) naive_synthesis.c \
	primitive_FST.c permroots.c indextables.c FFTcode.c \
	fft_grids.c FST_semi_memo.c

FSTSEMIOBJ = $(SEMIOBJ) naive_synthesis.o \
	primitive_FST.o permroots.o indextables.o FFTcode.o \
	fft_grids.o FST_semi_memo.o

#
# now code to massage the coefficients
#

# test rotation via euler angles
ROTATESO3SRC = $(FSTSEMISRC) rotate_so3.c
ROTATESO3OBJ = $(FSTSEMIOBJ) rotate_so3.o

# test rotation via euler angles - a little friendlier on
# the memory
ROTATESO3MEMSRC = $(FSTSEMISRC) rotate_so3_mem.c
ROTATESO3MEMOBJ = $(FSTSEMIOBJ) rotate_so3_mem.o


# for makedepend

SRC = FFTcode.c FST_semi_memo.c OURmods.c OURperms.c \
	cospmls.c csecond.c fft_grids.c fft_grids_so3.c \
	indextables.c makeWigner.c naive_synthesis.c newFCT.c \
	oddweights.c permroots.c primitive.c primitive_FST.c \
	rotate_so3.c rotate_so3_mem.c seminaive.c so3_correlate_sym.c \
	so3_correlate_fftw.c \
	soft.c soft_fftw.c soft_fftw_pc.c soft_fftw_wo.c soft_sym.c \
	test_s2_rotate.c test_s2_rotate_mem.c test_Wigner_Analysis.c \
	test_Wigner_Naive.c test_Wigner_Synthesis.c test_genWig.c \
	test_soft.c test_soft_fftw.c test_soft_fftw_correlate2.c \
	test_soft_fftw_pc.c test_soft_fftw_wo.c test_soft_for.c \
	test_soft_inv.c test_soft_sym.c test_soft_sym_correlate.c \
	test_soft_sym_correlate2.c test_soft_sym_for.c \
	test_soft_sym_inv.c test_wigSpec.c utils_so3.c \
	utils_vec_cx.c weights.c wignerTransforms.c \
	wignerTransforms_fftw.c wignerTransforms_sym.c

###################################################################
##################################################################
######
######              things that can be made
######
##################################################################
##################################################################

# things that can be made

depend:
	makedepend -I. ${FFTWINC} $(SRC)

clean:
	rm *.o

allnonfftw:
	make \
	test_wigSpec \
	test_genWig \
	test_Wigner_Analysis \
	test_Wigner_Synthesis \
	test_Wigner_Naive \
	test_soft \
	test_soft_sym \
	test_soft_for \
	test_soft_sym_for \
	test_soft_inv \
	test_soft_sym_inv \
	test_soft_sym_correlate \
	test_soft_sym_correlate2 \
	test_s2_rotate \
	test_s2_rotate_mem

all:
	make \
	allnonfftw \
	test_soft_fftw \
	test_soft_fftw_pc \
	test_soft_fftw_wo \
	test_soft_fftw_correlate2


# now the make definitions for the individual executables

test_wigSpec: $(GENWIGOBJ) test_wigSpec.o
	$(CC) $(CFLAGS) $(GENWIGOBJ) test_wigSpec.o \
	$(LDFLAGS) -o test_wigSpec

test_genWig: $(GENWIGOBJ) test_genWig.o
	$(CC) $(CFLAGS) $(GENWIGOBJ) test_genWig.o \
	$(LDFLAGS) -o test_genWig

test_Wigner_Analysis: $(WIGOBJ) csecond.o test_Wigner_Analysis.o
	$(CC) $(CFLAGS) $(WIGOBJ) csecond.o test_Wigner_Analysis.o \
	$(LDFLAGS) -o test_Wigner_Analysis

test_Wigner_Synthesis: $(WIGOBJ) csecond.o test_Wigner_Synthesis.o
	$(CC) $(CFLAGS) $(WIGOBJ) csecond.o test_Wigner_Synthesis.o \
	$(LDFLAGS) -o test_Wigner_Synthesis

test_Wigner_Naive: $(WIGOBJ) csecond.o test_Wigner_Naive.o
	$(CC) $(CFLAGS) $(WIGOBJ) csecond.o test_Wigner_Naive.o \
	$(LDFLAGS) -o test_Wigner_Naive

test_soft: $(SOFTOBJ) csecond.o test_soft.o
	$(CC) $(CFLAGS) $(SOFTOBJ) csecond.o test_soft.o \
	$(LDFLAGS) -o test_soft

test_soft_fftw: $(SOFTFFTWOBJ) csecond.o test_soft_fftw.o
	$(CC) $(CFLAGS) $(SOFTFFTWOBJ) csecond.o test_soft_fftw.o \
	${FFTWLIB} $(LDFLAGS) -o test_soft_fftw

test_soft_fftw_pc: $(SOFTFFTWPCOBJ) csecond.o test_soft_fftw_pc.o
	$(CC) $(CFLAGS) $(SOFTFFTWPCOBJ) csecond.o test_soft_fftw_pc.o \
	${FFTWLIB} $(LDFLAGS) -o test_soft_fftw_pc

test_soft_fftw_wo: $(SOFTFFTWWOOBJ) csecond.o test_soft_fftw_wo.o
	$(CC) $(CFLAGS) $(SOFTFFTWWOOBJ) csecond.o test_soft_fftw_wo.o \
	${FFTWLIB} $(LDFLAGS) -o test_soft_fftw_wo

test_soft_sym: $(SOFTSYMOBJ) csecond.o test_soft_sym.o
	$(CC) $(CFLAGS) $(SOFTSYMOBJ) csecond.o test_soft_sym.o \
	$(LDFLAGS) -o test_soft_sym

test_soft_for: $(SOFTOBJ) csecond.o test_soft_for.o
	$(CC) $(CFLAGS) $(SOFTOBJ) csecond.o test_soft_for.o \
	$(LDFLAGS) -o test_soft_for

test_soft_sym_for: $(SOFTSYMOBJ) csecond.o test_soft_sym_for.o
	$(CC) $(CFLAGS) $(SOFTSYMOBJ) csecond.o test_soft_sym_for.o \
	$(LDFLAGS) -o test_soft_sym_for

test_soft_inv: $(SOFTOBJ) csecond.o test_soft_inv.o
	$(CC) $(CFLAGS) $(SOFTOBJ) csecond.o test_soft_inv.o \
	$(LDFLAGS) -o test_soft_inv

test_soft_sym_inv: $(SOFTSYMOBJ) csecond.o test_soft_sym_inv.o
	$(CC) $(CFLAGS) $(SOFTSYMOBJ) csecond.o test_soft_sym_inv.o \
	$(LDFLAGS) -o test_soft_sym_inv

test_soft_sym_correlate: $(COROBJ) csecond.o test_soft_sym_correlate.o
	$(CC) $(CFLAGS) $(COROBJ) csecond.o test_soft_sym_correlate.o \
	$(LDFLAGS) -o test_soft_sym_correlate

test_soft_sym_correlate2: $(COR2OBJ) csecond.o test_soft_sym_correlate2.o
	$(CC) $(CFLAGS) $(COR2OBJ) csecond.o test_soft_sym_correlate2.o \
	$(LDFLAGS) -o test_soft_sym_correlate2

test_soft_fftw_correlate2: $(COR2FFTWOBJ) csecond.o test_soft_fftw_correlate2.o
	$(CC) $(CFLAGS) $(COR2FFTWOBJ) csecond.o test_soft_fftw_correlate2.o \
	${FFTWLIB} $(LDFLAGS) -o test_soft_fftw_correlate2

test_s2_rotate: $(ROTATESO3OBJ) csecond.o test_s2_rotate.o
	$(CC) $(CFLAGS) $(ROTATESO3OBJ) csecond.o test_s2_rotate.o \
	$(LDFLAGS) -o test_s2_rotate

test_s2_rotate_mem: $(ROTATESO3MEMOBJ) csecond.o test_s2_rotate_mem.o
	$(CC) $(CFLAGS) $(ROTATESO3MEMOBJ) csecond.o test_s2_rotate_mem.o \
	$(LDFLAGS) -o test_s2_rotate_mem


# and now for LOTS OF dependencies ...

# DO NOT DELETE THIS LINE -- make depend depends on it.

FFTcode.o: indextables.h permroots.h
FST_semi_memo.o: cospmls.h fft_grids.h
FST_semi_memo.o: naive_synthesis.h primitive.h primitive_FST.h seminaive.h
cospmls.o: newFCT.h primitive.h
fft_grids.o: primitive_FST.h FFTcode.h
fft_grids_so3.o: FFTcode.h
makeWigner.o: utils_so3.h weights.h
naive_synthesis.o: csecond.h weights.h
newFCT.o: OURperms.h OURmods.h newFCT.h
rotate_so3.o: complex.h primitive_FST.h rotate_so3.h FST_semi_memo.h
rotate_so3_mem.o: complex.h primitive_FST.h
rotate_so3_mem.o: rotate_so3_mem.h FST_semi_memo.h
seminaive.o: cospmls.h csecond.h
seminaive.o: newFCT.h oddweights.h weights.h
so3_correlate_sym.o: complex.h primitive_FST.h
so3_correlate_sym.o: utils_so3.h so3_correlate_sym.h
so3_correlate_fftw.o: complex.h primitive_FST.h
so3_correlate_fftw.o: utils_so3.h so3_correlate_fftw.h
soft.o: utils_so3.h fft_grids_so3.h makeWigner.h wignerTransforms.h
soft_fftw.o: complex.h utils_so3.h utils_vec_cx.h fft_grids_so3.h
soft_fftw.o: makeWigner.h wignerTransforms_fftw.h
soft_fftw_pc.o: complex.h utils_so3.h utils_vec_cx.h fft_grids_so3.h
soft_fftw_pc.o: makeWigner.h wignerTransforms_fftw.h
soft_fftw_wo.o: complex.h utils_vec_cx.h utils_so3.h fft_grids_so3.h
soft_fftw_wo.o: makeWigner.h wignerTransforms_fftw.h
soft_sym.o: utils_so3.h fft_grids_so3.h
soft_sym.o: makeWigner.h wignerTransforms.h wignerTransforms_sym.h
test_s2_rotate.o: complex.h FST_semi_memo.h csecond.h cospmls.h
test_s2_rotate.o: primitive_FST.h seminaive.h rotate_so3.h
test_s2_rotate_mem.o: FST_semi_memo.h csecond.h cospmls.h primitive_FST.h
test_s2_rotate_mem.o: seminaive.h rotate_so3_mem.h complex.h
test_Wigner_Analysis.o: utils_so3.h makeWigner.h wignerTransforms.h
test_Wigner_Naive.o: utils_so3.h makeWigner.h
test_Wigner_Naive.o: wignerTransforms.h csecond.h
test_Wigner_Synthesis.o: utils_so3.h makeWigner.h wignerTransforms.h
test_genWig.o: utils_so3.h makeWigner.h
test_soft.o: utils_so3.h soft.h csecond.h
test_soft_fftw.o: complex.h utils_so3.h soft_fftw.h csecond.h
test_soft_fftw_correlate2.o: complex.h csecond.h so3_correlate_fftw.h soft_fftw.h
test_soft_fftw_correlate2.o: FST_semi_memo.h cospmls.h primitive_FST.h
test_soft_fftw_correlate2.o: seminaive.h
test_soft_fftw_pc.o: complex.h makeWigner.h utils_so3.h
test_soft_fftw_pc.o: soft_fftw_pc.h csecond.h
test_soft_fftw_wo.o: complex.h utils_so3.h soft_fftw_wo.h csecond.h
test_soft_for.o: utils_so3.h soft.h csecond.h
test_soft_inv.o: utils_so3.h soft.h csecond.h
test_soft_sym.o: utils_so3.h soft_sym.h csecond.h
test_soft_sym_correlate.o: complex.h csecond.h so3_correlate_sym.h soft_sym.h
test_soft_sym_correlate2.o: complex.h so3_correlate_sym.h soft_sym.h FST_semi_memo.h
test_soft_sym_correlate2.o: csecond.h cospmls.h primitive_FST.h seminaive.h
test_soft_sym_for.o: utils_so3.h soft_sym.h csecond.h
test_soft_sym_inv.o: utils_so3.h soft_sym.h csecond.h
test_wigSpec.o: utils_so3.h makeWigner.h
utils_so3.o: utils_so3.h
utils_vec_cx.o: complex.h utils_vec_cx.h
wignerTransforms.o: utils_so3.h weights.h
wignerTransforms_fftw.o: complex.h utils_so3.h utils_vec_cx.h weights.h
wignerTransforms_sym.o: utils_so3.h weights.h
