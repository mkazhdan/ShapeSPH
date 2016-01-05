##############################################################################
##############################################################################
##  
##                           S2kit 1.0
##
##          A lite version of Spherical Harmonic Transform Kit
##
##  Peter Kostelec, Dan Rockmore
##  {geelong,rockmore}@cs.dartmouth.edu
##  
##  Contact: Peter Kostelec
##           geelong@cs.dartmouth.edu
##  
##  Copyright 2004 Peter Kostelec, Dan Rockmore
##
##  This file is part of S2kit.
##
##  S2kit is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 2 of the License, or
##  (at your option) any later version.
##
##  S2kit is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with S2kit; if not, write to the Free Software
##  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
##
##  See the accompanying LICENSE file for details.
##
##
############################################################################
############################################################################

CC = cc

# define WALLCLOCK on the CFLAGS line if want to time
# walltime and not cpu time (cpu time is default); also
# define optimization flags suitable for your platform

# define for fftw
## FFTWINC = -I/net/misc/geelong/local/linux/include
## FFTWLIB = -L/net/misc/geelong/local/linux/lib -lfftw3
FFTWINC =
FFTWLIB =


## CFLAGS = -O3 -pg ${FFTWINC}
## CFLAGS = -g -Wall ${FFTWINC}
CFLAGS = -O3 ${FFTWINC}

LDFLAGS = -lm

# for naive
NAIVESRC = primitive.c pmls.c naive_synthesis.c \
	makeweights.c csecond.c

NAIVEOBJ = primitive.c pmls.o naive_synthesis.o \
	makeweights.o csecond.o

# for semi-naive
SEMISRC = pmls.c cospmls.c seminaive.c \
	csecond.c primitive.c makeweights.c

SEMIOBJ = pmls.o cospmls.o seminaive.o \
	csecond.o primitive.o makeweights.o

# seminaive spherical transform and convolution
FSTSEMISRC = $(SEMISRC) naive_synthesis.c

FSTSEMIOBJ = $(SEMIOBJ) naive_synthesis.o

CONVSEMISRC = $(FSTSEMISRC)

CONVSEMIOBJ = $(FSTSEMIOBJ)

ALLSRC = FST_semi_fly.c FST_semi_memo.c cospmls.c csecond.c \
	makeweights.c naive_synthesis.c pmls.c primitive.c \
	seminaive.c test_conv_semi_fly.c test_conv_semi_memo.c \
	test_naive.c test_s2_semi_fly.c test_s2_semi_memo.c \
	test_s2_semi_memo_for.c test_s2_semi_memo_inv.c \
	test_semi.c


###################################################################
##################################################################
######
######              things that can be made
######
##################################################################
##################################################################

# first some shortcuts

all:
	make \
	legendre \
	sphere


legendre:
	make \
	test_naive \
	test_semi

sphere:
	make \
	test_s2_semi_memo \
	test_s2_semi_memo_for \
	test_s2_semi_memo_inv \
	test_s2_semi_fly \
	test_conv_semi_memo \
	test_conv_semi_fly

depend:
	makedepend ${FFTWINC} ${ALLSRC}

clean: 
	rm *.o

# now the make definitions for the individual executables

test_naive: $(NAIVEOBJ) test_naive.o
	$(CC) $(CFLAGS) $(NAIVEOBJ) test_naive.o \
	$(LDFLAGS) -o test_naive

test_semi: $(SEMIOBJ) test_semi.o
	$(CC) $(CFLAGS) $(SEMIOBJ) test_semi.o \
	${FFTWLIB} $(LDFLAGS) -o test_semi

test_s2_semi_memo: $(FSTSEMIOBJ) FST_semi_memo.o test_s2_semi_memo.o
	$(CC) $(CFLAGS) $(FSTSEMIOBJ) FST_semi_memo.o test_s2_semi_memo.o \
	${FFTWLIB} $(LDFLAGS) -o test_s2_semi_memo

test_s2_semi_memo_for: $(FSTSEMIOBJ) FST_semi_memo.o test_s2_semi_memo_for.o
	$(CC) $(CFLAGS) $(FSTSEMIOBJ) FST_semi_memo.o test_s2_semi_memo_for.o \
	${FFTWLIB} $(LDFLAGS) -o test_s2_semi_memo_for

test_s2_semi_memo_inv: $(FSTSEMIOBJ) FST_semi_memo.o test_s2_semi_memo_inv.o
	$(CC) $(CFLAGS) $(FSTSEMIOBJ) FST_semi_memo.o test_s2_semi_memo_inv.o \
	${FFTWLIB} $(LDFLAGS) -o test_s2_semi_memo_inv

test_s2_semi_fly: $(FSTSEMIOBJ) FST_semi_fly.o test_s2_semi_fly.o
	$(CC) $(CFLAGS) $(FSTSEMIOBJ) FST_semi_fly.o test_s2_semi_fly.o \
	${FFTWLIB} $(LDFLAGS) -o test_s2_semi_fly

test_conv_semi_memo: $(CONVSEMIOBJ) FST_semi_memo.o test_conv_semi_memo.o
	$(CC) $(CFLAGS) $(CONVSEMIOBJ) FST_semi_memo.o test_conv_semi_memo.o \
	${FFTWLIB} $(LDFLAGS) -o test_conv_semi_memo

test_conv_semi_fly: $(CONVSEMIOBJ) FST_semi_fly.o test_conv_semi_fly.o
	$(CC) $(CFLAGS) $(CONVSEMIOBJ) FST_semi_fly.o test_conv_semi_fly.o \
	${FFTWLIB} $(LDFLAGS) -o test_conv_semi_fly


# and now for LOTS OF dependencies ...

# DO NOT DELETE THIS LINE -- make depend depends on it.

FST_semi_fly.o: primitive.h makeweights.h pmls.h cospmls.h naive_synthesis.h
FST_semi_fly.o: seminaive.h FST_semi_fly.h

FST_semi_memo.o: makeweights.h FST_semi_memo.h
FST_semi_memo.o: cospmls.h primitive.h naive_synthesis.h seminaive.h

cospmls.o: primitive.h pmls.h

pmls.o: primitive.h

seminaive.o: cospmls.h

test_conv_semi_fly.o: FST_semi_fly.h

test_conv_semi_memo.o: FST_semi_memo.h cospmls.h

test_naive.o: pmls.h makeweights.h naive_synthesis.h csecond.h

test_s2_semi_fly.o: makeweights.h FST_semi_fly.h csecond.h

test_s2_semi_memo.o: makeweights.h cospmls.h FST_semi_memo.h csecond.h

test_s2_semi_memo_for.o: makeweights.h cospmls.h FST_semi_memo.h csecond.h

test_s2_semi_memo_inv.o: makeweights.h cospmls.h FST_semi_memo.h csecond.h

test_semi.o: makeweights.h cospmls.h primitive.h seminaive.h csecond.h
