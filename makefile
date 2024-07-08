#makefile for lrslib-072    2020.5.4 

# C compiler requires __int128 support for 128bit arithmetic (eg. gcc v. 4.6.0 or later for 128bit integer support) 
# otherwise use %make lrs64 to compile
# add -DTIMES if your compiler does *not* support ptimes()
# add -DSIGNALS if your compiler does *not* support <signal.h> <unistd.h>

#try uncommenting next line if cc is the default C compiler
CC = g++      # or gcc7
CPPFLAGS = -Wno-write-strings -fpermissive
default: lrs lrsgmp lrsnash checkpred inedel 

#choose line below instead if __int128 not supported
#default: lrs64 lrsgmp 

#make lrs               lrs,lrsgmp       hybrid and gmp versions 
#make lrs64             lrs,lrsgmp    compilers without 128 bit support
#make mplrs             mplrs,mplrsgmp hybrid and gmp versions,  make sure mpicc and an MPI library is installed
#make mplrs64           mplrs,mplrsgmp for compilers without 128 bit support

#make flint             lrs and mplrs with FLINT arithmetic 
#make single            makes lrs with various arithmetic packages (depending on compiler),lrsnash 
#make singlemplrs        makes mplrs with various arithmetic packages (depending on compiler)
#make allmp             uses native mp and long arithmetic
#make demo              various demo programs for lrslib     
#make alllrsnash        Nash equilibria for 2-person games: lrsnash (gmp), lrsnash1 (64bit), lrsnash2 (128bit), 2nash, setupnash
#make clean             removes binaries                                      

#INCLUDEDIR = /usr/include
#LIBDIR     = /usr/lib

#Kyoto machines usage
INCLUDEDIR = /usr/local/include
LIBDIR     = /usr/local/lib

CFLAGS     ?= -O3 -Wall 
#CFLAGS     = -g -Wall 

#use this if you want only output file contain data between begin/end lines
#CFLAGS     = -O3 -Wall -DLRS_QUIET

SHLIB_CFLAGS = -fPIC
mpicxx=mpicc


# for 32 bit machines

# BITS=
# MPLRSOBJ2=

# for 64 bit machines
BITS=-DB128
MPLRSOBJ2=lrslib2-mplrs.o lrslong2-mplrs.o


LRSOBJ=lrs.o lrslong1.o lrslong2.o lrslib1.o lrslib2.o lrslibgmp.o lrsgmp.o lrsdriver.o
LRSOBJMP=lrs.o lrslong1.o lrslong2.o lrslib1.o lrslib2.o lrslibmp.o lrsmp.o lrsdriver.o
MPLRSOBJ=lrslong1-mplrs.o lrslib1-mplrs.o ${MPLRSOBJ2} lrslibgmp-mplrs.o lrsgmp-mplrs.o lrsdriver-mplrs.o mplrs.o

LRSOBJ64=lrs64.o lrslong1.o lrslib1.o lrslibgmp.o lrsgmp.o lrsdriver.o
MPLRSOBJ64=lrslong1-mplrs.o lrslib1-mplrs.o lrslibgmp-mplrs.o lrsgmp-mplrs.o lrsdriver-mplrs.o mplrs64.o

lrs: ${LRSOBJ}
	$(CC) ${CFLAGS} -DMA ${BITS} -L${LIBDIR} -o lrs ${LRSOBJ} -lgmp
	$(CC) -O3 hvref.cc -o hvref
	ln -s -f lrs redund

lrsMP: ${LRSOBJMP}
	$(CC) ${CFLAGS} -DMA ${BITS} -o lrsMP ${LRSOBJMP} 
	$(CC) -O3 hvref.cc -o hvref
	ln -s -f lrs redund

lrs64: ${LRSOBJ64}
	$(CC) ${CFLAGS} -DMA -L${LIBDIR} -o lrs ${LRSOBJ64} -lgmp

lrs.o: lrs.cc
	$(CC) ${CFLAGS} -DMA ${BITS} -c -o lrs.o lrs.cc

lrs64.o: lrs.cc
	$(CC) ${CFLAGS} -DMA -c -o lrs64.o lrs.cc

lrslong1.o: lrslong.cc lrslong.h
	$(CC) ${CFLAGS} -DMA -DSAFE -DLRSLONG -c -o lrslong1.o lrslong.cc

lrslong2.o: lrslong.cc lrslong.h
	$(CC) ${CFLAGS} -DMA -DSAFE ${BITS} -DLRSLONG -c -o lrslong2.o lrslong.cc

lrslib1.o: lrslib.cc lrslib.h
	$(CC) ${CFLAGS} -DMA -DSAFE -DLRSLONG -c -o lrslib1.o lrslib.cc

lrslib2.o: lrslib.cc lrslib.h
	$(CC) ${CFLAGS} -DMA -DSAFE ${BITS} -DLRSLONG -c -o lrslib2.o lrslib.cc

lrslibgmp.o: lrslib.cc lrslib.h
	$(CC) ${CFLAGS}  -DMA -DGMP -I${INCLUDEDIR} -c -o lrslibgmp.o lrslib.cc

lrslibmp.o: lrslib.cc lrslib.h
	$(CC) ${CFLAGS}  -DMA -DMP -c -o lrslibmp.o lrslib.cc

lrsgmp.o: lrsgmp.cc lrsgmp.h
	$(CC) ${CFLAGS} -DMA -DGMP -I${INCLUDEDIR} -c -o lrsgmp.o lrsgmp.cc

lrsmp.o: lrsmp.cc lrsmp.h
	$(CC) ${CFLAGS} -DMA -DMP -c -o lrsmp.o lrsmp.cc

inedel: inedel.cc lrsgmp.h lrsgmp.cc
	$(CC) ${CFLAGS} -I${INCLUDEDIR} -L${LIBDIR} -DGMP -o inedel inedel.cc lrsgmp.cc -lgmp

checkpred: checkpred.cc lrsgmp.h lrsgmp.cc
	$(CC) $(CFLAGS) -I${INCLUDEDIR} -L${LIBDIR} -DGMP -o checkpred checkpred.cc lrsgmp.cc -lgmp

lrslong1-mplrs.o: lrslong.cc lrslong.h
	$(mpicxx) ${CFLAGS} -DTIMES -DSIGNALS -DMA -DSAFE -DLRSLONG -DPLRS -c -o lrslong1-mplrs.o lrslong.cc

lrslong2-mplrs.o: lrslong.cc lrslong.h
	$(mpicxx) ${CFLAGS} -DTIMES -DSIGNALS -DMA -DSAFE ${BITS} -DLRSLONG -DPLRS -c -o lrslong2-mplrs.o lrslong.cc

lrslib1-mplrs.o: lrslib.cc lrslib.h
	$(mpicxx) ${CFLAGS} -DTIMES -DSIGNALS -DMA -DSAFE -DLRSLONG -DPLRS -c -o lrslib1-mplrs.o lrslib.cc

lrslib2-mplrs.o: lrslib.cc lrslib.h
	$(mpicxx) ${CFLAGS} -DTIMES -DSIGNALS -DMA -DSAFE ${BITS} -DLRSLONG -DPLRS -c -o lrslib2-mplrs.o lrslib.cc

lrslibgmp-mplrs.o: lrslib.cc lrslib.h
	$(mpicxx) ${CFLAGS} -DMA -DTIMES -DSIGNALS -DGMP -DPLRS -I${INCLUDEDIR} -c -o lrslibgmp-mplrs.o lrslib.cc

lrsgmp-mplrs.o: lrsgmp.cc lrsgmp.h
	$(mpicxx) ${CFLAGS} -DMA -DTIMES -DSIGNALS -DGMP -DPLRS -I${INCLUDEDIR} -c -o lrsgmp-mplrs.o lrsgmp.cc

lrsdriver-mplrs.o: lrsdriver.cc lrsdriver.h lrslib.h
	$(mpicxx) $(CFLAGS) -c -o lrsdriver-mplrs.o lrsdriver.cc

mplrs.o: mplrs.cc mplrs.h lrslib.h lrsgmp.h
	$(mpicxx) ${CFLAGS} -I${INCLUDEDIR} -DMA -DPLRS -DTIMES ${BITS} -DSIGNALS -D_WITH_GETLINE -c -o mplrs.o mplrs.cc

mplrs64.o: mplrs.cc mplrs.h lrslib.h lrsgmp.h
	$(mpicxx) ${CFLAGS} -I${INCLUDEDIR} -DMA -DPLRS -DTIMES -DSIGNALS -D_WITH_GETLINE -c -o mplrs64.o mplrs.cc

mplrs: ${MPLRSOBJ} mplrsgmp
	$(mpicxx) ${CFLAGS} -DTIMES -DSIGNALS -D_WITH_GETLINE -DPLRS -DMA ${BITS} -L${LIBDIR} -o mplrs ${MPLRSOBJ} -lgmp

mplrs64: ${MPLRSOBJ64} mplrsgmp
	$(mpicxx) ${CFLAGS} -DTIMES -DSIGNALS -D_WITH_GETLINE -DPLRS -DMA -L${LIBDIR} -o mplrs ${MPLRSOBJ64} -lgmp

mplrsgmp: mplrs.cc mplrs.h lrslib.cc lrslib.h lrsgmp.cc lrsgmp.h lrsdriver.h lrsdriver.cc
	$(mpicxx) ${CFLAGS} -DTIMES -DSIGNALS -D_WITH_GETLINE -DPLRS -DGMP -I${INCLUDEDIR} mplrs.cc lrslib.cc lrsgmp.cc lrsdriver.cc -L${LIBDIR} -o mplrsgmp -lgmp

mplrs1: mplrs.cc mplrs.h lrslib.cc lrslib.h lrslong.cc lrslong.h lrsdriver.h lrsdriver.cc
	$(mpicxx) ${CFLAGS} -DTIMES -DSIGNALS -D_WITH_GETLINE -DPLRS -DSAFE -DLRSLONG mplrs.cc lrslib.cc lrslong.cc lrsdriver.cc -o mplrs1

mplrs2: mplrs.cc mplrs.h lrslib.cc lrslib.h lrslong.cc lrslong.h lrsdriver.h lrsdriver.cc
	$(mpicxx) ${CFLAGS} -DTIMES -DSIGNALS -D_WITH_GETLINE -DPLRS -DSAFE -DLRSLONG ${BITS} mplrs.cc lrslib.cc lrslong.cc lrsdriver.cc -o mplrs2

mplrsmp: mplrs.cc mplrs.h lrslib.cc lrslib.h lrsmp.cc lrsmp.h lrsdriver.h lrsdriver.cc
	$(mpicxx) ${CFLAGS} -DMP -DTIMES -DSIGNALS -D_WITH_GETLINE -DPLRS mplrs.cc lrslib.cc lrsmp.cc lrsdriver.cc -o mplrsmp

singlemplrs: mplrsgmp mplrs1 mplrs2

flint:	 	lrs.cc lrslib.cc lrslib.h lrsgmp.cc lrsgmp.h
		@test -d  ${INCLUDEDIR}/flint || { echo ${INCLUDEDIR}/flint not found; exit 1; }
		$(CC) -O3 -DFLINT -I/usr/local/include/flint lrs.cc lrslib.cc lrsgmp.cc lrsdriver.cc -L/usr/local/lib -Wl,-rpath=/usr/local/lib -lflint -o lrsflint -lgmp
#		$(CC) -O3 -DFLINT -I${INCLUDEDIR} -I${INCLUDEDIR}/flint lrs.cc lrsdriver.cc lrslib.cc lrsgmp.cc -L${LIBDIR} -lflint -o lrsflint -lgmp

mplrsflint:	mplrs.cc mplrs.h lrslib.cc lrslib.h lrsgmp.cc lrsgmp.h lrsdriver.cc lrsdriver.h
	${mpicxx} ${CFLAGS} -DTIMES -DSIGNALS -D_WITH_GETLINE -DFLINT -I${INCLUDEDIR}/flint -DPLRS -o mplrsflint mplrs.cc lrsdriver.cc lrslib.cc lrsgmp.cc -L${LIBDIR} -lflint -lgmp

#comment out lines with ${BITS} if __int128 not supported by your C compiler

lrsgmp:		lrs.cc lrslib.cc lrslib.h lrsgmp.cc lrsgmp.h lrsdriver.h lrsdriver.cc 
		$(CC)  ${CFLAGS}  -DGMP -I${INCLUDEDIR} -o lrsgmp lrs.cc lrslib.cc lrsgmp.cc lrsdriver.cc -L${LIBDIR}  -lgmp
		# ln -s -f lrsgmp redundgmp

single:		lrs.cc lrslong.cc lrslong.h lrslib.cc lrslib.h lrsgmp.cc lrsgmp.h lrsdriver.h lrsdriver.cc
		$(CC)  ${CFLAGS}  -DSAFE  -DLRSLONG -o lrs1 lrs.cc lrslib.cc lrslong.cc lrsdriver.cc
		$(CC)  ${CFLAGS} ${BITS} -DSAFE  -DLRSLONG -o lrs2 lrs.cc lrslib.cc lrslong.cc lrsdriver.cc
		$(CC)  ${CFLAGS} -DMP -o lrsmp lrs.cc lrslib.cc lrsdriver.cc lrsmp.cc

		# ln -s -f lrs1 redund1
		# ln -s -f lrs2 redund2

allmp:		lrs.cc lrslib.cc lrslib.h lrsmp.cc lrsmp.h lrsdriver.h lrsdriver.cc
		$(CC) -Wall -O3 -DMP  -o lrsmp lrs.cc lrslib.cc lrsdriver.cc lrsmp.cc
		$(CC) -Wall -O3  -DSAFE -DLRSLONG -o lrs1 lrs.cc lrslib.cc lrsdriver.cc lrslong.cc
		$(CC) -Wall -O3  -DSAFE -DLRSLONG ${BITS} -o lrs2 lrs.cc lrslib.cc lrsdriver.cc lrslong.cc
		$(CC) -O3 -DMP -DLRS_QUIET   -o lrsnash lrsnash.cc lrsnashlib.cc lrslib.cc lrsdriver.cc lrsmp.cc -static
		$(CC) -O3 -DMP -o setupnash setupnash.cc lrslib.cc lrsdriver.cc lrsmp.cc
		$(CC) -O3 -DMP -o setupnash2 setupnash2.cc lrslib.cc lrsdriver.cc lrsmp.cc
		$(CC) -O3  -o 2nash 2nash.cc

demo:	lpdemo1.cc lrslib.cc lrsdriver.cc lrslib.h lrsgmp.cc lrsgmp.h
	$(CC) -O3   -I${INCLUDEDIR} -L${LIBDIR} -o lpdemo1 lpdemo1.cc lrslib.cc lrsdriver.cc lrsgmp.cc -lgmp -DGMP
	$(CC) -O3   -I${INCLUDEDIR} -L${LIBDIR} -o lpdemo lpdemo.cc lrslib.cc lrsdriver.cc lrsgmp.cc -lgmp -DGMP
	$(CC) -O3   -I${INCLUDEDIR} -L${LIBDIR} -o lpdemo2 lpdemo2.cc lrslib.cc lrsdriver.cc lrsgmp.cc -lgmp -DGMP
	$(CC) -O3   -I${INCLUDEDIR} -L${LIBDIR} -o vedemo  vedemo.cc lrslib.cc lrsdriver.cc lrsgmp.cc -lgmp -DGMP
	$(CC) -O3   -I${INCLUDEDIR} -L${LIBDIR} -o chdemo  chdemo.cc lrslib.cc lrsdriver.cc lrsgmp.cc -lgmp -DGMP

lrsnash:	lrsnash.cc lrsnashlib.cc lrslib.cc lrsnashlib.h lrslib.h lrsgmp.cc lrsgmp.h lrsdriver.h lrsdriver.cc
		$(CC) -O3   -I${INCLUDEDIR} -L${LIBDIR} -o lrsnash lrsnash.cc lrsnashlib.cc lrslib.cc lrsgmp.cc lrsdriver.cc  -lgmp -DGMP

mp: lrsnash.cc lrsnashlib.cc lrslib.cc lrsnashlib.h lrslib.h lrslong.h lrsdriver.h lrsdriver.cc lrsmp.cc lrsmp.h
		$(CC) ${CPPFLAGS} -O3 -DMP -DLRS_QUIET   -o lrsnash lrsnash.cc lrsnashlib.cc lrslib.cc lrsdriver.cc lrsmp.cc -static

gmp: lrsnash.cc lrsnashlib.cc lrslib.cc lrsnashlib.h lrslib.h lrsgmp.cc lrsgmp.h lrsdriver.h lrsdriver.cc
		$(CC) ${CPPFLAGS} -O3 ${CPPFLAGS}  -I${INCLUDEDIR} -L${LIBDIR} -o lrsnashgmp lrsnash.cc lrsnashlib.cc lrslib.cc lrsgmp.cc lrsdriver.cc  -lgmp -DGMP

1: lrsnash.cc lrsnashlib.cc lrslib.cc lrsnashlib.h lrslib.h lrslong.h lrsdriver.h lrsdriver.cc
		$(CC) ${CPPFLAGS} -O3 ${CPPFLAGS}  -I${INCLUDEDIR} -L${LIBDIR} -o lrsnash1 lrsnash.cc lrsnashlib.cc lrslib.cc lrslong.cc lrsdriver.cc -DLRSLONG -DSAFE

2: lrsnash.cc lrsnashlib.cc lrslib.cc lrsnashlib.h lrslib.h lrslong.h lrsdriver.h lrsdriver.cc
		$(CC) ${CPPFLAGS} -O3 ${CPPFLAGS}  -I${INCLUDEDIR} -L${LIBDIR} -o lrsnash2 lrsnash.cc lrsnashlib.cc lrslib.cc lrslong.cc lrsdriver.cc -DLRSLONG -DSAFE ${BITS}

simple: lrsnash.cc lrsnashlib.cc lrslib.cc lrsnashlib.h lrslib.h lrsgmp.cc lrsgmp.h lrslong.h lrsdriver.h lrsdriver.cc lrsmp.cc lrsmp.h
		$(CC) ${CPPFLAGS} -O3 -DMP -DLRS_QUIET   -o lrsnash lrsnash.cc lrsnashlib.cc lrslib.cc lrsdriver.cc lrsmp.cc -static
		$(CC) ${CPPFLAGS} -O3 ${CPPFLAGS}  -I${INCLUDEDIR} -L${LIBDIR} -o lrsnashgmp lrsnash.cc lrsnashlib.cc lrslib.cc lrsgmp.cc lrsdriver.cc  -lgmp -DGMP
		$(CC) ${CPPFLAGS} -O3 ${CPPFLAGS}  -I${INCLUDEDIR} -L${LIBDIR} -o lrsnash1 lrsnash.cc lrsnashlib.cc lrslib.cc lrslong.cc lrsdriver.cc -DLRSLONG -DSAFE
		$(CC) ${CPPFLAGS} -O3 ${CPPFLAGS}  -I${INCLUDEDIR} -L${LIBDIR} -o lrsnash2 lrsnash.cc lrsnashlib.cc lrslib.cc lrslong.cc lrsdriver.cc -DLRSLONG -DSAFE ${BITS}

######################################################################
# From here on the author is David Bremner <bremner@unb.cca> to whom you should turn for help             
#
# Shared library variables
SONAME ?=liblrs.so.1
SOMINOR ?=.0.0
SHLIB ?=$(SONAME)$(SOMINOR)
SHLINK ?=liblrs.so

SHLIBOBJ2=lrslib2-shr.o lrslong2-shr.o

# for 32 bit machines

# SHLIBOBJ2=

SHLIBOBJ=lrslong1-shr.o lrslib1-shr.o  \
	lrslibgmp-shr.o lrsgmp-shr.o lrsdriver-shr.o \
	${SHLIBOBJ2}

SHLIBBIN=lrs-shared lrsnash-shared

# Building (linking) the shared library, and relevant symlinks.

${SHLIB}: ${SHLIBOBJ}
	$(CC) -shared -Wl,-soname=$(SONAME) $(SHLIBFLAGS) -o $@ ${SHLIBOBJ} -lgmp

${SONAME}: ${SHLIB}
	ln -sf ${SHLIB} ${SONAME}

${SHLINK}: ${SONAME}
	ln -sf $< $@

# binaries linked against the shared library

all-shared: ${SHLIBBIN}

lrs-shared: ${SHLINK} lrs-shared.o
	$(CC) $^ -o $@ -L . -llrs


lrsnash-shared: ${SHLINK}  lrsnash.cc
	$(CC) ${CFLAGS} -DGMP -DMA lrsnash.cc  lrsnashlib.cc -I${INCLUDEDIR} -o $@ -L . -llrs -lgmp

# driver object files

lrs-shared.o: lrs.cc
	$(CC) ${CFLAGS} -DMA ${BITS} -L${LIBDIR} -c -o $@ lrs.cc

# build object files for the shared library

lrslib1-shr.o: lrslib.cc lrslib.h
	$(CC) ${CFLAGS} ${SHLIB_CFLAGS} -DMA -DSAFE -DLRSLONG -c -o $@ lrslib.cc

lrsdriver-shr.o: lrsdriver.cc
	$(CC) ${CFLAGS} ${SHLIB_CFLAGS} -c -o $@ $<

lrslong1-shr.o: lrslong.cc lrslong.h
	$(CC) ${CFLAGS} ${SHLIB_CFLAGS} -DMA -DSAFE -DLRSLONG -c -o $@ lrslong.cc

lrslong2-shr.o: lrslong.cc lrslong.h
	$(CC) ${CFLAGS} ${SHLIB_CFLAGS} -DMA -DSAFE ${BITS} -DLRSLONG -c -o $@ lrslong.cc

lrslibgmp-shr.o: lrslib.cc lrslib.h
	$(CC) ${CFLAGS} ${SHLIB_CFLAGS} -DMA -DGMP -I${INCLUDEDIR} -c -o $@ lrslib.cc

lrsgmp-shr.o: lrsgmp.cc lrsgmp.h
	$(CC) ${CFLAGS} ${SHLIB_CFLAGS} -DMA -DGMP -I${INCLUDEDIR} -c -o $@ lrsgmp.cc

lrslib2-shr.o: lrslib.cc lrslib.h
	$(CC) ${CFLAGS} ${SHLIB_CFLAGS} -DMA -DSAFE ${BITS} -DLRSLONG -c -o $@ lrslib.cc

######################################################################
# install targets
# where to install binaries, libraries, include files
prefix ?= /usr/local
INSTALL_INCLUDES=lrslib.h lrsdriver.h lrsgmp.h lrslong.h lrsmp.h lrsrestart.h

install: all-shared install-common
	mkdir -p $(DESTDIR)${prefix}/bin
	for file in ${SHLIBBIN}; do cp $${file} $(DESTDIR)${prefix}/bin/$$(basename $$file -shared); done
	mkdir -p $(DESTDIR)${prefix}/lib
	install -t $(DESTDIR)${prefix}/lib $(SHLIB)
	cd $(DESTDIR)${prefix}/lib && ln -sf $(SHLIB) $(SHLINK)
	cd $(DESTDIR)${prefix}/lib && ln -sf $(SHLIB) $(SONAME)

install-common:
	mkdir -p $(DESTDIR)${prefix}/include/lrslib
	install -t $(DESTDIR)${prefix}/include/lrslib ${INSTALL_INCLUDES}

######################################################################
clean:		
	rm -f  lrs lrs1 lrsgmp lrs1n lpdemo lpdemo1 lpdemo2 mplrs1 mplrs mplrsmp  mplrsgmp lrs2 mplrs2 lrsflint mplrsflint *.o *.exe *.so
	rm -f  lrsmp lrsMP hvref setupnash setupnash2 lrsnashgmp lrsnash lrsnash1 lrsnash2 nashdemo 2nash vedemo checkpred inedel
	rm -f ${LRSOBJ} ${LRSOBJ64} ${SHLIBOBJ} ${SHLIB} ${SONAME} ${SHLINK}
	rm -f ${SHLIBBIN}
