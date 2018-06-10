#
# makefile of /home/treiber/trafficSim/sources/master
#
# benoetigt auch "general.cc" aus BINDIR

CFLAGS=-Wall -O
CC=g++

OBJDIR=.
LIBDIR=.
#LIBDIR=~/versionedProjects/lib/trunk
BINDIR=~/bin

# LIBDIR=.

OBJECTS=master.cc ${LIBDIR}/InOut.o Mac3phases.o FPE.o FPE_innov.o VMM.o SGM.o 

# beliebige .o Files sollen aus den jeweiligen .cc Files erzeugt werden:
.cc.o:
	${CC}  -I ${LIBDIR} ${CFLAGS} -c $<


master: $(OBJECTS)
	${CC} ${CFLAGS} -o ${BINDIR}/master -O3 $(OBJECTS) -lm 
#	rm master.o; ${CC} ${CFLAGS} -o ${BINDIR}/master -O3 master.cc -lm 



#
# Misc.
#
clean:
	rm *~ *.aux *.dvi *.log

