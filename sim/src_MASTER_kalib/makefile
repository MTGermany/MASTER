#
# makefile of general c++ projects
# Achtung! Auch ohne .h File muss man bei $OBJECTS immer auch das 
# File mit der Main-Methode dazunehmen!
# 

#compiler
CC=g++ -Wall -O3

#output:
OUTNAME_TEMPL=MASTER_kalib

#include path for headers
# !!! ACHTUNG keine Comments hinter Befehle/Defs => obskurer Fehler!!

LIBDIR=~/versionedProjects/lib/trunk
BINDIR=~/bin
#STATDIR=~/trafficSim/sources/statistics  # bis jan09
#STATDIR=~/trafficSim/sources/lib


#objects/other programs for linking

OBJECTS_TEMPLMAIN=MASTER_kalib.o
LIBOBJECTS=general.o InOut.o RandomUtils.o Statistics.o


###########################################################
# suffix regel: mache aus *.cpp ein *.o, und zwar fuer eingabedatei $<
.cpp.o:
	${CC} -I ${LIBDIR} -I ${STATDIR} -c $<

############################################################

#wenn man nur h-file veraendert, sollte make alles neu compilieren
#irgendwie muss das gehen!


#o-files werden in "." als target verzeichnis geschrieben!
#dann wird alles von hier gelinkt. sehr gut so 
#(weil InOut.o ja zu diesem build gehoert)
general.o: ${LIBDIR}/general.cpp
	${CC} -c ${LIBDIR}/general.cpp -o general.o

InOut.o: ${LIBDIR}/InOut.cpp
	${CC} -c ${LIBDIR}/InOut.cpp -o InOut.o

Statistics.o: ${LIBDIR}/Statistics.cpp
	${CC} -c ${LIBDIR}/Statistics.cpp -o Statistics.o

Math.o: ${LIBDIR}/Math.cpp
	${CC} -c ${LIBDIR}/Math.cpp -o Math.o

RandomUtils.o: ${LIBDIR}/RandomUtils.cpp
	${CC} -c ${LIBDIR}/RandomUtils.cpp -o RandomUtils.o

##########################################################


#target micNew: abhaengigkeiten in Macro $(OBJECTS_TEMPLMAIN) zusammengefasst 
#wichtig ist dass naechste befehlszeile mit TAB anfaengt:
MASTER_kalib: $(OBJECTS_TEMPLMAIN)  $(LIBOBJECTS)
	${CC} -o ${BINDIR}/${OUTNAME_TEMPL} -O3 $(OBJECTS_TEMPLMAIN) $(LIBOBJECTS) -lm

makedistr:  $(OBJECTS_DIST) 
	${CC} -o ${BINDIR}/makedistr -O3 $(OBJECTS_DIST) -lm 


#
# Misc.
#
clean:
	rm -v *.o


#-I directory
#Mit der -I-Compileroption kann man beim Übersetzen eines Programmes die Liste der Verzeichnisse erweitern, in denen nach einer Datei gesucht wird.
# gcc -Iinc hello.c
#sucht nach stdio.h zuerst als inc/stdio.h, und erst dann als /usr/include/stdio.h.


#The -I option or the INCLUDE_PATH variable described below should always be used to set the include path for header files.