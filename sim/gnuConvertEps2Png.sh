#!/bin/bash


# change plotting .eps files => plotting .gnu files 
# in gnuplot files and in calling files 
# notice: size options "large", "giant" etc DOS; also set size need to be
# not larger than 1,1 since, in contrast to eps, otherwise bug

perl -i -p -e 's/set term post eps.*$/set term png enhanced font Helvetica 14 xffffff/g' `find . -name "*.gnu"`  
perl -i -p -e 's/\.eps/\.png/g' `find . -name "*.gnu"`
perl -i -p -e 's/\.eps/\.png/g' `find . -name "*.plot"`
perl -i -p -e 's/\.eps/\.png/g' `find . -name "*.run"`

# reduces sizes larger than 1,x 
# sizes smaller than 1 can remain, e.g., set size 1,0.4

perl -i -p -e 's/set size [2-9].*$/set size 1,1/g' `find . -name "*.gnu"`
perl -i -p -e 's/set size 1\.[1-9].*$/set size 1,1/g' `find . -name "*.gnu"`
perl -i -p -e 's/set size 1\,1\.2/set size 0\.85\,1/g' `find . -name "*.gnu"`

# implement robust ls 1 (each terminal has other confusing color choice)
perl -i -p -e 's/style line 1 .*$/style line 1 lt 1 lw 2 pt 7 ps 1\.9  lc rgb \"#000000\"/g' `find . -name "*.gnu"`

# change gv eps viewer into xv png viewer

perl -i -p -e 's/--orientation=PORTRAIT //g' `find . -name "*.plot"`
perl -i -p -e 's/--orientation=PORTRAIT //g' `find . -name "*.run"`
perl -i -p -e 's/gv /xv /g' `find . -name "*.plot"`
perl -i -p -e 's/gv /xv /g' `find . -name "*.run"`

# general change file settings

chmod ugo-x `find . -type f`
chmod u+x `find . -name "*.run"`
chmod u+x `find . -name "*.plot"`
chmod u+x `find . -name "*.sh"`

# remove old latex options and transform to pdflatex

perl -i -p -e 's/.*usepackage\{german\}.*\n//g' `find . -name "*.tex"`
perl -i -p -e 's/.*usepackage\{umlaute\}.*\n//g' `find . -name "*.tex"`
perl -i -p -e 's/\.eps/\.png/g'  `find . -name "*.tex"`















