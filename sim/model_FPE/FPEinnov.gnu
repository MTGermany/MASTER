
set style line 1 lt 1 lw 2 pt 7 ps 1.9  lc rgb "#000000"
set style line 2 lt 1 lw 5 pt 5 ps 1.5 #rot
set style line 3 lt 8 lw 5 pt 5 ps 1.5 #blassrot
set style line 4 lt 6 lw 5 pt 5 ps 1.5 #gelb  (screen: tausche mit lt 1=orange-ocker)
set style line 5 lt 2 lw 5 pt 5 ps 1.5 #gruen
set style line 6 lt 5 lw 5 pt 5 ps 1.5 #blasstuerkisblau
set style line 7 lt 3 lw 5 pt 5 ps 1.5 #blau
set style line 8 lt 4 lw 5 pt 5 ps 1.5 #lila

##############################################

set term png enhanced font Helvetica 14 xffffff
set title "FPEinnov, A=,B="

set out "./FPEinnov.1.png"

set view 30,70

set auto z

set xlabel "t"
set ylabel "y"
set zlabel "P(x,t)"
splot "FPEinnov.dat" u 2:($1):3 w l lt 1


#############################################################

set out "./FPEinnov.3.png"
set auto

set key at screen 0.4,0.8

set xlabel "x"
set ylabel "P(x,t)"
set auto y

plot\
 "FPEinnov.t0" u ($1):3 t "t=0" w l ls 1,\
 "FPEinnov.t1" u ($1):3 t "t=1" w l ls 2,\
 "FPEinnov.t2" u ($1):3 t "t=2" w l ls 3,\
 "FPEinnov.t3" u ($1):3 t "t=3" w l ls 4,\
 "FPEinnov.t4" u ($1):3 t "t=4" w l ls 5

#############################################################

set out "./FPEinnov.2.png"
set auto

set key at screen 0.4,0.8

set auto x
set auto y
set xlabel "Time (a.u.)"
set ylabel "<x>, Var(x), Skewness(x)"

plot\
 "FPEinnov.stat" u 1:(1.0*($2)) t "<Delta x>" w l ls 1,\
 "FPEinnov.stat" u 1:3 t "Var(x)" w l ls 2,\
 "FPEinnov.stat" u 1:(10*$4) t "10*Skewness" w l ls 3


#############################################################

set out "./FPEinnov.4.png"

set key at screen 0.55,0.8

set auto

set xlabel "t"
set ylabel "P(x,t)"

plot\
 "FPEinnov.x10"  u 2:3 t "x=10" w l ls 1,\
 "FPEinnov.x20"  u 2:3 t "x=20" w l ls 2,\
 "FPEinnov.x30"  u 2:3 t "x=30" w l ls 3,\
 "FPEinnov.x40"  u 2:3 t "x=40" w l ls 4,\
 "FPEinnov.x50"  u 2:3 t "x=50" w l ls 5,\
 "FPEinnov.x60"  u 2:3 t "x=60" w l ls 6,\
 "FPEinnov.x70"  u 2:3 t "x=70" w l ls 7,\
 "FPEinnov.x80"  u 2:3 t "x=80" w l ls 8

