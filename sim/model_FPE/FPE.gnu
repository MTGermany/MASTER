
set style line 1 lt 1 lw 2 pt 7 ps 1.9  lc rgb "#000000"
set style line 2 lt 1 lw 5 pt 5 ps 1.5 #rot
set style line 3 lt 8 lw 5 pt 5 ps 1.5 #blassrot
set style line 4 lt 6 lw 5 pt 5 ps 1.5 #gelb  (screen: tausche mit lt 1=orange-ocker)
set style line 5 lt 2 lw 5 pt 5 ps 1.5 #gruen
set style line 6 lt 5 lw 5 pt 5 ps 1.5 #blasstuerkisblau
set style line 7 lt 3 lw 5 pt 5 ps 1.5 #blau
set style line 8 lt 4 lw 5 pt 5 ps 1.5 #lila

##############################################

set title "FPE, A=,B="


set term png enhanced font Helvetica 14 xffffff

set out "./FPE.1.png"

set view 50,100


set xlabel "t"
set ylabel "y"
set zlabel "P(x,t)"
splot "FPE.dat" u 2:($1):3 w l lt 1

#############################################################

set term png enhanced font Helvetica 14 xffffff
set out "./FPE.2.png"

set key at screen 0.55,0.8


set xlabel "t"
set ylabel "P(x,t)"

plot\
 "FPE.x10"  u 2:3 t "x=10" w l ls 1,\
 "FPE.x20"  u 2:3 t "x=20" w l ls 2,\
 "FPE.x30"  u 2:3 t "x=30" w l ls 3,\
 "FPE.x40"  u 2:3 t "x=40" w l ls 4,\
 "FPE.x50"  u 2:3 t "x=50" w l ls 5,\
 "FPE.x60"  u 2:3 t "x=60" w l ls 6,\
 "FPE.x70"  u 2:3 t "x=70" w l ls 7,\
 "FPE.x80"  u 2:3 t "x=80" w l ls 8

#############################################################

set out "./FPE.3.png"
set auto

set key at screen 0.4,0.8


set xlabel "x"
set ylabel "P(x,t)"


plot\
 "FPE.t0" u ($1):3 t "t=0" w l ls 1,\
 "FPE.t1" u ($1):3 t "t=1" w l ls 2,\
 "FPE.t2" u ($1):3 t "t=2" w l ls 3,\
 "FPE.t3" u ($1):3 t "t=3" w l ls 4,\
 "FPE.t4" u ($1):3 t "t=4" w l ls 5,\
 "FPE.t5" u ($1):3 t "t=5" w l ls 6

