
#set term png enhanced font Helvetica 14 xffffff
set term png enhanced font Helvetica 14 xffffff


##################################################
set out "GKTring-class1b.rho3d.png"

tmin=0    #h
tmax=3.3  #h

set view 30, 20
set xlabel "Location (km)"
set ylabel "Time (h)"
set title
set zlabel "{/Symbol r} [veh/km]" offset 5,3

set auto x
set yrange [tmin:]
#set xrange [-10:6]
splot "GKTring-class1b.dat" u ($1-15):($2/60):3 w l lt 2


##################################################
set out "GKTring-class1b.Qe.png"

set auto x
set xlabel  "Density (veh/km/lane)"
set auto y
set ylabel "Equilibrium flow (veh/h/lane)"
plot\
  "GKTring-class1b.tab"   u 1:4 t "Equilibrium flow"  w l lt 1 lw 3

quit

#********************************************************
# Sections x=const ("loop data")
#********************************************************

set style line 1 lt 1 lw 2 pt 7 ps 1.9  lc rgb "#000000"
set style line 2 lt 5 lw 8
set style line 3 lt 2 lw 8
set style line 4 lt 6 lw 8
set style line 5 lt 1 lw 8
set style line 6 lt 4 lw 8
set style line 10 lt 0 lw 4
set style line 11 lt 7 lw 8

set out "GKTring-class1b.tseriesUpDown.png"

set auto x
set xrange [-3:]
set xlabel "Time (min)"
set yrange [0:160]
set ylabel "V (km/h)"
set key

plot[t=0:1]\
  "GKTring-class1b.x12000"  u ($2-60):4 t "10 km upstream  " w l ls 1,\
  "GKTring-class1b.x21000" u ($2-60):4 t "1 km upstream   " w l ls 2,\
  "GKTring-class1b.x23000" u ($2-60):4 t "1 km downstream " w l ls 3,\
  "GKTring-class1b.x25000" u ($2-60):4 t "3 km downstream " w l ls 3,\
  0, 110*t t "begin/end of obstruction" w l ls 10,\
  30, 110*t t "" w l ls 10

set out "GKTring-class1b.tseriesUp.png"
plot[t=0:1]\
  "GKTring-class1b.x21000" u ($2-60):4 t " 1 km upstream " w l ls 1,\
  "GKTring-class1b.x19000"  u ($2-60):4 t " 3 km upstream " w l ls 2,\
  "GKTring-class1b.x17000"  u ($2-60):4 t " 5 km upstream " w l ls 3,\
  "GKTring-class1b.x12000"  u ($2-60):4 t "10 km upstream " w l ls 5,\
  0, 110*t t "begin/end of obstruction" w l ls 10,\
  30, 110*t t "" w l ls 10

set auto y
set out "GKTring-class1b.aUp.png"

plot\
  "GKTring-class1b.x12000" u ($2-60):6 t "10 km upstream  " w l ls 1,\
  "GKTring-class1b.x21000" u ($2-60):6 t "1 km upstream   " w l ls 2,\
  "GKTring-class1b.x23000" u ($2-60):6 t "1 km downstream " w l ls 3,\
  "GKTring-class1b.x25000" u ($2-60):6 t "3 km downstream " w l ls 4

set auto y
set out "tmp.png"

plot\
  "GKTring-class1b.x12000"  u ($2-60):5 t "10 km upstream  " w l ls 1,\
  "GKTring-class1b.x21000" u ($2-60):5 t "1 km upstream   " w l ls 2,\
  "GKTring-class1b.x23000" u ($2-60):5 t "1 km downstream " w l ls 3,\
  "GKTring-class1b.x25000" u ($2-60):5 t "3 km downstream " w l ls 4

#********************************************************
# Dynamic trajectories in Q-rho space with 
# equilibrium fundamental diagram (along x=const and t=const)
#********************************************************

set out "GKTring-class1b.Qrho.png"

set style line 1 lt 1 lw 2 pt 7 ps 1.9  lc rgb "#000000"
set style line 2 lt 5 lw 3
set style line 3 lt 2 lw 3
set style line 4 lt 6 lw 3
set style line 5 lt 1 lw 3
set style line 6 lt 4 lw 3

set xlabel "Density (veh/km/h)"
set ylabel "Flow (veh/h/lane)"
set key

tmin=40  #min
rhostart=19
Qstart=1800
rho(rhoin,t)=(t>tmin)?rhoin:rhostart
Q(Qin,t)=(t>tmin)?Qin:Qstart

plot\
  "GKTring-class1b.x19000"  u (rho($3,$2)):(Q($5,$2)) t "3 km upstream " w l ls 5,\
  "GKTring-class1b.tab"   u 1:4 t "Equilibrium flow"  w l ls 11








