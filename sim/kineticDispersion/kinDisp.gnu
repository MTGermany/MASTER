
set term png enhanced font Helvetica 14 xffffff

#####################################
print "plotting kinDisp_V3d.png"
set out "kinDisp_V3d.png"
#####################################

set view 15, 70
set xlabel "Location (km)"
set auto x

set ylabel "Time (min)"
set zlabel "V (km/h)" offset 0,2

set title

splot "kinDisp.dat" u ($1):2:4 w l lt 2

#####################################
print "plotting kinDisp_rho3d.png"
set out "kinDisp_rho3d.png"
#####################################

set zlabel "{/Symbol r} [veh/km]" offset 0,2
splot "kinDisp.dat" u ($1):2:3 w l lt 2




set auto

#********************************************************
# Sections x=const ("loop data")
#********************************************************

set style line 1 lt 1 lw 2 pt 7 ps 1.9  lc rgb "#000000"
set style line 2 lt 5 lw 8
set style line 3 lt 2 lw 8
set style line 4 lt 6 lw 8
set style line 5 lt 1 lw 8
set style line 6 lt 4 lw 8
set style line 10 lt 0 lw 8
set style line 11 lt 7 lw 8

set out "kinDisp_tseries.png"
print "plotting kinDisp_tseries.png"

tmin=0   #min
tmax=10  #min
set xrange [tmin:tmax]
set xlabel "Time (min)"
set yrange [0:160]
set ylabel "V (km/h)"
set key at screen 0.7,0.4
v0=120   #km/h

#plot[t=0:1]\
#  "kinDisp.x45000" u ($2):4 t "40 km downstream " w l ls 6,\
#  "kinDisp.x25000" u ($2):4 t "20 km downstream " w l ls 5,\
#  "kinDisp.x15000" u ($2):4 t "10 km downstream " w l ls 4,\

plot[t=0:1]\
  "kinDisp.x15000" u ($2):4 t "10 km downstream " w l ls 6,\
  "kinDisp.x12000" u ($2):4 t "7 km downstream " w l ls 5,\
  "kinDisp.x10000" u ($2):4 t "5 km downstream " w l ls 4,\
  "kinDisp.x8000"  u ($2):4 t "3 km downstream  " w l ls 3,\
  "kinDisp.x7000"  u ($2):4 t "2 km downstream  " w l ls 2,\
  "kinDisp.x6000"  u ($2):4 t "1 km downstream  " w l ls 1,\
tmin+t*tmax, v0 t "Avg. desired elocity" w l ls 11


quit

#********************************************************
# Dynamic trajectories in Q-rho space with 
# equilibrium fundamental diagram (along x=const and t=const)
#********************************************************

set out "kinDisp_Qrho.png"


set xlabel "Density (veh/km/h)"
set ylabel "Flow (veh/h/lane)"

set key
plot\
  "kinDisp.x8000"  u 3:5 t "Trajectory 3 km downstream"  w linesp ls 2,\
  "kinDisp.x25000" u 3:5 t "Trajectory 20 km downstream"   w linesp ls 5,\
  "kinDisp.tab"   u 1:4 t "Equilibrium flow"            w l ls 10








