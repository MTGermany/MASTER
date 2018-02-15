
set style line 1 lt 1 lw 2 pt 7 ps 1.9  lc rgb "#000000"
set style line 2 lt 1 lw 2 pt 2 ps 1.5  lc rgb "#CC0022" #rot, solid, Kreuz
set style line 3 lt 8 lw 2 pt 4 ps 1.2  lc rgb "#FF3300"#orange, offenes Quadrat
set style line 4 lt 6 lw 2 pt 4 ps 1.5  lc rgb "#FFAA00"  #gelb, offenes Quadrat
set style line 5 lt 1 lw 2 pt 5 ps 1.5  lc rgb "#00DD22"  #gruen,solid,ClosedBox
set style line 6 lt 5 lw 2 pt 4 ps 1.5  lc rgb "#00AAAA" #offenes Quadrat
set style line 7 lt 3 lw 2 pt 4 ps 2.0  lc rgb "#1100FF"  #blau,gepunktet,offenes Quadrat
set style line 8 lt 4 lw 2 pt 8 ps 1.5  lc rgb "#220088"
set style line 9 lt 7 lw 2 pt 9 ps 1.5  lc rgb "#999999"  #grau, aufr. gschl. Dreieck

set style line 11 lt 1 lw 6 pt 7 ps 1.9  lc rgb "#000000" #schwarz,solid,bullet
set style line 12 lt 1 lw 6 pt 2 ps 1.5  lc rgb "#CC0022" #rot, dash, Kreuz
set style line 13 lt 8 lw 6 pt 4 ps 1.2  lc rgb "#FF3300"#orange, offenes Quadrat
set style line 14 lt 6 lw 6 pt 4 ps 1.5  lc rgb "#FFAA00"  #gelb, offenes Quadrat
set style line 15 lt 1 lw 6 pt 5 ps 1.5  lc rgb "#00DD22"  #gruen,solid,ClosedBox
set style line 16 lt 5 lw 6 pt 7 ps 1.5  lc rgb "#00AAAA" #offener Kreis
set style line 17 lt 1 lw 6 pt 7 ps 1.5  lc rgb "#1100FF"  #blau,solid,Bullet
set style line 18 lt 4 lw 6 pt 8 ps 1.5  lc rgb "#220088"
set style line 19 lt 7 lw 6 pt 9 ps 1.5  lc rgb "#999999"  #grau, aufr. gschl. Dreieck

set term png enhanced font Helvetica 14 xffffff

xoffset_km=15  # pos of blockage
tbegin_s=3600  # begin blockage

########################################################
print "plotting incidentConvexFD.rho3d.png"
set out "incidentConvexFD.rho3d.png"
########################################################


set view 15, 70
set xlabel "Location (km)"
set xrange [-10:6]
set auto x
set ylabel "Time (s)"
set yrange [-300:2400]
set title
set zlabel "Density (veh/km/lane)" offset 0,2

splot "incidentConvexFD.dat" u ($1-xoffset_km):(60*$2-tbegin_s):3 w l lt 2


########################################################
set out "incidentConvexFD.tseries_V.png"
print "plotting incidentConvexFD.tseries_V.png"
########################################################

set auto x
set xrange [-3:]
set xlabel "Time (min)"
set yrange [0:160]
set ylabel "V (km/h)"
set key

plot[t=0:1]\
  "incidentConvexFD.x12000"  u ($2-60):4 t "10 km upstream  " w l ls 1,\
  "incidentConvexFD.x21000" u ($2-60):4 t "1 km upstream   " w l ls 2,\
  "incidentConvexFD.x23000" u ($2-60):4 t "1 km downstream " w l ls 3,\
  "incidentConvexFD.x25000" u ($2-60):4 t "3 km downstream " w l ls 4,\
  0, 110*t t "begin/end of obstruction" w l ls 11,\
  30, 110*t t "" w l ls 11


set auto y
########################################################
set out "incidentConvexFD.tseries_acc.png"
print "plotting incidentConvexFD.tseries_acc.png"
########################################################

set ylabel "acc=(d/dt+V d/dx) V"
plot\
  "incidentConvexFD.x12000" u ($2-60):6 t "10 km upstream  " w l ls 1,\
  "incidentConvexFD.x21000" u ($2-60):6 t "1 km upstream   " w l ls 2,\
  "incidentConvexFD.x23000" u ($2-60):6 t "1 km downstream " w l ls 3,\
  "incidentConvexFD.x25000" u ($2-60):6 t "3 km downstream " w l ls 4

set auto y
set out "incidentConvexFD.tseries_Q.png"
print "plotting incidentConvexFD.tseries_Q.png"
set ylabel "Q [veh/h]"

plot\
  "incidentConvexFD.x12000"  u ($2-60):5 t "10 km upstream  " w l ls 1,\
  "incidentConvexFD.x21000" u ($2-60):5 t "1 km upstream   " w l ls 2,\
  "incidentConvexFD.x23000" u ($2-60):5 t "1 km downstream " w l ls 3,\
  "incidentConvexFD.x25000" u ($2-60):5 t "3 km downstream " w l ls 4


########################################################
set out "incidentConvexFD.snapshots_rho.png"
print "plotting incidentConvexFD.snapshots_rho.png"
########################################################

#xoffset_km=15  # pos of blockage
#tbegin_s=3600  # begin blockage

set xlabel "x [m]"
set xrange [-10000:10000]
set ylabel "Density [veh/km]"
set yrange [0:160]

set key top center box

plot\
  "incidentConvexFD.t4200"  u (1000*($1-xoffset_km)):3 t "600 s after block"  w l ls 11,\
  "incidentConvexFD.t4800"  u (1000*($1-xoffset_km)):3 t "1200 s=lift time"  w l ls 2,\
  "incidentConvexFD.t5000"  u (1000*($1-xoffset_km)):3 t "200 s  after lift"  w l ls 3,\
  "incidentConvexFD.t5200"  u (1000*($1-xoffset_km)):3 t "400 s \" "  w l ls 6

set auto x
set auto y


########################################################
set out "incidentConvexFD.Qrho.png"
print "plotting incidentConvexFD.Qrho.png"
########################################################

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
  "incidentConvexFD.x19000"  u (rho($3,$2)):(Q($5,$2)) t "3 km upstream " w l ls 5,\
  "incidentConvexFD.tab"   u 1:4 t "Equilibrium flow"  w l ls 11








