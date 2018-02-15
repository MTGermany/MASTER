
set style line 1 lt 1 lw 2 pt 7 ps 1.9  lc rgb "#000000"
set style line 11 lt 7 lw 1 #schwarz (screen: tausche mit ls4=braun auf screen); 
set style line 2 lt 1 lw 5 #rot
set style line 12 lt 1 lw 1 #rot
set style line 3 lt 8 lw 5 #blassrot
set style line 13 lt 8 lw 1 #blassrot
set style line 4 lt 6 lw 5 #gelb
set style line 14 lt 6 lw 1 #gelb
set style line 5 lt 2 lw 5 #gruen
set style line 15 lt 2 lw 1 #gruen
set style line 6 lt 5 lw 5 #blasstuerkisblau
set style line 16 lt 5 lw 1 #blasstuerkisblau
set style line 7 lt 3 lw 5 #blau
set style line 17 lt 3 lw 1 #blau
set style line 8 lt 4 lw 5 #lila
set style line 18 lt 4 lw 1 #lila


set term png enhanced font Helvetica 14 xffffff

#####################################################
set palette defined ( 0 "#6600aa", 0.4  "blue", 0.6 "green",\
0.65 "yellow", 0.7 "orange", 0.8 "red", 1 "#cc00ff")

set auto z
set pm3d                     # color coded 3d surface

unset surface      # pm3d off: kein Gitternetz gemacht
                            # pm3d on: Schnelles Plotten mit/ohne Netz
                            # (je nach pm3d Einstellung) mit Artefakten
set nogrid    # sonst Bug: durchscheinendes Koordinatengitter
set isosamples 100,100
set view 45,30
set size 1,1
##################################################

set out "GKTring-class2b.rho3d.png"
print "plotting GKTring-class2b.rho3d.png"

tmin=0    #min
tmax=20  #min
xmin=-5  #km
xmax=5   #km

set xlabel "Location (km)"
set ylabel "Time (h)"
set title
set zlabel "Density (veh/km/lane)" offset 0,2

set yrange [tmin:]
set xrange [xmin:xmax]
splot "GKTring-class2b.dat" u ($1-6):($2):3 w l lt 2

set out "GKTring-class2b.v3d.png"
print "plotting GKTring-class2b.v3d.png"
set zlabel "V(km/h)" offset 0,2
splot "GKTring-class2b.dat" u ($1-6):($2):4 w l lt 2

##################################################

set size 1,1
set grid
set out "GKTring-class2b.Qe.png"
print "plotting GKTring-class2b.Qe.png"

set auto x
set xlabel  "Density (veh/km/lane)"
set auto y
set ylabel "Equilibrium flow (veh/h/lane)"
plot\
  "GKTring-class2b.tab"   u 1:4 t "Equilibrium flow"  w l ls 2

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

set out "GKTring-class2b.tseriesUpDown.png"

set auto x
set xrange [-3:]
set xlabel "Time (min)"
set yrange [0:160]
set ylabel "V (km/h)"
set key

plot[t=0:1]\
  "GKTring-class2b.x12000"  u ($2-60):4 t "10 km upstream  " w l ls 1,\
  "GKTring-class2b.x21000" u ($2-60):4 t "1 km upstream   " w l ls 2,\
  "GKTring-class2b.x23000" u ($2-60):4 t "1 km downstream " w l ls 3,\
  "GKTring-class2b.x25000" u ($2-60):4 t "3 km downstream " w l ls 3,\
  0, 110*t t "begin/end of obstruction" w l ls 10,\
  30, 110*t t "" w l ls 10

set out "GKTring-class2b.tseriesUp.png"
plot[t=0:1]\
  "GKTring-class2b.x21000" u ($2-60):4 t " 1 km upstream " w l ls 1,\
  "GKTring-class2b.x19000"  u ($2-60):4 t " 3 km upstream " w l ls 2,\
  "GKTring-class2b.x17000"  u ($2-60):4 t " 5 km upstream " w l ls 3,\
  "GKTring-class2b.x12000"  u ($2-60):4 t "10 km upstream " w l ls 5,\
  0, 110*t t "begin/end of obstruction" w l ls 10,\
  30, 110*t t "" w l ls 10

set auto y
set out "GKTring-class2b.aUp.png"

plot\
  "GKTring-class2b.x12000" u ($2-60):6 t "10 km upstream  " w l ls 1,\
  "GKTring-class2b.x21000" u ($2-60):6 t "1 km upstream   " w l ls 2,\
  "GKTring-class2b.x23000" u ($2-60):6 t "1 km downstream " w l ls 3,\
  "GKTring-class2b.x25000" u ($2-60):6 t "3 km downstream " w l ls 4

set auto y
set out "tmp.png"

plot\
  "GKTring-class2b.x12000"  u ($2-60):5 t "10 km upstream  " w l ls 1,\
  "GKTring-class2b.x21000" u ($2-60):5 t "1 km upstream   " w l ls 2,\
  "GKTring-class2b.x23000" u ($2-60):5 t "1 km downstream " w l ls 3,\
  "GKTring-class2b.x25000" u ($2-60):5 t "3 km downstream " w l ls 4

#********************************************************
# Dynamic trajectories in Q-rho space with 
# equilibrium fundamental diagram (along x=const and t=const)
#********************************************************

set out "GKTring-class2b.Qrho.png"

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
  "GKTring-class2b.x19000"  u (rho($3,$2)):(Q($5,$2)) t "3 km upstream " w l ls 5,\
  "GKTring-class2b.tab"   u 1:4 t "Equilibrium flow"  w l ls 11








