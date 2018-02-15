
set term png enhanced font Helvetica 14 xffffff

tmin=60    #minutes
tmax=180  #minutes
set view 15, 60
set size 1,1

##################################################
set out "rampsGen.rho3d.png"
print "plotting rampsGen.rho3d.png ..."


set ylabel "t (min)"
set yrange [tmin:]

set xlabel "x (km)"
set auto x

#set title "{/Symbol r}(vehicles/km/lane)" 
#set zlabel "{/Symbol r}(vehicles/km/lane)" offset 0,2
set zlabel "V (km/h)" offset 0,2
set ztics 50
set zrange [:] reverse


set cntrparam levels discrete 10,20,40 # freely set lines
unset clabel  # dann lauter gleiche Kontourlinien; 
                     # Farbe und Typ mit "w l ls" beim splot-Kommando
#set palette defined ( 0 "#2200aa", 0.2  "blue", 0.4 "green",0.5 "yellow", 0.7 "orange", 0.9 "red", 1 "#cc00ff")
set palette defined (0 "#cc00ff", 0.15 "red",  0.25 "orange", 0.35 "yellow",0.5 "green",0.7  "blue", 1 "#2200aa")

#set cbrange [0:120]  #[min:max] of color-coded  range (default zrange)

set pm3d                     # color coded 3d surface
#set pm3d     map                # color coded 3d surface
#set contour surface      # Aktiviert Kontourlinien auf 3D-Flaeche
unset surface      # pm3d off: kein Gitternetz gemacht
                            # pm3d on: Schnelles Plotten mit/ohne Netz
                            # (je nach pm3d Einstellung) mit Artefakten


splot "rampsGen.dat" u 1:($2):3 w l lt 2 

##################################################
set out "rampsGen.v3d.png"
set title "V (km/h)"
#set zlabel "V (km/h)"

print "plotting rampsGen.v3d.png ..."

splot "rampsGen.dat" u 1:($2):4 w l lt 2 

##################################################
#set out "rampsGen.Q3d.png"
#set zrange [1500:]
#set zlabel "Flow (veh/h/lane)" offset 0,2
#splot "rampsGen.dat" u ($1-15):2:5 w l lt 2 

##################################################
set out "rampsGen.Qe.png"
print "plotting rampsGen.Qe.png ..."
set auto x
set xlabel  "Density (veh/km/lane)"
set auto y
set ylabel "Equilibrium flow (veh/h/lane)"
plot\
  "rampsGen.tab"   u 1:4 t "Equilibrium flow"  w l lt 1 lw 3

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

set out "rampsGen.tseriesUpDown.png"

set auto x
set xrange [-3:]
set xlabel "Time (min)"
set yrange [0:160]
set ylabel "V (km/h)"
set key

plot[t=0:1]\
  "rampsGen.x12000"  u ($2-60):4 t "10 km upstream  " w l ls 1,\
  "rampsGen.x21000" u ($2-60):4 t "1 km upstream   " w l ls 2,\
  "rampsGen.x23000" u ($2-60):4 t "1 km downstream " w l ls 3,\
  "rampsGen.x25000" u ($2-60):4 t "3 km downstream " w l ls 3,\
  0, 110*t t "begin/end of obstruction" w l ls 10,\
  30, 110*t t "" w l ls 10

set out "rampsGen.tseriesUp.png"
plot[t=0:1]\
  "rampsGen.x21000" u ($2-60):4 t " 1 km upstream " w l ls 1,\
  "rampsGen.x19000"  u ($2-60):4 t " 3 km upstream " w l ls 2,\
  "rampsGen.x17000"  u ($2-60):4 t " 5 km upstream " w l ls 3,\
  "rampsGen.x12000"  u ($2-60):4 t "10 km upstream " w l ls 5,\
  0, 110*t t "begin/end of obstruction" w l ls 10,\
  30, 110*t t "" w l ls 10

set auto y
set out "rampsGen.aUp.png"

plot\
  "rampsGen.x12000" u ($2-60):6 t "10 km upstream  " w l ls 1,\
  "rampsGen.x21000" u ($2-60):6 t "1 km upstream   " w l ls 2,\
  "rampsGen.x23000" u ($2-60):6 t "1 km downstream " w l ls 3,\
  "rampsGen.x25000" u ($2-60):6 t "3 km downstream " w l ls 4

set auto y
set out "tmp.png"

plot\
  "rampsGen.x12000"  u ($2-60):5 t "10 km upstream  " w l ls 1,\
  "rampsGen.x21000" u ($2-60):5 t "1 km upstream   " w l ls 2,\
  "rampsGen.x23000" u ($2-60):5 t "1 km downstream " w l ls 3,\
  "rampsGen.x25000" u ($2-60):5 t "3 km downstream " w l ls 4

#********************************************************
# Dynamic trajectories in Q-rho space with 
# equilibrium fundamental diagram (along x=const and t=const)
#********************************************************

set out "rampsGen.Qrho.png"

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
  "rampsGen.x19000"  u (rho($3,$2)):(Q($5,$2)) t "3 km upstream " w l ls 5,\
  "rampsGen.tab"   u 1:4 t "Equilibrium flow"  w l ls 11








