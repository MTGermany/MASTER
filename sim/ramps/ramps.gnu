
set term png enhanced font Helvetica 14 xffffff

tmin=20    #minutes
tmax=120  #minutes
xmin=5   # km
xmax=18 #km
set view 12, 60
set size 1,1

##################################################
set out "ramps.rho3d.png"
print "plotting ramps.rho3d.png ..."



set xlabel "x (km)"
set xrange [xmin:xmax]

set ylabel "t (min)"
set yrange [tmin:]
set ytics 20

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


splot "ramps.dat" u 1:($2):3 w l lt 2

##################################################
set out "ramps.v3d.png"
#set title "V (km/h)"
#set zlabel "V (km/h)"

print "plotting ramps.v3d.png ..."

splot "ramps.dat" u 1:($2):4 w l lt 2

##################################################
#set out "ramps.Q3d.png"
#set zrange [1500:]
#set zlabel "Flow (veh/h/lane)" offset 0,2
#splot "ramps.dat" u ($1-15):2:5 w l lt 2

##################################################
set out "ramps.Qe.png"
print "plotting ramps.Qe.png ..."
set auto x
set xlabel  "Density (veh/km/lane)"
set auto y
set ylabel "Equilibrium flow (veh/h/lane)"
plot\
  "ramps.tab"   u 1:4 t "Equilibrium flow"  w l lt 1 lw 3



#********************************************************
# Sections x=const ("loop data")
#********************************************************

set style line 1 lt 1 lw 2 pt 7 ps 1.9  lc rgb "#000000"
set style line 2 lt 5 lw 3
set style line 3 lt 2 lw 3
set style line 4 lt 6 lw 3
set style line 5 lt 1 lw 3
set style line 6 lt 4 lw 3
set style line 10 lt 0 lw 2
set style line 11 lt 7 lw 3

set out "ramps.tseries_V.png"
print "plotting ramps.tseries_V.png"

set auto x
set xrange [-3:]
set xlabel "Time (min)"
set yrange [0:160]
set ylabel "V (km/h)"
set key

plot[t=0:1]\
  "ramps.x12000"  u ($2-60):4 t "2 km upstream  " w l ls 1,\
  "ramps.x14000" u ($2-60):4 t "at oramp   " w l ls 2,\
  "ramps.x16000" u ($2-60):4 t "2 km downstream " w l ls 3



#********************************************************
# Dynamic trajectories in Q-rho space with 
# equilibrium fundamental diagram (along x=const and t=const)
#********************************************************

set out "ramps.Qrho.png"
print "plotting ramps.Qrho.png"
set style line 1 lt 1 lw 2 pt 7 ps 1.9  lc rgb "#000000"
set style line 2 lt 5 lw 1
set style line 3 lt 2 lw 1
set style line 4 lt 6 lw 1
set style line 5 lt 1 lw 1
set style line 6 lt 4 lw 1
set style line 11 lt 1 lw 2 lc rgb "#000000"

set xlabel "Density (veh/km/h)"
set ylabel "Flow (veh/h/lane)"
set ytics 500
set key

tmin=40  #min
rhostart=19
Qstart=1800
rho(rhoin,t)=(t>tmin)?rhoin:rhostart
Q(Qin,t)=(t>tmin)?Qin:Qstart
set auto

plot\
  "ramps.x12000"  u (rho($3,$2)):(Q($5,$2)) t "2 km upstream " w lp ls 5,\
  "ramps.tab"   u 1:4 t "Steady-state flow"  w l ls 11








