
####################################################
# von  ~/info/gnuColoredContour/coloredContourAnd3dCurves.gnu
#
# plottet farbige 3d- oder map-Oberflaechen und Kurven,
# z.B. Trajektorien auf V(x,t)-Plots
# Feature ermoeglicht durch "explicit" Option von set pm3d:
# "set pm3d explicit <weitere Optionen>"
#
# Neue gnuplot42-Syntax! 
# (ohne interpolate-Option geht auch gnuplot4, aber Bug bei multiplot)
#
# siehe auch
# ~/info/gnuplot,   
# ~/info/gnuTemplate.gnu,  
# ~/info/gnuColoredContour/*.gnu
####################################################


set style line 1 lt 1 lw 2 pt 7 ps 1.9  lc rgb "#000000"
set style line 11 lt 7 lw 8 pt 7 ps 0.2
set style line 2 lt 1 lw 4 pt 13 ps 0.5 #rot
set style line 3 lt 8 lw 3 pt 5 ps 2 #blassrot
set style line 4 lt 6 lw 3 pt 5 ps 2 #gelb  (screen: tausche mit ls 1=orange-ocker auf screen)
set style line 5 lt 2 lw 3 pt 5 ps 2 #gruen
set style line 6 lt 5 lw 3 pt 5 ps 2 #blasstuerkisblau
set style line 7 lt 3 lw 4 pt 7 ps 1 #blau
set style line 8 lt 4 lw 3 pt 5 ps 2 #lila
set style line 9 lt 9 lw 3 pt 5 ps 2
set style line 10 lt 10 lw 3 pt 5 ps 2
set style line 99 lt 7 lw 0 linecolor rgb "#000000" # beliebige Farben:Schwarz
set style line 98 lt 7 lw 4 linecolor rgb "#000000" # beliebige Farben:Schwarz


set term png enhanced font Helvetica 14 xffffff

######################################################
set out "LWR_platoonIC.fund.png"
print "plotting LWR_platoonIC.fund.png"
######################################################


set xlabel "Density {/Symbol r} [km^{-1}]"
set ylabel "Flow Q [h^{-1}]"
plot "LWR_platoonIC.tab" u 1:4 w l ls 2
set xlabel ""
set ylabel ""



######################################################
set out "LWR_platoonIC.rho3d.png"
print "plotting LWR_platoonIC.rho3d.png"
######################################################



rhomin=26

set nokey
set label "Location [km]" at screen 1.4,0.35  rotate by 33
set xrange [3:11.8] reverse
set xtics 2

set label "Time [minutes]" at screen 0.4,0.52  rotate by -30
#set ylabel "Time [minutes]" offset 0,2
set yrange [0:20]

set zlabel "{/Symbol r} [km^{-1}]" offset 7,7
set label  "{/Symbol r} [km^{-1}]" at screen 1.77,1.65
set zrange [rhomin:42]
set ztics 10

set ticslevel  0

set view 14, 134
set size 1,1
set nogrid

#set size 1,1


set pm3d 
#set pm3d; set pm3d interpolate 4,4
#set pm3d; set pm3d interpolate 2,2 hidden3d 99
unset surface
set hidden3d


# bug: contour geht nicht gleichzeitig mit color, entweder contour oder color!!
set contour surface      # Aktiviert Kontourlinien auf 3D-Flaeche
#set cntrparam levels discrete 10,20,40 # freely set lines
set cntrparam cubicspline
set cntrparam levels incremental 27.2, 0.8, 40.8
unset clabel  # dann lauter gleiche Kontourlinien; 
                     # Farbe und Typ mit "w l ls" beim splot-Kommando
                     # Achtung bug! explizit ("w l lt 3 lw 1" etc) gibt DOS!

#set palette defined ( 0 "red", 20 "orange", 40 "yellow",\
#      60 "green", 80 "blue",  100 "#dd00ff" ) 
set palette defined ( 0 "#2200aa", 0.2  "blue", 0.4 "green",0.5 "yellow",\
 0.7 "orange", 0.9 "red", 1 "#cc00ff")




splot "LWR_platoonIC.dat" u 1:2:3 w l ls 99

quit






##################################################
set out "LWR_platoonIC.v3d.png"
#set title "V (km/h)"
#set zlabel "V (km/h)"

print "plotting LWR_platoonIC.v3d.png ..."

splot "LWR_platoonIC.dat" u 1:($2):4 w l ls 99

##################################################
#set out "LWR_platoonIC.Q3d.png"
#set zrange [1500:]
#set zlabel "Flow (veh/h/lane)" offset 0,2
#splot "LWR_platoonIC.dat" u ($1-15):2:5 w l lt 2

##################################################
set out "LWR_platoonIC.Qe.png"
print "plotting LWR_platoonIC.Qe.png ..."
set auto x
set xlabel  "Density (veh/km/lane)"
set auto y
set ylabel "Equilibrium flow (veh/h/lane)"
plot\
  "LWR_platoonIC.tab"   u 1:4 t "Equilibrium flow"  w l lt 1 lw 3

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

set out "LWR_platoonIC.tseriesUpDown.png"

set auto x
set xrange [-3:]
set xlabel "Time (min)"
set yrange [0:160]
set ylabel "V (km/h)"
set key

plot[t=0:1]\
  "LWR_platoonIC.x12000"  u ($2-60):4 t "10 km upstream  " w l ls 1,\
  "LWR_platoonIC.x21000" u ($2-60):4 t "1 km upstream   " w l ls 2,\
  "LWR_platoonIC.x23000" u ($2-60):4 t "1 km downstream " w l ls 3,\
  "LWR_platoonIC.x25000" u ($2-60):4 t "3 km downstream " w l ls 3,\
  0, 110*t t "begin/end of obstruction" w l ls 10,\
  30, 110*t t "" w l ls 10

set out "LWR_platoonIC.tseriesUp.png"
plot[t=0:1]\
  "LWR_platoonIC.x21000" u ($2-60):4 t " 1 km upstream " w l ls 1,\
  "LWR_platoonIC.x19000"  u ($2-60):4 t " 3 km upstream " w l ls 2,\
  "LWR_platoonIC.x17000"  u ($2-60):4 t " 5 km upstream " w l ls 3,\
  "LWR_platoonIC.x12000"  u ($2-60):4 t "10 km upstream " w l ls 5,\
  0, 110*t t "begin/end of obstruction" w l ls 10,\
  30, 110*t t "" w l ls 10

set auto y
set out "LWR_platoonIC.aUp.png"

plot\
  "LWR_platoonIC.x12000" u ($2-60):6 t "10 km upstream  " w l ls 1,\
  "LWR_platoonIC.x21000" u ($2-60):6 t "1 km upstream   " w l ls 2,\
  "LWR_platoonIC.x23000" u ($2-60):6 t "1 km downstream " w l ls 3,\
  "LWR_platoonIC.x25000" u ($2-60):6 t "3 km downstream " w l ls 4

set auto y
set out "tmp.png"

plot\
  "LWR_platoonIC.x12000"  u ($2-60):5 t "10 km upstream  " w l ls 1,\
  "LWR_platoonIC.x21000" u ($2-60):5 t "1 km upstream   " w l ls 2,\
  "LWR_platoonIC.x23000" u ($2-60):5 t "1 km downstream " w l ls 3,\
  "LWR_platoonIC.x25000" u ($2-60):5 t "3 km downstream " w l ls 4

#********************************************************
# Dynamic trajectories in Q-rho space with 
# equilibrium fundamental diagram (along x=const and t=const)
#********************************************************

set out "LWR_platoonIC.Qrho.png"

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
  "LWR_platoonIC.x19000"  u (rho($3,$2)):(Q($5,$2)) t "3 km upstream " w l ls 5,\
  "LWR_platoonIC.tab"   u 1:4 t "Equilibrium flow"  w l ls 11








