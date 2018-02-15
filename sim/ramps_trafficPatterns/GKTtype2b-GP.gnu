
set style line 1 lt 1 lw 2 pt 7 ps 1.9  lc rgb "#000000"
set style line 2 lt 3 lw 6 pt 12 #blau
set style line 3 lt 4 lw 6 #lila
set style line 4 lt 1 lw 6 pt 10 #rot
set style line 5 lt 6 lw 6 #gelb  (screen: tausche mit ls 1=orange-ocker auf screen)
set style line 6 lt 2 lw 6 pt 14 ps 0.6 #gruen
set style line 7 lt 7 lw 6 #schwarz
set style line 8 lt 8 lw 6 #blassrot

set style line 11 lt 7 lw 1 #schwarz
set style line 12 lt 3 lw 1 pt 12 #blau
set style line 13 lt 4 lw 1 #lila
set style line 14 lt 1 lw 1 #lila


set term png enhanced font Helvetica 14 xffffff

tmin=30    #minutes
tmax=150  #minutes
xmin=5   # km
xmax=18 #km
#set view 12, 60
set view 0, 89.9
set size 1,1

##################################################
set out "GKTtype2b-GP.rho3d.png"
print "plotting GKTtype2b-GP.rho3d.png ..."



set ylabel "x (km)" offset 0.5, 0
set yrange [xmin:xmax]

set xlabel "t (min)" offset 0, 0.5
set xrange [tmin:tmax]
set xtics 30

#set title "{/Symbol r}(vehicles/km/lane)" 
#set cblabel "{/Symbol r}(vehicles/km/lane)" 

set noztics 
#set zrange [:]


set cntrparam levels discrete 10,20,40 # freely set lines
unset clabel  # dann lauter gleiche Kontourlinien; 
                     # Farbe und Typ mit "w l ls" beim splot-Kommando
#set palette defined (0 "#cc00ff", 0.15 "red",  0.25 "orange", 0.35 "yellow",0.5 "green",0.600 "blue", 1 "#2200aa")
set palette defined (0 "#cc00ff", 0.15 "red",  0.25 "orange", 0.35 "yellow",0.45 "green",0.6  "blue", 1 "#2200aa")

#set cbrange [0:120]  #[min:max] of color-coded  range (default zrange)

set pm3d                     # color coded 3d surface
set pm3d     map                # color coded 3d surface
#set contour surface      # Aktiviert Kontourlinien auf 3D-Flaeche
unset surface      # pm3d off: kein Gitternetz gemacht
                            # pm3d on: Schnelles Plotten mit/ohne Netz
                            # (je nach pm3d Einstellung) mit Artefakten


splot "GKTtype2b-GP.dat" u 2:1:3 w l lt 2

##################################################
set out "GKTtype2b-GP.v3d.png"
#set title "V (km/h)"
#set zlabel "V (km/h)"
set label 1 "V (km/h)" at screen 1.73,1.62
#set cbrange [0:100]

print "plotting GKTtype2b-GP.v3d.png ..."

splot "GKTtype2b-GP.dat" u 2:1:4 w l lt 2

##################################################
#set out "GKTtype2b-GP.Q3d.png"
#set zrange [1500:]
#set zlabel "Flow (veh/h/lane)" offset 0,2
#splot "GKTtype2b-GP.dat" u ($1-15):2:5 w l lt 2

##################################################

####################################
set term png enhanced font Helvetica 14 xffffff
set out "GKTtype2b-GP_fund.png"
print "plotting GKTtype2b-GP_fund.png ..."
set size 0.85,1
set xlabel "{/Symbol r}(vehicles/km)" offset 0,0.5
set xrange [0:130]
set ylabel "Q (vehicles/h)" offset 1,0
set auto y
set yrange [:] # "reverse" switched anscheinend hin- und her!
set nokey
set nogrid

rhoc1=33.
rhoc2=47.
rhocconv=47.
rhoc3=60.
rhoc4=62.

Qc1=1900.
Qc2=1760.
Qcconv=1760.
Qc3=1150.
Qc4=1100.

set ytics 500
rholine=40.
Qline=500.
set label 1 "{/Symbol r}_{1}" at rhoc1+rholine+5, 1000+Qline
set label 2 "{/Symbol r}_{2}" at rhoc2+rholine+5, 800+Qline
set label 3 "{/Symbol r}_{3}" at rhoc3+rholine+5, 600+Qline
set label 4 "{/Symbol r}_{4}" at rhoc4+rholine+5, 300+Qline

plot[t=0:1]\
 "GKTtype2b-GP.tab" u 1:4 t "" w l ls 1,\
 rhoc1, t*Qc1 t "Limit {/Symbol r}_{1}" w l ls 2,\
 rhoc2, t*Qc2 t "Limit {/Symbol r}_{2}" w l ls 3,\
 rhocconv, t*Qcconv t "Limit {/Symbol r}_{conv}" w l lt 0,\
 rhoc3, t*Qc3 t "Limit {/Symbol r}_{3}" w l ls 3,\
 rhoc4, t*Qc4 t "Limit {/Symbol r}_{4}" w l ls 2,\
 rhoc1+t*rholine, 1000+t*Qline w l lt 0,\
 rhoc2+t*rholine, 800+t*Qline w l lt 0,\
 rhoc3+t*rholine, 600+t*Qline w l lt 0,\
 rhoc4+t*rholine, 300+t*Qline w l lt 0

set nolabel 1

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

set out "GKTtype2b-GP.tseriesUpDown.png"

set auto x
set xrange [-3:]
set xlabel "Time (min)"
set yrange [0:160]
set ylabel "V (km/h)"
set key

plot[t=0:1]\
  "GKTtype2b-GP.x12000"  u ($2-60):4 t "10 km upstream  " w l ls 1,\
  "GKTtype2b-GP.x21000" u ($2-60):4 t "1 km upstream   " w l ls 2,\
  "GKTtype2b-GP.x23000" u ($2-60):4 t "1 km downstream " w l ls 3,\
  "GKTtype2b-GP.x25000" u ($2-60):4 t "3 km downstream " w l ls 3,\
  0, 110*t t "begin/end of obstruction" w l ls 10,\
  30, 110*t t "" w l ls 10

set out "GKTtype2b-GP.tseriesUp.png"
plot[t=0:1]\
  "GKTtype2b-GP.x21000" u ($2-60):4 t " 1 km upstream " w l ls 1,\
  "GKTtype2b-GP.x19000"  u ($2-60):4 t " 3 km upstream " w l ls 2,\
  "GKTtype2b-GP.x17000"  u ($2-60):4 t " 5 km upstream " w l ls 3,\
  "GKTtype2b-GP.x12000"  u ($2-60):4 t "10 km upstream " w l ls 5,\
  0, 110*t t "begin/end of obstruction" w l ls 10,\
  30, 110*t t "" w l ls 10

set auto y
set out "GKTtype2b-GP.aUp.png"

plot\
  "GKTtype2b-GP.x12000" u ($2-60):6 t "10 km upstream  " w l ls 1,\
  "GKTtype2b-GP.x21000" u ($2-60):6 t "1 km upstream   " w l ls 2,\
  "GKTtype2b-GP.x23000" u ($2-60):6 t "1 km downstream " w l ls 3,\
  "GKTtype2b-GP.x25000" u ($2-60):6 t "3 km downstream " w l ls 4

set auto y
set out "tmp.png"

plot\
  "GKTtype2b-GP.x12000"  u ($2-60):5 t "10 km upstream  " w l ls 1,\
  "GKTtype2b-GP.x21000" u ($2-60):5 t "1 km upstream   " w l ls 2,\
  "GKTtype2b-GP.x23000" u ($2-60):5 t "1 km downstream " w l ls 3,\
  "GKTtype2b-GP.x25000" u ($2-60):5 t "3 km downstream " w l ls 4

#********************************************************
# Dynamic trajectories in Q-rho space with 
# equilibrium fundamental diagram (along x=const and t=const)
#********************************************************

set out "GKTtype2b-GP.Qrho.png"

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
  "GKTtype2b-GP.x19000"  u (rho($3,$2)):(Q($5,$2)) t "3 km upstream " w l ls 5,\
  "GKTtype2b-GP.tab"   u 1:4 t "Equilibrium flow"  w l ls 11








