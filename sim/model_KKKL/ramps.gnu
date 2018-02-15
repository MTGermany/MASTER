set pointsize 1.2 #VOR definition der ls! 

#set style line 1 lt 1 lw 2 pt 7 ps 1.9  lc rgb "#000000"
set style line 1 lt 1 lw 2 pt 7 ps 1.9  lc rgb "#000000"
set style line 2 lt 3 lw 6 pt 12 #blau
set style line 3 lt 4 lw 6 #lila
set style line 4 lt 1 lw 6 pt 10 #rot
set style line 5 lt 6 lw 6 #gelb  (screen: tausche mit ls 1=orange-ocker auf screen)
set style line 6 lt 2 lw 6 pt 14 ps 0.6 #gruen
set style line 7 lt 7 lw 6 #schwarz
set style line 8 lt 8 lw 6 #blassrot
set style line 9 lt 1 lw 6 pt 1 ps 1.5  lc rgb "#888888" #grau,solid,plus sign

set style line 11 lt 7 lw 1 #schwarz
set style line 12 lt 3 lw 1 pt 12 ps 1.5#blau
set style line 13 lt 4 lw 1 ps 1.5 #lila
set style line 14 lt 1 lw 1 #rot
set style line 15 lt 2 lw 1 #gruen
set style line 16 lt 6 lw 1

max(x,y)=(x>y) ? x : y

set term png enhanced font Helvetica 14 xffffff
set size 1,1

tmin=60
tmax=180
xmin=2.2
xmax=20
######### 3d colored/contour


set cntrparam levels discrete 10,20,30,40,50,60,70,80,90,100 # freely set lines
unset clabel  # dann lauter gleiche Kontourlinien; 
                     # Farbe und Typ mit "w l ls" beim splot-Kommando
set palette defined ( 0 "#cc00ff", 0.20  "red", 0.27 "orange",\
0.35 "yellow", 0.45 "green", 0.600 "blue", 1 "#2200aa")
#set cbrange [0:120]  #[min:max] of color-coded  range (default zrange)

set pm3d                     # color coded 3d surface
set contour surface      # Aktiviert Kontourlinien auf 3D-Flaeche
set nocontour
unset surface      # pm3d off: kein Gitternetz gemacht
                            # pm3d on: Schnelles Plotten mit/ohne Netz
                            # (je nach pm3d Einstellung) mit Artefakten
set view 20,70
set nogrid

set out "ramps_v3d.png"
print "plotting ramps_v3d.png ..."
set xrange [xmin:xmax]
set yrange [tmin:tmax]
set ytics 30

set ztics 50
#set xlabel "x (km)" offset 1,0
set xlabel ""
#set xtics 5
set noxtics 
set ylabel "t (min)"
set zlabel "V (km/h)" offset 0,0.5
splot  "ramps.dat" u ($1):2:4   w l lt 6

######################
# Spatiotemporal velocity map
######################

set out "ramps_v3d_map_color.png"
print "plotting ramps_v3d_map_color.png ..."

tmin=30
tmax=150
xmin=2.2
xmax=18

set view map
set contour
set cntrparam levels discrete 10,20,30,40,50,60 # freely set lines
set palette defined ( 0 "#cc00aa", 0.20  "red", 0.27 "orange",\
0.38 "yellow", 0.50 "green", 0.700 "blue", 1 "#440077")

set xlabel "t (min)"
set xrange [tmin:tmax]
set xtics 30
set ylabel "x (km)"
set yrange [xmin:xmax]
set ytics 2
unset zlabel
set label 1 "V (km/h)" at screen 1.26,1.8

splot  "ramps.dat" u ($2):($1):(max(0,$4))   w l lt 6

######################
set out "ramps_v3d_map.png"
print "plotting ramps_v3d_map.png ..."
######################

#set palette defined ( 0  "orange",\
#0.20 "yellow", 0.35 "green", 0.55 "blue", 0.70 "#440077", 1 "#440022")
set palette defined ( 0  "#440044", 0.1  "#770055", 0.3 "red", 0.45 "orange",\
0.55 "yellow", 0.65  "#99ff55", 1 "#bbffbb")
splot  "ramps.dat" u ($2):($1):4   w l lt 6

######################
# Fundamental diagram
######################

set term png enhanced font Helvetica 14 xffffff
set size 0.85,1

unset pm3d

set out "ramps_fund.png"
print "plotting ramps_fund.png .."
set xlabel "Dichte {/Symbol r} (1/km)"
set xtics 20
set ylabel "Fluss Q (Fz/h)" offset 0.5,0
set ytics 500
set key
set xrange [0:120]
set auto y
# set yrange [0:2500]
plot\
  "ramps.tab"   u 1:4 t "Fundamentaldiagramm" w l ls 9,\
  "ramps.x14000"  u ($3):5 t "x=14000 m" w linesp ls 12,\
  "ramps.x12000"  u ($3):5 t "x=12000 m" w linesp ls 13,\
  "ramps.x10000"  u ($3):5 t "x=10000 m" w linesp ls 11,\
  "ramps.x6000"  u ($3):5 t "x=6000 m" w linesp ls 10

quit

set term png enhanced font Helvetica 14 xffffff
set auto
set size 1,0.4
set xlabel "t (min)" offset 0, 0.5
set auto x
set xtics 30
set ylabel "V (km/h)"
set auto y
set ytics 20

set out "ramps_v_x9000.png"
print "plotting ramps_v_x9000.png ..."
plot "ramps.x9000" u ($2) : 4 t "x9000" w l ls 1
