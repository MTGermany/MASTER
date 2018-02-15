set pointsize 1.2 #VOR definition der ls! 

set style line 1 lt 1 lw 2 pt 7 ps 1.9  lc rgb "#000000"
set style line 2 lt 1 lw 2 pt 7 ps 1.5  lc rgb "#CC0022" #rot, solid Kreuz
set style line 3 lt 8 lw 2 pt 4 ps 1.2  lc rgb "#FF3300" #orange, openSquare
set style line 4 lt 6 lw 2 pt 4 ps 1.5  lc rgb "#FFAA00" #gelb, openSquare
set style line 5 lt 1 lw 2 pt 5 ps 1.5  lc rgb "#00DD22" #gruen, closedBox
set style line 6 lt 5 lw 2 pt 4 ps 1.5  lc rgb "#00AAAA" #offenes Quadrat
set style line 7 lt 3 lw 2 pt 4 ps 2.0  lc rgb "#1100FF" #blau,dottedSquare
set style line 8 lt 4 lw 2 pt 8 ps 1.5  lc rgb "#220088"
set style line 9 lt 7 lw 2 pt 9 ps 1.5  lc rgb "#999999" #aufrClosedTriang. 

set style line 11 lt 1 lw 6 pt 7 ps 1.9  lc rgb "#000000" 
set style line 12 lt 1 lw 6 pt 2 ps 1.5  lc rgb "#CC0022" 
set style line 13 lt 8 lw 6 pt 4 ps 1.2  lc rgb "#FF3300"
set style line 14 lt 6 lw 6 pt 4 ps 1.5  lc rgb "#FFAA00"
set style line 15 lt 1 lw 6 pt 5 ps 1.5  lc rgb "#00DD22"
set style line 16 lt 5 lw 6 pt 7 ps 1.5  lc rgb "#00AAAA"
set style line 17 lt 1 lw 6 pt 7 ps 1.5  lc rgb "#1100FF"
set style line 18 lt 4 lw 6 pt 8 ps 1.5  lc rgb "#220088"
set style line 19 lt 7 lw 6 pt 9 ps 1.5  lc rgb "#999999"

set term png enhanced font Helvetica 14 xffffff
set size 1,1


tmin=0
tmax=120
xmin=2.2
xmax=25
######### 3d colored/contour


set cntrparam levels discrete 10,20,30,40,50,60,70,80,90,100 # freely set lines
unset clabel  # dann lauter gleiche Kontourlinien; 
                     # Farbe und Typ mit "w l ls" beim splot-Kommando
set palette defined ( 0 "#cc00ff", 0.20  "red", 0.27 "orange",\
0.35 "yellow", 0.45 "green", 0.600 "blue", 1 "#2200aa")
#set cbrange [0:100]  #[min:max] of color-coded  range (default zrange)

set pm3d                     # color coded 3d surface
set contour surface      # Aktiviert Kontourlinien auf 3D-Flaeche
set nocontour
unset surface      # pm3d off: kein Gitternetz gemacht
                            # pm3d on: Schnelles Plotten mit/ohne Netz
                            # (je nach pm3d Einstellung) mit Artefakten
set view 20,70
set nogrid

set out "demo_hystereticFD_v3d.png"
print "plotting demo_hystereticFD_v3d.png ..."
set xrange [xmin:xmax]
set yrange [tmin:tmax]
set ytics 30
#set zrange [0:100] reverse
set ztics 50
#set xlabel "x (km)" offset 1,0
set xlabel ""
#set xtics 5
set noxtics 
set ylabel "t (min)"
set zlabel "V (km/h)" offset 0,0.5
splot  "demo_hystereticFD.dat" u ($1):2:4   w l lt 6

######################
# Spatiotemporal velocity map
######################

set out "demo_hystereticFD_v3d_map.png"
print "plotting demo_hystereticFD_v3d_map.png ..."


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

splot  "demo_hystereticFD.dat" u ($2):($1):4   w l lt 6


######################
# Fundamental diagram
######################

set term png enhanced font Helvetica 14 xffffff
set size 0.9,1

unset pm3d

set out "demo_hystereticFD_fund.png"
print "plotting demo_hystereticFD_fund.png .."
set xlabel "Dichte {/Symbol r} (Fz/km)"
set xtics 20
set ylabel "Fluss Q (Fz/h)" offset 0.5,0
set ytics 500
set key
set xrange [0:120]
set yrange [0:2500]

plot\
  "demo_hystereticFD.tab"   u 1:4 t "Fundamentaldiagramm" w l ls 1,\
  "demo_hystereticFD.x16000"  u ($3):5 t "x=16000 m" w linesp ls 2,\
  "demo_hystereticFD.x14000"  u ($3):5 t "x=14000 m" w linesp ls 3,\
  "demo_hystereticFD.x5000"  u ($3):5 t "x=5000 m" w linesp ls 7



set term png enhanced font Helvetica 14 xffffff
set auto
set size 1,0.4
set xlabel "t (min)" offset 0, 0.5
set auto x
set xtics 30
set ylabel "V (km/h)"
set auto y
set ytics 20

set out "demo_hystereticFD_v_x16000.png"
print "plotting demo_hystereticFD_v_x16000.png ..."
plot "demo_hystereticFD.x16000" u ($2) : 4 t "x16000" w l ls 1

