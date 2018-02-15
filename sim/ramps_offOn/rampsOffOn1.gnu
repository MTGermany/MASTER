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

set style line 11 lt 7 lw 1 #schwarz
set style line 12 lt 3 lw 1 pt 12 #blau
set style line 13 lt 4 lw 1 #lila
set style line 14 lt 1 lw 1 #rot
set style line 15 lt 2 lw 1 #gruen
set style line 16 lt 6 lw 1

set term png enhanced font Helvetica 14 xffffff
set size 1,1

######### 3d colored/contour


set cntrparam levels discrete 10,20,30 # freely set lines
unset clabel  # dann lauter gleiche Kontourlinien; 
                     # Farbe und Typ mit "w l ls" beim splot-Kommando
set palette defined ( 0 "#cc00ff", 0.15  "red", 0.25 "orange",\
0.35 "yellow", 0.55 "green", 0.8 "blue", 1 "#2200aa")
#set cbrange [0:120]  #[min:max] of color-coded  range (default zrange)

set pm3d                     # color coded 3d surface
set contour surface      # Aktiviert Kontourlinien auf 3D-Flaeche
unset surface      # pm3d off: kein Gitternetz gemacht
                            # pm3d on: Schnelles Plotten mit/ohne Netz
                            # (je nach pm3d Einstellung) mit Artefakten
set view 20,70
set nogrid

set out "rampsOffOn1_v3d.png"
print "plotting rampsOffOn1_v3d.png ..."
set xrange [2:18]
set yrange [50:]
set zrange [:] reverse
set xlabel "x (km)"
set ylabel "t (min)"
set zlabel "V (km/h)"
splot  "rampsOffOn1.dat" u ($1):2:4   w l lt 6

######################
# Fundamental diagram
######################

set out "rampsOffOn1_fund.png"
print "plotting rampsOffOn1_fund.png .."
set xlabel "density (1/km)"
set ylabel "Q (1/h)" offset 0.5,0

set key
set xrange [0:80]
set auto y
set size 0.8,1
plot\
  "rampsOffOn1.tab"   u 1:4 t "Equilibrium" w l ls 1,\
  "rampsOffOn1.x15000"  u ($3):5 t "x=15000 m" w linesp ls 12,\
  "rampsOffOn1.x14000"  u ($3):5 t "x=14000 m" w linesp ls 13,\
  "rampsOffOn1.x12000"  u ($3):5 t "x=12000 m" w linesp ls 14,\
  "rampsOffOn1.x10000"  u ($3):5 t "x=10000 m" w linesp ls 15,\
  "rampsOffOn1.x5000"  u ($3):5 t "x=5000 m" w linesp ls 16


set term png enhanced font Helvetica 14 xffffff
set auto
set size 1,0.4
set xlabel "t (min)" offset 0,0.5
set ylabel "V (km/h)"
set ytics 20
set auto x
set auto y

set out "rampsOffOn1_v_x16000.png"
print "plotting rampsOffOn1_v_x16000.png ..."
plot "rampsOffOn1.x16000" u ($2) : 4 t "x16000" w l ls 1

set out "rampsOffOn1_v_x15500.png"
print "plotting rampsOffOn1_v_x15500.png ..."
plot "rampsOffOn1.x15500" u ($2) : 4 t "x15500" w l ls 1

set out "rampsOffOn1_v_x15000.png"
print "plotting rampsOffOn1_v_x15000.png ..."
plot "rampsOffOn1.x15000" u ($2) : 4 t "x15000" w l ls 1

set out "rampsOffOn1_v_x14000.png"
print "plotting rampsOffOn1_v_x14000.png ..."
plot "rampsOffOn1.x14000" u ($2) : 4 t "x14000" w l ls 1

set out "rampsOffOn1_v_x12000.png"
print "plotting rampsOffOn1_v_x12000.png ..."
plot "rampsOffOn1.x12000" u ($2) : 4 t "x12000" w l ls 1

set out "rampsOffOn1_v_x10000.png"
print "plotting rampsOffOn1_v_x10000.png ..."
plot "rampsOffOn1.x10000" u ($2) : 4 t "x10000" w l ls 1

set out "rampsOffOn1_v_x5000.png"
print "plotting rampsOffOn1_v_x5000.png ..."
plot "rampsOffOn1.x5000" u ($2) : 4 t "x5000" w l ls 1







