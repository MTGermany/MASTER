set pointsize 1.2 #VOR definition der ls! 

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
#
set term png enhanced font Helvetica 14 xffffff
set size 1,1

######### 3d colored/contour

unset contour              # no contour lines
set contour surface      # Aktiviert Kontourlinien auf 3D-Flaeche
#set contour base          # Aktiviert Kontourlinien auf xy-Ebene
#set contour both           # Aktiviert Kontourlinien auf beidem
unset clabel  # dann lauter gleiche Kontourlinien; 
                     # Farbe und Typ mit "w l ls" beim splot-Kommando
                     # Achtung bug! explizit ("w l lt 3 lw 1" etc) gibt DOS!

set cntrparam bspline 

set palette defined ( 0 "#cc00ff", 0.15  "red", 0.25 "orange",\
0.35 "yellow", 0.55 "green", 0.8 "blue", 1 "#2200aa")
#set cbrange [0:120]  #[min:max] of color-coded  range (default zrange)

set pm3d
set style line 99 lt 7 lw 0 linecolor rgb "#000000" # beliebige Farben:Schwarz
set pm3d  hidden3d 99                   # color coded 3d surface
set contour surface      # Aktiviert Kontourlinien auf 3D-Flaeche
#unset surface      # pm3d off: kein Gitternetz gemacht
                            # pm3d on: Schnelles Plotten mit/ohne Netz
                            # (je nach pm3d Einstellung) mit Artefakten
set view 20,70
set nogrid

set out "GKT_Marathon_rho3d.png"
print "plotting GKT_Marathon_rho3d.png ..."
set xrange [0:18]
set auto y
set auto z
set xlabel "x (km)"
set ylabel "t (min)"
set zlabel "{/Symbol r} (runners/m)"
splot  "GKT_Marathon.dat" u ($1):2:($3*0.001)   w l lt 6

set out "GKT_Marathon_v3d.png"
print "plotting GKT_Marathon_v3d.png ..."
set xrange [0:18]
set auto y
set auto z
set xlabel "x (km)"
set ylabel "t (min)"
set zlabel "V (km/h)"
splot  "GKT_Marathon.dat" u ($1):2:4   w l lt 6

######################
# Fundamental diagram
######################

unset colorbox
set size 1,1
set out "GKT_Marathon_fund.png"
print "plotting GKT_Marathon_fund.png .."
set xlabel "density (1/km)"
set auto x
set ylabel "Q (1/h)" offset 0.5,0

set key
set auto y

plot\
  "GKT_Marathon.tab"   u 1:4 t "Equilibrium" w l ls 1,\
  "GKT_Marathon.x3500"  u ($3):5 t "x=3500 m" w linesp ls 12,\
  "GKT_Marathon.x5000"  u ($3):5 t "x=5000 m" w linesp ls 13,\
  "GKT_Marathon.x10000"  u ($3):5 t "x=10000 m" w linesp ls 14




########################################
set out "GKT_Marathon_Q.png"
print "plotting GKT_Marathon_Q.png ..."
########################################

set auto
set size 1,0.9
set xlabel "t (min)" offset 0,0.5
set ylabel "Q (runners/h)"

set auto x
set auto y

plot\
 "GKT_Marathon.x500"   u ($2) : 5 t "x10"    w l ls 1,\
 "GKT_Marathon.x500"   u ($2) : 5 t "x500"   w l ls 2,\
 "GKT_Marathon.x1000"  u ($2) : 5 t "x1000"  w l ls 3,\
 "GKT_Marathon.x2500"  u ($2) : 5 t "x2500"  w l ls 4,\
 "GKT_Marathon.x5000"  u ($2) : 5 t "x5000"  w l ls 5,\
 "GKT_Marathon.x10000" u ($2) : 5 t "x10000" w l ls 6


########################################
set out "GKT_Marathon_rho.png"
print "plotting GKT_Marathon_rho.png ..."
########################################

set auto
set size 1,0.9
set xlabel "t (min)" offset 0,0.5
set ylabel "{/Symbol r} (runners/km)"

set auto x
set auto y

plot\
 "GKT_Marathon.x500"   u ($2) : 3 t "x10"    w l ls 1,\
 "GKT_Marathon.x500"   u ($2) : 3 t "x500"   w l ls 2,\
 "GKT_Marathon.x1000"  u ($2) : 3 t "x1000"  w l ls 3,\
 "GKT_Marathon.x2500"  u ($2) : 3 t "x2500"  w l ls 4,\
 "GKT_Marathon.x5000"  u ($2) : 3 t "x5000"  w l ls 5,\
 "GKT_Marathon.x10000" u ($2) : 3 t "x10000" w l ls 6







