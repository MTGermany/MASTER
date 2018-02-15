#http://t16web.lanl.gov/Kawano/gnuplot/gallery/func2.html

set style line 1 lt 7 lw 3 pt 7 ps 0.5

########### Aufloesung, Zahl der Contourlinien ####

#set isosample 40,40
#set samples 101
set isosamples 20



#################################
# Contourlinien 
#################################

#set contour surface # Aktiviert Kontourlinien auf 3D-Flaeche
#set contour base  # Aktiviert Kontourlinien auf xy-Ebene
set cntrparam bspline 
set cntrparam levels 10
unset clabel  # dann lauter gleiche Kontourlinien; Bug: nur lt 1=rot waehlbar!

#################################
# Color-Coding!
#################################

#set pm3d map  # Aktiviert color-coding auf Ebene
set pm3d           # Aktiviert color-coding auf 3D-Plot


set palette color
set palette model RGB
# set palette defined (z1, farbe1, ...), wobei zi immer automat. auf
# min/max normiert werden 
# tatsaechl. min/max. unabhaengig von zrange mit "set cbrange" steuerbar
set palette defined ( 0 "red", 15 "orange", 30 "yellow", 45 "green",\
                               70 "blue", 100 "#dd00ff" ) 
set cbrange [0:140]  #Minimum/Maximum der Palette

# set palette rgbformulae -33,-13,31 #obtained with show palette fit2rgbformulae


#################################
# Gitternetz
#################################

set style line 99 lt 7 lw 0
set pm3d hidden3d 99 # Color-coded und Gitternetz
                                       # braucht einen linestyle, hier 99

# Bug: nur mit "set surface" hidden3d bezueglich Koordinatengitters aktiv!
# Achtung: Braucht lang, v.a. in Verbindung mit Kontourlinien
# Das eigentliche Gitternetz wird von pm3d uenerschrieben, aber die
# stoerenden Koordinatenlinien sind weg 
# (klappt nicht mit Punkten, deshalb "w l ls x" notwendig)
# Um tatsaechlich Gitternetze auf die Farbcodierung zu zeichnen:
# set pm3d hidden3d 99   (siehe oben)

set surface         # setzt Gitternetz bzw. Punktewolke 
                          # wenn ohne "w l ls x"
                          # (beides wird allerdings von pm3d  ueberschrieben)


#################################
# eigentliches Plotten
#################################

unset surface  # falls schnelles Plotten mit Artefakten

set term png enhanced font Helvetica 14 xffffff
set out "rampPinch.v3d.png"
#set out "|awk -f pm3dConvertToImage.awk >coloredContourData.png"
print "plotting rampPinch.v3d.png ..."

set nogrid  #ansonsten 2D-Gitter auf xy-Ebene

#set size square
set size 1,1



set xlabel "x (km)"
set ylabel "Time (h)"
set zlabel "V (km/h)"
set xrange [4:18]
#set auto x
set yrange [40:160]
set zrange [0:]reverse
set ztics 50

set view 20,70

splot  "rampPinch.dat" u 1:2:4   w l ls 99
 # linestyles z.B. w l ls 3 ergeben DOS