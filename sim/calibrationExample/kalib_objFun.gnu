
set term png enhanced font Helvetica 14 xffffff
set out "kalib_objFun.png"
print "plotting kalib_objFun.png"

set style line 99 lt 7 lw 0 linecolor rgb "#000000" # beliebige Farben:Schwarz

unset pm3d                            # no color coding
set pm3d                                # color coded 3d surface
set pm3d  hidden3d 99         # color coded 3d surface with grid, ls 99
set view 20,330

set palette defined ( 0 "#dd00ff", 10 "blue", \
      20 "green", 30 "yellow", 80 "orange", 100 "red")

#danach

unset colorbox
min(x,y)=(x<y) ? x : y
splot "kalib_objFun.dat" u 1:2:(min($3,6))  w l ls 99

