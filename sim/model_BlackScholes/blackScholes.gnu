
set style line 1 lt 1 lw 2 pt 7 ps 1.9  lc rgb "#000000"
set style line 2 lt 1 lw 5 pt 5 ps 1.5 #rot
set style line 3 lt 8 lw 5 pt 5 ps 1.5 #blassrot
set style line 4 lt 6 lw 5 pt 5 ps 1.5 #gelb  (screen: tausche mit lt 1=orange-ocker)
set style line 5 lt 2 lw 5 pt 5 ps 1.5 #gruen
set style line 6 lt 5 lw 5 pt 5 ps 1.5 #blasstuerkisblau
set style line 7 lt 3 lw 5 pt 5 ps 1.5 #blau
set style line 8 lt 4 lw 5 pt 5 ps 1.5 #lila

############## Europ. option ################################

set title "Discountzertifikat, strike 100, {/Symbol s}=0.2/sqrt{(year),r0=2%"
xmin=0        # minimum stock price!
xminPlot=50   # minimum plotted stock price
xmaxPlot=150  # max. plotted stock price


set term png enhanced font Helvetica 14 xffffff

set out "./blackScholes.1.png"

set view 50,10
set xrange [xminPlot:xmaxPlot]


set xlabel "Stock Price"
set ylabel "Time(years)"
set zlabel "Option Price"
splot "blackScholes.dat" u ($1+xmin):2:3 w l lt 2

#############################################################

set term png enhanced font Helvetica 14 xffffff
set out "./blackScholes.2.png"

set key at screen 0.55,0.8

set xrange [0:2]
set auto y

set xlabel "Time(years)"
set ylabel "Option price"

plot\
 "blackScholes.x80"  u 2:3 t "stock price xmin+  80 Eur" w l ls 2,\
 "blackScholes.x100" u 2:3 t "stock price xmin+100 Eur" w l ls 5,\
 "blackScholes.x120" u 2:3 t "stock price xmin+120 Eur" w l ls 7,\
 "blackScholes.x160" u 2:3 t "stock price xmin+160 Eur" w l ls 8


#############################################################

set out "./blackScholes.3.png"
set auto

set key at screen 0.4,0.8

set xrange [xminPlot:xmaxPlot]
set xlabel "Stock price"
set ylabel "Option price"

plot\
 "blackScholes.t0" u ($1+xmin):3 t "tau=0 Years" w l ls 1,\
 "blackScholes.t1" u ($1+xmin):3 t "tau=1 Years" w l ls 2,\
 "blackScholes.t2" u ($1+xmin):3 t "tau=2 Years" w l ls 3

