
# KKKL


set term png enhanced font Helvetica 14 xffffff
set out "./KKKL.1.png"

set view 50,60
splot "KKKL.dat" u 1:2:3 w l lt 2
