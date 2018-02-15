
# VMM


set term png enhanced font Helvetica 14 xffffff
set out "./VMM.rho3d.png"

set view 50,60
splot "VMM.dat" u 1:2:3 w l lt 2
