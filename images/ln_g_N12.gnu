# ---  
clear
reset
set title '{/Symbol a} = 0.2, N = 12' font 'Times-Roman, 28'
set xlabel 'E/N' font 'Times-Roman, 26'
set ylabel 'ln(g(E))' font 'Times-Roman, 26'
set xtics font 'Times-Roman, 22'
set ytics font 'Times-Roman, 22'

set key inside top vertical right
set key font 'Times-Roman, 20'

set size 0.8, 1
N = 12

# -----------------------------------------------------------------------------

plot  'wl_N0012_a0.20log.csv' u (($1)/N):($2) notitle ls 1


