# ---  
clear
reset
set title 'Wang-Landau (tupla fractions, unnormalized)'  font 'Times-Roman, 26'
set xlabel 'E/N' font 'Times-Roman, 26'
set ylabel 'ln(g(E))' font 'Times-Roman, 26'
set xtics font 'Times-Roman, 22'
set ytics font 'Times-Roman, 22'

set key inside top vertical right
set key font 'Times-Roman, 20'

set size 0.8, 1
N = 100
N2 = N*N
# -----------------------------------------------------------------------------

set yrange [7050:7350]

plot  'wl_2_N0100_a0.20log.csv' u (($1)/N2):($2) notitle ls 1


