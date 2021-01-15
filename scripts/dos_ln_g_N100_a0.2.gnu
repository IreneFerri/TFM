# ---  
clear
reset
set title 'Exact Solution {/Symbol a} = 0.2, N = 100' font 'Times-Roman, 26'
set xlabel 'E/N' font 'Times-Roman, 26'
set ylabel 'ln(g(E))' font 'Times-Roman, 26'
set xtics font 'Times-Roman, 22'
set ytics font 'Times-Roman, 22'

set key inside top vertical right
set key font 'Times-Roman, 20'

set size 0.8, 1
N = 100

# -----------------------------------------------------------------------------

plot  'dos_N100_a0.2.csv' u (($1)/N**2):($2) notitle ls 1


