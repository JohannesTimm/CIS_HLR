set title "Strong Scaling Jacobi 960 Interlines"
set autoscale
set pointsize 2.5
set xlabel 'Prozesse'
set ylabel 'Knoten'
set zlabel 'Zeit'
set terminal postscript eps
set output "STRONG_SCALING_JA.eps"
set grid
splot "STRONG_SCALING_JA.dat" using 1:2:4 with linespoints
