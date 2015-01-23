set title "Weak Scaling Gauss-Seidel"
set autoscale
set pointsize 2.5
set xlabel 'Prozesse'
set ylabel 'Knoten'
set zlabel 'Zeit'
set cblabel "Interlines"
set palette defined ( 0 "blue", 3 "green", 6 "yellow", 10 "red" )
set terminal postscript eps
set output "WEAK_SCALING_GS.eps"
set grid
splot "WEAK_SCALING_GS.dat" using 1:2:4:3 with linespoints palette
