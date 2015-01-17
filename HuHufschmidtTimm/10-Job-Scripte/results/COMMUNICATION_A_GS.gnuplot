set title "Kommunikationstest A Gauss-Seidel 200 Interlines"
set autoscale
set pointsize 2.5
set xlabel 'Knoten'
set ylabel 'Zeit'
set terminal postscript eps
set output "COMMUNICATION_A_GS.eps"
set grid
plot "COMMUNICATION_A_GS.dat" using 2:4 with linespoints
