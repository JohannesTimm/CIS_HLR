#!/bin/sh 

gnuplot -persist STRONG_SCALING_JA.gnuplot &
gnuplot -persist WEAK_SCALING_JA.gnuplot &
gnuplot -persist COMMUNICATION_A_JA.gnuplot &
gnuplot -persist STRONG_SCALING_GS.gnuplot &
gnuplot -persist WEAK_SCALING_GS.gnuplot &
gnuplot -persist COMMUNICATION_A_GS.gnuplot &
wait

