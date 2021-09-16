#set term postscript eps enhanced color dl 0.3 blacktext "Helvetica" 20
set term png size 640,480
set output 'homogeneous.png'
unset key
set size square
set xtics 500
set xlabel "Time (in ms)"
set ylabel "Neuron index"

plot 'spikeraster.dat' u 2:1 w p pt 7 ps 0.15 lc rgb 'black'

unset output
