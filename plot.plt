filename = "./spikeraster.dat"  # need the same file several times

stats filename             # get number of blocks
show variables             # check STATS_blocks
set term png size 1024,768
unset key
set xlabel "Time (ms)" font ",20"
set ylabel "Neuron Index" font ",20"
set tics font ",16"

set pointsize 0.18
set style line 1 lc rgb '#0000FF' pt 7
set style line 2 lc rgb '#FF0000' pt 7
set output "stim_perfect.png"
plot for [b=0:STATS_blocks-2] filename u 2:1 index b notitle w p ls (b%2+1)
unset output
