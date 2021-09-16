# Heterogenous-cortical-network
This is based on my project on firing statistics of heterogenous cluster sizes in cortical networks.
Written in C++ (also works with C) and plotting is done with gnuplot.

This project has 4 main C++/C codes:
1. homegeneous.cpp - Produces datafile for plotting a raster plot with homogeneous cluster sizes.
2. balanced_heterogeneous.cpp - Produces firing statistics for perfectly balanced network with heterogeneous cluster sizes.
3. partial_balance.cpp - Produces partially balanced network and its firing statistics.
4. hierarchy.cpp - Produces a partially balanced network with hierarchy.

The community sizes for the exciters are defined in the popexp.dat file. The sizes are exponentially distributed with mean = 80 (on average 50 clusters).

ran2.cpp is the random number generator (RAN2) used in all the scripts.

plot.plt plots the raster plot for the spiking time of each neuron with each community distinguished by alternative coloring scheme. This shows the variability/uniformity of firing of the neurons.  
