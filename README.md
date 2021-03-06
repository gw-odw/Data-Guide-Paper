# Data-Guide-Paper
Codes used in the paper "Guide to  LIGO-Virgo detector noise and extraction of transient gravitational wave signals"

The codes and scripts in this repository make Figures 2, 3, 4, 5, 6, 7, 8, 9, 13, 14, 16, 17, 18, 19 of the data guide paper.

There is one maser script called dataguide.sh that compiles and runs all the codes and plotting scripts. Assuming that you have the gcc compiler and gnuplot installed you can simply type

source dataguide.sh

Before running the script you will need to download public LIGO data from the GWOSC website. The files are:

https://www.gw-openscience.org/GW150914data/H-H1_LOSC_4_V2-1126257414-4096.txt.gz
https://www.gw-openscience.org/GW150914data/L-L1_LOSC_4_V2-1126257414-4096.txt.gz
https://www.gw-openscience.org/GW150914data/P150914/fig2-unfiltered-waveform-H.txt
https://www.gw-openscience.org/GW150914data/P150914/fig1-observed-H.txt
https://www.gw-openscience.org/GW150914data/P150914/fig1-observed-L.txt
https://www.gw-openscience.org/GW150914data/P150914/fig1-residual-H.txt
https://www.gw-openscience.org/GW150914data/P150914/fig1-residual-L.txt

You also have to download 256s of H1 data centered on GPS time 1165067917 and save it to an ASCII file called frame_256_1165067917_1.dat, 256s of L1 data centered on GPS time 1186741733 and save it to an ASCII file called frame_256_1186741733_1.dat and 128 seconds of L1 data centered on GPS time 1166358283 and save it to an ASCII file called frame_128_1166358283_1.dat
