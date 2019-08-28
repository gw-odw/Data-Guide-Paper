#!/bin/bash

#These codes and scripts make Figures 2, 3, 4, 5, 6, 7, 8, 9, 13, 14, 16, 17, 18, 19

#These codes do not make Figues 1, 10, 11, 12, 15.
#Figues 12 and 15 can be made using LALinference https://github.com/lscsoft/lalsuite/tree/master/lalinference

echo "Compiling the codes"

g++ -O2 -o gauss_test gauss_test.cc Biorthogonal.cc Haar.cc Meyer.cc Symlet.cc WaveDWT.cc Wavelet.cc wavearray.cc wavecomplex.cc wseries.cc wavefft.cc
g++ -O2 -o stationary stationary.cc Biorthogonal.cc Haar.cc Meyer.cc Symlet.cc WaveDWT.cc Wavelet.cc wavearray.cc wavecomplex.cc wseries.cc wavefft.cc
gcc -o dataguide dataguide.c -lm -lgsl
gcc -o fig1_corr fig1_corr.c -lm -lgsl
gcc -o CBCresiduals CBCresiduals.c -lm -lgsl
gcc -o Qscans Qscans.c -lm -lgsl
gcc -o Gauss_Test_Freq Gauss_Test_Freq.c -lm -lgsl
gcc -o ReIm_kde ReIm_kde.c -lm -lgsl
gcc -o whiten_data whiten_data.c -lm -lgsl


# Uses GWOSC data and template for GW150914. Illustrates windowing, bandpass, maximum likelihood
# subtraction and computes correlations in the data and residuals
# Need to first download and unzip the files
# https://www.gw-openscience.org/GW150914data/H-H1_LOSC_4_V2-1126257414-4096.txt.gz
# https://www.gw-openscience.org/GW150914data/L-L1_LOSC_4_V2-1126257414-4096.txt.gz
# https://www.gw-openscience.org/GW150914data/P150914/fig2-unfiltered-waveform-H.txt
./dataguide

# makes Figure 2
gnuplot timeH1.gnu

# makes Figure 3
gnuplot spectraH1.gnu

# makes Figure 4
gnuplot phasesH1.gnu

# makes Figure 16
gnuplot templatesNR.gnu

# makes Figure 17
gnuplot residuals.gnu

# makes Figure 6 Using the output from several BayesWave runs. Code at https://git.ligo.org/lscsoft/bayeswave
gnuplot psds.gpi

# read in LIGO data and use a median + line catcher to compute PSD
# This can take a few minutes for 256 seconds of data since not using a running median
# Have to first download 256s of L1 data centered on GPS time 1186741733 and save it to
# an ASCII file called frame_256_1186741733_1.dat
# Basic instructions given here: https://www.gw-openscience.org/archive/O1/
echo " "
echo "Processing the data to make Figure 5 takes a while. You might want to go and get a coffee"
echo " "
./whiten_data 1 256 1186741733

# Test Gaussiantiy of whitened data
./Gauss_Test_Freq freq_1_256_1186741733.dat

# KDE of 2-d real and imaginary parts
./ReIm_kde

# makes Figure 5
gnuplot Gauss_check.gnu

# Uses the GWOSC data for Figure 1 to computes correlations in the data and residuals
# Also uses the redisuals and whitened data produced by dataguide.c
# Need to first download the GWOSC files
#https://www.gw-openscience.org/GW150914data/P150914/fig1-observed-H.txt
#https://www.gw-openscience.org/GW150914data/P150914/fig1-observed-L.txt
#https://www.gw-openscience.org/GW150914data/P150914/fig1-residual-H.txt
#https://www.gw-openscience.org/GW150914data/P150914/fig1-residual-L.txt
./fig1_corr

# makes Figure 18 using output from dataguide.c and fig1_corr.c
gnuplot data_correlations.gnu

# makes Figure 19 using output from dataguide.c and fig1_corr.c
gnuplot correlations.gnu

# reads in the original CBC residuals used for the BayesWave TGR residual SNR test
# (files are clean_frequency_data_199.dat.0/1)
./CBCresiduals

# makes Figure 14 using the output of the CBCresiduals.c code
gnuplot Gauss_Norm.gnu

# reads in the original CBC residuals used for the BayesWave TGR residual SNR test
# (files are clean_frequency_data_199.dat.0/1 and clean_frequency_residual_199.dat.0,1)
# and the original data (with the signal), which are whitened by the BW PSD.
# The whitend data is then Q-transformed for plotting
./Qscans

# The following gnuplot calls uses the output from Qscans.c to make the files H1raw.png, L1raw.png,
# Hresidual.png, Lresidual.png, which are then stiched together to produce Figure 13
gnuplot H1_raw.gnu
gnuplot L1_raw.gnu
gnuplot H1_residual.gnu
gnuplot L1_residual.gnu

echo " "
echo "Making figure 7,8,9. These take a while, you might want to go get lunch then come back"
echo " "

#makes figure 7
#Have to first download 256s of H1 data centered on GPS time 1165067917 and save it to
# an ASCII file called frame_256_1165067917_1.dat
source H1_1165067917.sh

#Have to first download 256s of L1 data centered on GPS time 1186741733 and save it to
# an ASCII file called frame_256_1186741733_1.dat
#makes figure 8
source L1_1186741733.sh

#Have to first download 128s of L1 data centered on GPS time 1166358283 and save it to
# an ASCII file called frame_128_1166358283_1.dat
#makes figure 9
source L1_1166358283.sh



