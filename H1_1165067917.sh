#!/bin/bash

./whiten_data 0 256 1165067917

./gauss_test time_0_256_1165067917.dat 32.0 256.0

for ((j=10;j<16;j++))
do
echo $j
cp plot_$j.dat Binary.dat
foo=$(printf "Spectogram_%d_1165067917.dat" $j)
cp plot_$j.dat $foo
done

for ((j=10;j<16;j++))
do
echo $j
./stationary time_0_256_1165067917.dat $j 32.0 256.0
foo=$(printf "PowerTime_H1_1165067917_%d.dat" $j)
mv PowerTime.dat $foo
done

gnuplot power_fluc_1165067917.gnu

