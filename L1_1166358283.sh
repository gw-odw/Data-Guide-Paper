#!/bin/bash

./whiten_data 1 128 1166358283

./gauss_test time_1_128_1166358283.dat 32.0 256.0

for ((j=11;j<17;j++))
do
echo $j
cp plot_$j.dat Binary.dat
foo=$(printf "Spectogram_%d_1166358283.dat" $j)
cp plot_$j.dat $foo
done

for ((j=11;j<17;j++))
do
echo $j
./stationary time_1_128_1166358283.dat $j 32.0 256.0
foo=$(printf "PowerTime_L1_1166358283_%d.dat" $j)
mv PowerTime.dat $foo
done

gnuplot power_normal_1166358283.gnu

