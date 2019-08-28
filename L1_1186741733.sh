#!/bin/bash

#./whiten_data 1 256 1186741733  Already done in main script

./gauss_test time_1_256_1186741733.dat 16.0 512.0

for ((j=10;j<16;j++))
do
echo $j
cp plot_$j.dat Binary.dat
foo=$(printf "Spectogram_%d_1186741733.dat" $j)
cp plot_$j.dat $foo
gnuplot discrete_amp.gnu
gnuplot discrete_power.gnu
foo=$(printf "L1_Amp_%d_1186741733.png" $j)
cp Amp.png $foo
foo=$(printf "L1_Pow_%d_1186741733.png" $j)
cp Pow.png $foo
done

for ((j=10;j<16;j++))
do
echo $j
./stationary time_1_256_1186741733.dat $j 16.0 510.0
foo=$(printf "PowerTime_L1_1186741733_%d.dat" $j)
mv PowerTime.dat $foo
done

gnuplot power_normal_1186741733.gnu
