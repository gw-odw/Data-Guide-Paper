set term png enhanced truecolor crop font 'Helvetica,20'  size 1200,800
set output "spectraH1.png"
set format "%g"
set logscale
set key top right
set xtics (32,64,128,256)
set ylabel "S_n(f)"
set xlabel "f (Hz)"
set xrange [20:400]
set yrange [1e-48:1e-41]
plot "freqnowindow_0_4_1126259462.dat" using 1:2 title "No window" with lines, "pspec_0_4_1126259462.dat" using 1:2 title "Tukey window" with lines, "Welch_0_4_1126259462.dat" using 1:2 title "Welch average" with lines lt -1, "Welch_0_4_1126259462.dat" using 1:(1.2e-40/($1*$1)) title "1/f^2" with lines lc rgb "red" lw 2



