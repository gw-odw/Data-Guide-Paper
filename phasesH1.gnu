set term png enhanced truecolor crop font Helvetica 20  size 1200,800
set output "phasesH1.png"
set format "%g"
set key box opaque
#set xtics (32,64,128,256)
set ylabel "Phase"
set xlabel "f (Hz)"
set xrange [20:400]
set yrange [0:6.2831853071795865]
plot "freqnowindowres_0_4_1126259462.dat" using 1:3 title "No window" with points pt 9 ps 0.7 lc rgb "blue", "freqres_0_4_1126259462.dat" using 1:3 title "Tukey window" with points pt 5 ps 0.7 lc rgb "red"