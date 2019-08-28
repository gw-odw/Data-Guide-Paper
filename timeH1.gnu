set term pngcairo enhanced truecolor crop font 'Helvetica,24'  size 1600,1600
set output "timeH1.png"
set xrange [-2:2]
unset xtics
set ylabel "d(t)"
unset xtics
set multiplot
set label "Hanford Data" at graph 0.02,0.92
set label "{/Symbol \264} 10^{-19}" at graph 0,1.065 left
set size 1, 0.24
set origin 0.0, 0.72
set yrange [-7:7]
set ytics (-6, -3, 0, 3, 6)
plot "raw_0_4_1126259462.dat" using 1:($2*1e19) notitle with lines lc rgb "#1180D1" lw 2
unset label
set size 1, 0.24
set origin 0.0, 0.49
set label "Windowed" at graph 0.02,0.92
set label "{/Symbol \264} 10^{-19}" at graph 0,1.065 left
plot "windowed_0_4_1126259462.dat" using 1:($2*1e19) notitle with lines lc rgb "#319947" lw 2
unset label
set ylabel "d_w(t)"
set size 1, 0.24
set origin 0.0, 0.265
set label "Whitened" at graph 0.02,0.92
set yrange [-6:6]
set ytics (-6,-3, 0, 3, 6)
plot "time_0_4_1126259462.dat" using 1:2 notitle with lines lc rgb "#D1112E" lw 2
unset label
set size 1, 0.28
set origin 0.0, 0.0
set label "Bandpassed" at graph 0.02,0.92
set yrange [-8:8]
set ytics (-6, -3, 0, 3, 6)
set xtics
set ylabel "d_w(t)"
set xlabel "t (sec)"
plot "timebp_0_4_1126259462.dat" using 1:2 notitle with lines lc rgb "black" lw 2





