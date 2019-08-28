set term pdf enhanced color font 'Helvetica,11' 
set output "correlations.pdf" 
set key left
set size 1.0,1.0
set origin 0,0
set multiplot
set size 1.0, 0.52
set origin 0, 0.503
set xrange [-20:20]
set yrange [-1.0:1.0]
set ytics (-0.9,-0.6, -0.3, 0, 0.3, 0.6, 0.9)
set ylabel 'C({/Symbol t})'
unset xtics
set label "Fig 1 PRL" at graph 0.8,0.9
plot "fig1_first_corr.dat" using 1:2 title "First half" with lines lw 2 lc rgb "black", "fig1_last_corr.dat" using 1:2 title "Last half" with lines lw 2  lc rgb "#D1112E", "fig1_full_corr.dat" using 1:2 title "All 0.2s" with lines lw 2  lc rgb "#1180D1", "fig1_small_corr.dat" using 1:2 title "0.39s-0.43s" with lines lw 2  lc rgb "#319947", "line2.dat" using 1:2 notitle with lines lw 2 lc rgb "black" dashtype 3
set origin 0.0, 0.0
set xtics
unset label
set size 1.0, 0.57
set label "Whiten + Maximized" at graph 0.7,0.9
set xlabel '{/Symbol t} (ms)'
plot "res_corr_first.dat" using ($1*1000):2 title "First half" with lines lw 2 lc rgb "black", "res_corr_last.dat" using ($1*1000):2 title "Last half" with lines lw 2  lc rgb "#D1112E", "res_corr_full.dat" using ($1*1000):2 title "All 0.2s" with lines lw 2  lc rgb "#1180D1", "res_corr_small.dat" using ($1*1000):2 title "0.39s-0.43s" with lines lw 2  lc rgb "#319947", "line2.dat" using 1:2 notitle with lines lw 2 lc rgb "black" dashtype 3
unset multiplot
