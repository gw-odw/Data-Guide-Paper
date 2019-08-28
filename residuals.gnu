set term png enhanced truecolor crop font Helvetica 24  size 1600,1600
set output "residuals.png"
set xrange [0.25:0.46]
set key top left
unset xtics
set ylabel "h(t)"
unset xtics
set multiplot
set yrange [-8:8]
set ytics (-6, -3, 0, 3, 6)
set size 1.0, 0.27
set origin 0.0, 0.74
plot "res.dat" using 1:2 title "H1 data" with lines lc rgb "blue" lw 3, "res.dat" using 1:3 title "H1 model" with lines lc rgb "red" lw 3
set size 1.0, 0.27
set origin 0.0, 0.51
set yrange [-8:8]
set ytics (-6, -3, 0, 3, 6)
plot "res.dat" using 1:4 title "L1 data" with lines lc rgb "blue" lw 3, "res.dat" using 1:5 title "L1 model" with lines lc rgb "red" lw 3
set size 1.0, 0.27
set origin 0.0, 0.28
set yrange [-4:4]
set ytics (-3, 0, 3)
plot "res.dat" using 1:($2-$3) title "H1 residual" with lines lc rgb "black" lw 3
set size 1.0, 0.32
set origin 0.0, 0.0
set yrange [-4:4]
set ytics (-3, 0, 3)
set xtics (0.3,0.35, 0.4, 0.45)
set xlabel "t (sec)"
plot "res.dat" using 1:($4-$5) title "L1 residual" with lines lc rgb "black" lw 3





