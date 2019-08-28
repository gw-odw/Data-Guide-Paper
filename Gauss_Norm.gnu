set term pdf enhanced truecolor font "Helvetica,12"
set output "GaussNorm.pdf"
set size 1.0, 1.0
set multiplot
set xrange [-4:4]
set logscale y
set format "%g"
set yrange [1e-4:]
set size 0.54,0.55
set title "Hanford"
set ylabel "PDF"
unset xtics
set format y '10^{%T}'
set origin 0.0, 0.49
plot "hist_freq_H.dat" u 1:($3*(1.0+3.0/sqrt($4))):($3*(1.0-3.0/sqrt($4))) w filledcurves lc rgb '#ffcccc' notitle, "hist_freq_H.dat" using 1:2 notitle with histeps lc rgb "blue", "hist_freq_H.dat" using 1:3 notitle with lines lc rgb "red"
unset ylabel
unset title
set title "Livingston"
set size 0.51,0.55
set origin 0.51, 0.49
plot "hist_freq_L.dat" u 1:($3*(1.0+3.0/sqrt($4))):($3*(1.0-3.0/sqrt($4))) w filledcurves lc rgb '#ffcccc' notitle, "hist_freq_L.dat" using 1:2 notitle with histeps lc rgb "blue", "hist_freq_L.dat" using 1:3 notitle with lines lc rgb "red"
unset title
set origin 0.012, 0.0
set xtics
set size 0.527,0.55
set xlabel "Quantile"
set ylabel "Observed Quantile"
set format y "%g"
unset logscale y
set yrange [-4:4]
plot "PP.dat" using 1:2 notitle with points lc rgb "blue" pt 6 ps 0.3, "PPref.dat" using 1:2 notitle with lines lc rgb "black"
unset ylabel
set size 0.497,0.55
set origin 0.522, 0.0
plot "PP.dat" using 1:3 notitle with points  lc rgb "blue" pt 6 ps 0.3, "PPref.dat" using 1:2 notitle with lines lc rgb "black"


