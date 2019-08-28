set term pngcairo enhanced truecolor font "Helvetica,14" size 850,400
set output "Gauss_check.png"
set multiplot
set size 0.5,1.12
set ylabel "Im"
set xlabel "Re"
set xrange [-3.5:3.5]
set yrange [-3.5:3.5]
# colorbrewer palette colors
set palette defined (0 '#f7fbff', 1 '#deebf7', 2 '#c6dbef', 3 '#9ecae1', 4 '#6baed6', 5 '#4292c6', 6 '#2171b5',7 '#084594')
set pm3d map impl
set contour
set origin -0.02,-0.02
do for [i=1:18] {set style line i lc rgb "black" lw 1}
set style increment userstyle
set cntrlabel onecolor
#set view map
splot 'freq_kde.dat' u 1:2:3 with pm3d nocontour notitle, \
    'freq_kde.dat' u 1:2:3 w l lc rgb "black" nosurface notitle
reset
set size 0.5,0.95
set origin 0.52,0.04
set ylabel "PDF"
set xlabel "Deviation"
set xrange [-6:6]
set logscale y
set format y '10^{%T}'
set yrange [1e-3:]
plot "hist_freq.dat" using 1:2 notitle with lines lc rgb "blue", \
"hist_freq.dat" using 1:3 notitle with lines lc rgb "black", \
"hist_freq.dat" using 1:($3*(1.0+3.0/sqrt($4))) notitle with lines \
dashtype 2 lc rgb "black", "hist_freq.dat" using 1:($3*(1.0-3.0/sqrt($4)))\
 notitle with lines dashtype 2 lc rgb "black"




