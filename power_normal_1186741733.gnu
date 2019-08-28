set term png enhanced truecolor crop font 'Helvetica,20'  size 1600,800
set output "power_normal_1186741733.png"
unset xtics
set ylabel "s"
set multiplot
set size 0.825, 0.39
set yrange [-4:4]
set xrange [-120:120]
set ytics (-3,-2,-1,0,1,2,3)
set origin 0.046, 0.56
plot "PowerTime_L1_1186741733_12.dat" using 1:2 notitle with lines, "PowerTime_L1_1186741733_13.dat" using 1:2 notitle with lines, "PowerTime_L1_1186741733_14.dat" using 1:2 notitle with lines, "PowerTime_L1_1186741733_15.dat" using 1:2 notitle with lines
#, "PowerTime_L1_1186741733_12.dat" using 1:(-3) notitle with lines lw 1 lc rgb "black" lt 2, "PowerTime_L1_1186741733_12.dat" using 1:(3) notitle with lines lw 1 lc rgb "black" dashtype 1
set size 0.95, 0.7
set origin 0.0, 0.0
set xtics
set pm3d map corners2color c1
set ylabel "f (Hz)"
set xlabel "t (secs)"
set ytics (100,200,300,400,500)
set yrange [20:510]
set cbrange [-3:3]
#unset cbtics
# colorbrewer palette colors
set palette defined (0 '#b2182b', 1 '#ef8a62', 2 '#fddbc7', 3 '#f7f7f7', 4 '#d1e5f0', 5 '#67a9cf', 6 '#2166ac')
splot "Spectogram_13_1186741733.dat" using 1:2:3 notitle




