set term png enhanced truecolor crop font Helvetica 20  size 1600,800
set output "powerfluc_1165067917.png"
set xrange [-120:120]
unset xtics
set ylabel "Average Power"
set multiplot
set size 0.845, 0.4
set yrange [0.5:2.5]
set ytics (0.5,1,1.5,2,2.5)
set origin 0.0255, 0.57
plot "PowerTime_H1_1165067917_12.dat" using 1:4 notitle with lines
set size 0.95, 0.7
set origin 0.0, 0.0
set xtics
set pm3d map corners2color c1
set ylabel "f (Hz)"
set xlabel "t (secs)"
set ytics (50,100,150,200)
set yrange [16:250]
set cbrange [-3:3]
#unset cbtics
# colorbrewer palette colors
set palette defined (0 '#b2182b', 1 '#ef8a62', 2 '#fddbc7', 3 '#f7f7f7', 4 '#d1e5f0', 5 '#67a9cf', 6 '#2166ac')
splot "Spectogram_12_1165067917.dat" using 1:2:3 notitle




