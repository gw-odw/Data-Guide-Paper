set term png enhanced truecolor crop font Helvetica 20  size 1600,800
set output "Amp.png"
set pm3d map corners2color c1
set ylabel "f (Hz)"
set xlabel "t (secs)"
#set logscale y
set xrange [-120:120]
set yrange [16:250]
set cbrange [-3:3]
#unset cbtics
# colorbrewer palette colors
set palette defined (0 '#b2182b', 1 '#ef8a62', 2 '#fddbc7', 3 '#f7f7f7', 4 '#d1e5f0', 5 '#67a9cf', 6 '#2166ac')
splot "Binary.dat" using 1:2:3 notitle




