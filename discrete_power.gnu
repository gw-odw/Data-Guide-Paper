set term png enhanced truecolor crop font Helvetica 20  size 1600,800
set output "Pow.png"
set pm3d map corners2color c1
set ylabel "f (Hz)"
set xlabel "t (secs)"
#set logscale y
set xrange [-120:120]
set yrange [16:250]
set cbrange [0:9]
#unset cbtics
# colorbrewer palette colors
set palette defined (0 '#fff7ec', 1'#fee8c8', 2 '#fdd49e', 3 '#fdbb84', 4 '#fc8d59', 5 '#ef6548', 6 '#d7301f',7 '#990000')
splot "Binary.dat" using 1:2:($3*$3) notitle




