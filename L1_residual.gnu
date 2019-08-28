set term png enhanced truecolor crop font Helvetica 24  size 1600,800
set output "L1residual.png"
set pm3d map corners2color c1
set ylabel "f (Hz)"
set xlabel "t (s)"
set logscale y
#set xrange [15.5:16.5]
set xrange [-1.5:1.5]
set yrange [32:512]
set ytics (8,16,32,64,128,248,512,1024,2048)
set cbrange [0:9]
set cbtics (0,3,6,9)
# colorbrewer palette colors
set palette defined (0 '#fff7ec', 1'#fee8c8', 2 '#fdd49e', 3 '#fdbb84', 4 '#fc8d59', 5 '#ef6548', 6 '#d7301f',7 '#990000')
splot "residualL1.dat" using 1:2:3 notitle




