set term png enhanced truecolor crop font Helvetica 25  size 1600,1200
set output "templatesNR.png"
set xrange [0.2:0.5]
unset xtics
set ylabel "h(t)"
unset xtics
set multiplot
set size 0.937, 0.35
set origin 0.063, 0.62
set yrange [-1.5:1.5]
set ytics (-1, 0, 1)
set label "Reference Template" at graph 0.04,0.9
plot "ref_template.dat" using 1:($2*1e21) notitle with lines lc rgb "black" lw 2
set size 1.0, 0.33
set origin 0.0, 0.345
set yrange [-1.2e-21:1.2e-21]
set ytics (-1e-21, -5e-22, 0, 5e-22,1e-21)
unset label
set label "Hanford" at graph 0.04,0.9
plot "templates.dat" using 1:2 notitle with lines lc rgb "blue" lw 2
set size 1.0, 0.4
set origin 0.0, 0.0
unset label
set xtics
set xlabel "t (s)"
set label "Livingston" at graph 0.04,0.9
plot "templates.dat" using 1:3 notitle with lines lc rgb "blue" lw 2