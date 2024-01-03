ef=  5.56322474
set xrange [ -0.12 : 0.12]
set terminal pdfcairo enhanced font "DejaVu"  transparent fontscale 1 size 5.00in, 5.00in
set output "band.pdf"
set border
unset xtics
unset ytics
set encoding iso_8859_1
set size ratio 0 1.0,1.0
set yrange [-0.4: 0.4 ]
unset key
set mytics 2
set parametric
set trange [-10:10]
set multiplot
plot "band.dat" every 4 u 1:($2-ef):(column(3)*4) with points pt 7 ps variable lc rgb "royalblue"
plot "band.dat" every 4 u 1:($2-ef):(column(4)*4) with points pt 7 ps variable lc rgb "light-red"
plot "band.dat" every 4 u 1:($2-ef):(column(5)*4) with points pt 7 ps variable lc rgb "forest-green"
plot "band.dat" u 1:($2-ef) with l lt 1 lw 3.5 lc rgb "yellow"
unset multiplot

