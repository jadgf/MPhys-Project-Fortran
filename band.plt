ef=  5.99799786
set xtics ( \
"L"   -0.837034, \
"A"    0.000000, \
"H"    0.966523 )
set xrange [ -0.15 : 0.15]
set terminal pdfcairo enhanced font "DejaVu"  transparent fontscale 1 size 5.00in, 7.50in
set output "band.pdf"
set encoding iso_8859_1
set size ratio 0 1.0,1.0
set ylabel "E-E_{CBM} (eV)"
set yrange [-2 : 2 ]
unset key
set ytics 1.0 scale 1 nomirror out
set mytics 2
set parametric
set trange [-10:10]
set multiplot
plot "band.dat" u 1:($2-ef):(column(3)*1.5) with points pt 7 ps variable lc rgb "blue"
plot "band.dat" u 1:($2-ef):(column(4)*1.5) with points pt 7 ps variable lc rgb "red"
plot "band.dat" u 1:($2-ef):(column(5)*1.5) with points pt 7 ps variable lc rgb "green"
plot "band.dat" u 1:($2-ef) with l lt 1 lw 1.5 lc rgb "yellow",\
    0.000000,t with l lt 2  lc -1,\
   t,0 with l lt 2  lc -1
unset multiplot
