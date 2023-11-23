#set title "Spin Projection"
set title "Orbital Angular Momentum Projection"
set terminal pdfcairo enhanced font "DejaVu" transparent fontscale 1 size 8.00in, 8.00in
set output "pdfs/kmesh.pdf"
set encoding iso_8859_1
set size ratio 0 1.0,1.0
 
set xlabel "k_x"
set ylabel "k_y"
set xrange [-0.1 : 0.1]
set yrange [ -0.1 : 0.1 ]
unset key
set ytics 0.05 scale 1 nomirror out
set mytics 2
set multiplot
 
#Connect two rings separately
R(x,y) = sqrt(x*x + y*y)
theta(x,y) = atan2(y,x)
gap = 0.05
set style data linespoints
set datafile missing NaN
 
#Palette for arrows
set palette defined ( 0 "blue", 0.5 "white", 1 "red" )
set cblabel "m_z"
set style arrow 1 head filled size screen 0.02,10,45 lt 1 lc palette
 

#Plot two rings
plot "kmesh.dat" u (R($1,$2) < gap ? $1 : NaN) : ($2) : (theta($1,$2)) pt 7 ps 0.01 lt 5 lw 7 smooth zsort, \
     "kmesh.dat" u (R($1,$2) > gap ? $1 : NaN) : ($2) : (theta($1,$2)) pt 7 ps 0.01 lt 5 lw 7 smooth zsort, \
     "kmesh.dat" u 1:2:3:4:5 with vectors arrowstyle 1  #Spin Projection
#"kmesh.dat" u 1:2:6:7:8 with vectors arrowstyle 1  #Angular Momentum Projection
unset multiplot

