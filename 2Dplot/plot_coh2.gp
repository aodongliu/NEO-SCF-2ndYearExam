set terminal pdf enhanced size 8in,18in font "Times,30"
set xlabel "Electronic Basis Set"
set ylabel "Protonic Basis Set"
set size ratio 1
set palette defined (0 "white", 1 "#145DA0")
set cbrange[0:1.0]
set cblabel "Computational Cost"
set output "coh2_2D.pdf"
set xtics ("cc-pvdz" 0, "cc-pvtz" 1, "cc-pvqz" 2)
set ytics ("pb4d" 0, "pb4f1" 1, "pb4f2" 2)
set size 1.0,1
set origin 0.0,0.0
plot "coh2.map" u ($1):($2):3 w image t ""
unset multiplot
