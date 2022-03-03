set terminal pdf enhanced size 8in,18in font "Times,30"
set xlabel "Electronic Basis Set"
set ylabel "Protonic Basis Set"
set size ratio 1
set palette defined (0 "white", 1 "slateblue1")
set cbrange[0:1.0]
set cblabel "Normalized Computational Cost"
set output "small_system_2D.pdf"
set multiplot
set size 1.0,0.333
set origin 0.0,0.667
set label at 20,245 "(A)" front font "Times BOLD, 30"
set xtics ("cc-pvdz" 0, "cc-pvtz" 1, "cc-pvqz" 2)
set ytics ("pb4d" 0, "pb4f1" 1, "pb4f2" 2)
set title "COH2"
plot "coh2.map" u ($1):($2):3 w image t ""
set size 1.0,0.333
set origin 0.0,0.333
unset label
set label at 20,245 "(B)" front font "Times BOLD, 30"
set xtics ("cc-pvdz" 0, "cc-pvtz" 1, "cc-pvqz" 2)
set ytics ("pb4d" 0, "pb4f1" 1, "pb4f2" 2)
set title "HCN"
plot "hcn.map" u ($1):($2):3 w image t ""
set size 1.0,0.333
set origin 0.0,0.0
unset label
set label at 20,245 "(C)" front font "Times BOLD, 30"
set xtics ("cc-pvdz" 0, "cc-pvtz" 1, "cc-pvqz" 2)
set ytics ("pb4d" 0, "pb4f1" 1, "pb4f2" 2)
set title "FHF-"
plot "fhf.map" u ($1):($2):3 w image t ""
unset multiplot
