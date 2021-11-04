set xrange [0:100]
set yrange [0:500]

set key above
set xlabel "t"
set ylabel ""

#Metoda Trapezow z metoda Picarda 
plot "uP.txt" lc 7 lw 2 with lines title "u(t)",\
"zP.txt" lc 18 lw 2 with lines title "z(t)"

#Metoda Trapezow z iteracja Newtona
plot "uN.txt" lc 7 lw 2 with lines title "u(t)",\
"zN.txt" lc 18 lw 2 with lines title "z(t)"

#Metoda niejawna RK2
plot "uRK2.txt" lc 20 lw 2 with lines title "u(t)",\
"zRK2.txt" lc 9 lw 2 with lines title "z(t)"