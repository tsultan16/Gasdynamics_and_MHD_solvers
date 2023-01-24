set terminal png
set output 'L1.png'

set key font ",10"

set title "Order of Accuracy"
set xlabel "log(dx)"
set ylabel "log(L1 error)"
plot "error_weno.txt" using 1:2 with points pointtype 7 lc rgb "blue" title "WENO"
#"error_vl.txt" using 1:2 with points pointtype 5 lc rgb "red" title "Van Leer",\
#"error_vl.txt" using 1:1 with lines lc rgb "red" title "dx"
unset output
