set terminal png
set output 'eno-weno.png'

set key font ",10"


nx=50

  set xlabel "x"
  set ylabel "v(x)"
  set yrange [-1.1:1.1]
  set pointsize 0.6
  plot "output_adv.txt" using 1:2 with points pointtype 6 lc rgb "red" title "WENO5",\
 "output_adv.txt" using 1:3 with points pointtype 3 lc rgb "blue" title "Exact"

unset output
