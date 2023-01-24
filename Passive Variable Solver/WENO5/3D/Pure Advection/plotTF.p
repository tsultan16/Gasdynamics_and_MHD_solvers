set terminal png
set output 'eno-weno.png'

set key font ",10"

nt=188
nx=300

i=nt-1

  set title "Time Step = #".(i+1)
  #set title "Piecewise Polynomial Reconstruction"
  set xlabel "x"
  set ylabel "v(x)"
  set pointsize 0.6
  plot "output_adv.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:2 with points pointtype 6 lc rgb "red" title "WENO5",\
 "output_adv.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:3 with lines lc rgb "blue" title "Exact"

unset output
