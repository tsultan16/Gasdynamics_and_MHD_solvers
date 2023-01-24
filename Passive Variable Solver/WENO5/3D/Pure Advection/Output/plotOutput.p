#set terminal png

set terminal gif animate delay 10 size 1280, 680
set output 'diagslice.gif'

set key font ",10"

nx=150
nt=100

do for [i=0:nt-1] {
  set title "Time Step = #".(i+1)
  #set title "Piecewise Polynomial Reconstruction"
  set xlabel "x"
  set ylabel "v(x)"
  set yrange [-1.1:1.1]
  set pointsize 0.6
  plot "output_adv1.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:2 with points pointtype 6 lc rgb "red" title "WENO5",\
 "output_adv1.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:3 with lines lc rgb "green" title "Exact"

}


unset output
