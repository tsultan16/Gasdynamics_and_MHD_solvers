#set terminal png

set terminal gif animate delay 1 size 1280, 680
set output 'eno-weno.gif'

set key font ",10"

nx=200
nt=500

do for [i=0:nt-1] {
  set title "Time Step = #".(i+1)
  #set title "Piecewise Polynomial Reconstruction"
  set xlabel "x"
  set ylabel "v(x)"
  set pointsize 0.6
  plot "output.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:2 with linespoints pointtype 6 lc rgb "blue" title "WENO5" ,\
  "output.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:3 with lines lc rgb "red" title "Exact"
}


unset output
