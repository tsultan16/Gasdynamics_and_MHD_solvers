#set terminal png

set terminal gif animate delay 10 size 1280, 680
set output 'weno_advection.gif'

set key font ",10"

nx=200
nt=200

do for [i=0:nt-1] {
  set title "Time Step = #".(i+1)
  #set title "Piecewise Polynomial Reconstruction"
  set xlabel "x"
  set ylabel "v(x)"
  set yrange [-1.1:1.1]
  set pointsize 0.6
  plot "output_adv.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:2 with linespoints pointtype 6 lc rgb "red" title "WENO5"#,\
# "output_adv.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:3 with lines lc rgb "green" title "Exact"

}


unset output
