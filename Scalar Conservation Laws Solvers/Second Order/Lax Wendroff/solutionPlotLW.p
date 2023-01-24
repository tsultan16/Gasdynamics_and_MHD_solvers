#Solution Plot
#plot from file: output_harten.txt

set terminal gif animate delay 1 size 1280, 680
set output 'u(x,t)_LW.gif'

TOP=0.90
DY = 0.23
tSteps=100
nx=500

set key font ",10"

set autoscale
do for [i=0:tSteps-1] {
  set title "Time Step = #".(i+1)
  set xlabel "x"
  set ylabel "u"
  #set yrange[-0.2:1.2]
  
  plot "output_LW.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:2 with linespoint pointtype 5 
 
}

unset output
