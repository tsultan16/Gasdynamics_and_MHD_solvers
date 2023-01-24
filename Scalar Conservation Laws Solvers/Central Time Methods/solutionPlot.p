#Solution Plot
#plot from file: outputCTFS.txt,outputCTBS.txt,outputCTCS.txt

set terminal gif animate delay 70 size 1280, 680
set output 'u(x,t).gif'

TOP=0.90
DY = 0.23
tSteps=20
nx=128

set key font ",10"

set autoscale
do for [i=0:tSteps-1] {
  set title "Time Step = #".(i+1)
  set multiplot layout 3,1 rowsfirst
  #plot 1: CTFS
  set label 1 'CTFS' at graph 0.92,0.9 font ',8'
  #set xlabel "x"
  set ylabel "u"
  set tmargin at screen TOP-0*DY
  set bmargin at screen TOP-1*DY
  plot "outputCTFS.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:2 with linespoints notitle
  unset title
  #unset xlabel
  unset ylabel
  unset yrange
  #plot 2: CTBS
  set label 1 'CTBS' at graph 0.92,0.9 font ',8'
  set tmargin at screen TOP-1*DY-0.05
  set bmargin at screen TOP-2*DY-0.05
  #set xlabel "x"
  set ylabel "u"
  plot "outputCTBS.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:2 with linespoints notitle
  #unset xlabel
  unset ylabel
  unset yrange
  #plot 3: CTCS
  set label 1 'CTCS' at graph 0.92,0.9 font ',8'
  set tmargin at screen TOP-2*DY-0.1
  set bmargin at screen TOP-3*DY-0.1
  set xlabel "x"
  set ylabel "u"
  plot "outputCTCS.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:2 with linespoints notitle
  unset xlabel
  unset ylabel
  unset yrange
  unset multiplot
}

unset output
