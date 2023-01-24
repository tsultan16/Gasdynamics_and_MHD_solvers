set pm3d map
set palette color

set terminal gif animate delay 20 size 1280, 680
set output '2DMHDFluidState.gif'


set dgrid3d 200,200

TOP=0.90
DY = 0.23
tSteps=50
tSkip=100

set key font ",10"

set autoscale

do for [i=0:tSteps-1] { 
  j=i*tSkip
  print "Time Steps Completed =".i
  filename="t=".j.".txt"
  set title "Time Step = #".i
  set multiplot layout 1,2  
  #plot 1: fluid density
  #set label 1 'Density' at graph 0.92,0.9 font ',8'
  set xlabel "x"
  set ylabel "y"
  #set tmargin at screen TOP-0*DY
  #set bmargin at screen TOP-1*DY
  splot filename using 1:2:3 
  unset title
  unset xlabel
  unset ylabel
  #plot 2: Pressure
  #set label 2 'By' at graph 0.92,0.9 font ',8'
  #set tmargin at screen TOP-1*DY-0.05
  #set bmargin at screen TOP-2*DY-0.05
  set xlabel "x"
  set ylabel "y" 
  splot filename using 1:2:9 
  unset xlabel
  unset ylabel
  unset multiplot

}
unset output
