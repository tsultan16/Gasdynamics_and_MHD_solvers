set pm3d map
set palette color

set terminal gif animate delay 100 size 1280, 680
set output 'fluidStateMHD.gif'


set dgrid3d 150,150

TOP=0.90
DY = 0.23
tSteps=50
tSkip=1000

set key font ",10"

#set autoscale

do for [i=0:tSteps-1] { 
  j=i*tSkip
  print "Fluid Pressure, Time Steps Completed =".i
  filename="t=".j.".txt"
  set title "Time Step = #".j
  #plot 1: fluid density
  set xlabel "x"
  set ylabel "y"
  #set tmargin at screen TOP-0*DY
  #set bmargin at screen TOP-1*DY
  splot filename using 1:2:6
  unset title
  unset xlabel
  unset ylabel
  

}
unset output
