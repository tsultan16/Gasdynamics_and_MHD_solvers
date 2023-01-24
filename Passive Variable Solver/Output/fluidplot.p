set pm3d map
set palette color
#set palette rgb 33,13,10

set terminal gif animate delay 60 size 1280, 680
set output 'fluidState.gif'


set dgrid3d 300,300

TOP=0.90
DY = 0.23
tSteps=100
tSkip=200

set key font ",5"

set autoscale

print "Total time steps= ".(tSteps*tSkip)

do for [i=0:tSteps-1] { 
  j=i*tSkip
  filename="t=".j.".txt"
  set title "Time Step = #".i
  set multiplot layout 1,2 rowsfirst  
  #plot 1: fluid density
  set label 1 'Fluid Density' at graph 0.92,0.9 font ',8'
  set xlabel "x"
  set ylabel "y"
  splot filename using 1:2:3
  unset title
  unset xlabel
  unset ylabel

  #plot 2:Tracer density
  set label 2 'Tracer Density' at graph 0.92,0.9 font ',8'
  set xlabel "x"
  set ylabel "y" 
  splot filename using 1:2:5
  unset xlabel
  unset ylabel

  unset multiplot

  print "% Complete =".((i+1)*100/tSteps)

}
unset output
