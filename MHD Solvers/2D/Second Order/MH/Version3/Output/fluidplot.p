set pm3d map
set palette color

set terminal gif animate delay 100 size 1280, 680
set output 'fluidStateMHD.gif'


set dgrid3d 500,500

TOP=0.90
DY = 0.23
tSteps=120
tSkip=1000

set key font ",10"

set autoscale

do for [i=0:tSteps-1] { 
  j=i*tSkip
  print "Time Steps Completed =".i
  filename="KH_500x500_120000 timesteps/t=".j.".txt"
  set title "Time Step = #".j
  set multiplot layout 1,2 rowsfirst  
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
  #plot 2: fluid velocity
  #set label 2 'Pressure' at graph 0.92,0.9 font ',8'
  set tmargin at screen TOP-1*DY-0.05
  set bmargin at screen TOP-2*DY-0.05
  set xlabel "x"
  set ylabel "y" 
  splot filename using 1:2:7
  unset xlabel
  unset ylabel
  unset multiplot

}
unset output
