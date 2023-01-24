set term png  

set pm3d map
set palette color



set output 'FluidState_t=0.png'


set dgrid3d 150,150

TOP=0.90
DY = 0.23

set key font ",10"

set autoscale

  j=0
  
  filename="t=".j.".txt"
  set title "Time Step = #".j
  set multiplot layout 2,1   
  #plot 1: Pressure
  set xlabel "x"
  set ylabel "y"
  splot filename using 1:2:3
  unset title
  #unset xlabel
  unset ylabel
  #plot 2: Vx
  set xlabel "x"
  set ylabel "y" 
  splot filename using 1:2:5
  unset xlabel
  unset ylabel
  unset multiplot

unset output
