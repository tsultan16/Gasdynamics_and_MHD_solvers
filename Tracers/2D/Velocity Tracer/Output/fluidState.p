set term png  

set pm3d map
set palette color



set output 'FluidState_tf.png'


set dgrid3d 150,150

TOP=0.90
DY = 0.23

set key font ",10"

set autoscale

  j=5000
  
  filename="t=".j.".txt"
  set title "Time Step = #".j
  set multiplot layout 1,2 rowsfirst  
  #plot 1: fluid density
  set label 1 'Fluid Density' at graph 0.92,0.9 font ',8'
  set xlabel "x"
  set ylabel "y"
  #set tmargin at screen TOP-0*DY
  #set bmargin at screen TOP-1*DY
  splot filename using 1:2:3 
  unset title
  #unset xlabel
  unset ylabel
  #plot 2: fluid velocity
  set label 2 'Tracer' at graph 0.92,0.9 font ',8'
  #set tmargin at screen TOP-1*DY-0.05
  #set bmargin at screen TOP-2*DY-0.05
  set xlabel "x"
  set ylabel "y" 
  splot filename using 1:2:5
  unset xlabel
  unset ylabel
  unset multiplot

unset output
