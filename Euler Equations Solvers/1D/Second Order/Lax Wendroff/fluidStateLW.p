#Fluid State Plot
#plot from file: ouput_euler_godunov.txt

set terminal gif animate delay 50 size 1280, 680
set output 'LWFluidState.gif'

TOP=0.90
DY = 0.23
tSteps=100
nx=256

set key font ",10"

set autoscale
do for [i=0:tSteps-1] {
  set title "Time Step = #".(i+1)
  set multiplot layout 3,1 rowsfirst
  #plot 1: fluid density
  set label 1 'Fluid Density' at graph 0.92,0.9 font ',8'
  #set xlabel "x(m)"
  set ylabel "Density(kg/m^3)"
  #set yrange [-0.5:1.05]
  set tmargin at screen TOP-0*DY
  set bmargin at screen TOP-1*DY
  plot "output_euler_LW.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:2 with linespoint pointtype 6 lc rgb "blue"
  unset title
  #unset xlabel
  unset ylabel
  unset yrange
  #plot 2: fluid velocity
  set label 1 'Fluid Velocity' at graph 0.92,0.9 font ',8'
  set tmargin at screen TOP-1*DY-0.05
  set bmargin at screen TOP-2*DY-0.05
  #set xlabel "x(m)"
  set ylabel "Velocity(m/s)"
  #set yrange [-5:305]
  plot "output_euler_LW.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:3 with linespoint pointtype 6 lc rgb "blue"
  #unset xlabel
  unset ylabel
  unset yrange
  #plot 3: fluid pressure
  set label 1 'Fluid Pressure' at graph 0.92,0.9 font ',8'
  set tmargin at screen TOP-2*DY-0.1
  set bmargin at screen TOP-3*DY-0.1
  set xlabel "x"
  set ylabel "Pressure(N/m^2)"
  #set yrange [-500:100500]
  plot "output_euler_LW.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:4 with linespoint pointtype 6 lc rgb "blue"
  unset xlabel
  unset ylabel
  unset yrange
  unset multiplot
}

unset output
