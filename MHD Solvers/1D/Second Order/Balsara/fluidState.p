#Fluid State Plot
#plot from file: output.txt

set terminal gif animate delay 10 size 1280, 680
set output 'FluidState_Roe.gif'

TOP=0.90
DY = 0.23
tSteps=300
nx=200

set key font ",10"

set autoscale
do for [i=0:tSteps-1] {
  set title "Time Step = #".(i+1)
  set multiplot layout 2,3 
  #plot 1: fluid density
  #set label 1 'Fluid Density' at graph 0.92,0.9 font ',8'
  #set xlabel "x"
  set ylabel "Density"
  plot "output_fluid.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:2 with linespoint pointtype 7 lc rgb "blue" notitle
  unset title
  #unset xlabel
  unset ylabel
  #plot 2: fluid velocity
  #set label 1 'Fluid Velocity' at graph 0.92,0.9 font ',8'
  #set xlabel "x"
  set ylabel "Vx"
  plot "output_fluid.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:3 with linespoint pointtype 7 lc rgb "blue" notitle 
  #unset xlabel
  unset ylabel
 #plot 3: fluid velocity
  #set label 1 'Fluid Velocity' at graph 0.92,0.9 font ',8'
  #set xlabel "x"
  set ylabel "Vy"
  plot "output_fluid.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:4 with linespoint pointtype 7 lc rgb "blue" notitle 
  #unset xlabel
  unset ylabel
 #plot 4: fluid velocity
  #set label 1 'Fluid Velocity' at graph 0.92,0.9 font ',8'
  #set xlabel "x"
  set ylabel "Vz"
  plot "output_fluid.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:5 with linespoint pointtype 7 lc rgb "blue" notitle 
  #unset xlabel
  unset ylabel
  #plot 5: fluid pressure
  #set label 1 'Fluid Pressure' at graph 0.92,0.9 font ',8'
  set xlabel "x"
  set ylabel "Pressure"
  plot "output_fluid.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:6 with linespoint pointtype 7 lc rgb "blue" notitle  
  unset xlabel
  unset ylabel
  unset multiplot
}

unset output
