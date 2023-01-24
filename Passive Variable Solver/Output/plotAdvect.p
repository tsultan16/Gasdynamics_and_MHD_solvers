#Fluid State Plot
#plot from file: output.txt

set terminal gif animate delay 40 size 1280, 680
set output 'FluidState_y=0.5.gif'

TOP=0.90
DY = 0.23

nx=200
tSteps=100
tSkip=1

set key font ",10"

set autoscale
do for [i=0:tSteps-1] {
    j=i*tSkip
  set title "Time Step = #".(i+1)
  #set multiplot layout 3,1 rowsfirst
  #plot 1: fluid density
  #set label 1 'Fluid Density' at graph 0.92,0.9 font ',8'
  #set xlabel "x"
  set ylabel "Density"
  set yrange [0:55]
  #set tmargin at screen TOP-0*DY
  #set bmargin at screen TOP-1*DY
  plot "horcut_fluid.txt" every ::j*nx+1::nx+(j*nx)-1 using 1:2 with lines lc rgb "blue" notitle  ,\
  "horcut_fluid.txt" every ::j*nx+1::nx+(j*nx)-1 using 1:(0.044*$4) with points pointtype 7 lc rgb "red" 
  unset title
  unset xlabel
  unset ylabel
  unset yrange
  #plot 2: fluid velocity
  #set label 1 'Fluid Velocity' at graph 0.92,0.9 font ',8'
  #set tmargin at screen TOP-1*DY-0.05
  #set bmargin at screen TOP-2*DY-0.05
  #set xlabel "x"
  #set ylabel "Tracer Density"
  #set yrange [-5:305]
  #plot "horcut_fluid.txt" every ::j*nx+1::nx+(j*nx)-1 using 1:4 with linespoint pointtype 7 lc rgb "blue" notitle 
  #unset xlabel
  #unset ylabel
  #unset yrange

  #unset multiplot
}

unset output
