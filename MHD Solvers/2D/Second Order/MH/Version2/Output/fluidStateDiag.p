#Fluid State Plot

set terminal gif animate delay 10 size 1280, 680
set output 'FluidState_x=y.gif'

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
  set xlabel "x"
  set ylabel "Density"
  set tmargin at screen TOP-0*DY
  set bmargin at screen TOP-1*DY
  plot "diagcut_fluid.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:2 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset title
  unset xlabel
  unset ylabel
  #plot 2: fluid velocity
  set tmargin at screen TOP-1*DY-0.05
  set bmargin at screen TOP-2*DY-0.05
  set xlabel "x"
  set ylabel "Vx"
  plot "diagcut_fluid.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:3 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  #plot 3: fluid velocity
  set tmargin at screen TOP-1*DY-0.05
  set bmargin at screen TOP-2*DY-0.05
  set xlabel "x"
  set ylabel "Vy"
  plot "diagcut_fluid.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:4 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  #plot 4: fluid velocity
  set tmargin at screen TOP-1*DY-0.05
  set bmargin at screen TOP-2*DY-0.05
  set xlabel "x"
  set ylabel "Vz"
  plot "diagcut_fluid.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:5 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  #plot 5: fluid pressure
  set tmargin at screen TOP-2*DY-0.1
  set bmargin at screen TOP-3*DY-0.1
  set xlabel "x"
  set ylabel "Pressure"
  plot "diagcut_fluid.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:6 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  unset multiplot
}
unset output
