#Fluid State Plot

set terminal gif animate delay 10 size 1280, 680
set output 'MagneticFieldState_x=y.gif'

TOP=0.90
DY = 0.23
tSteps=200
nx=120

set key font ",10"

set autoscale
do for [i=0:tSteps-1] {
  set title "Time Step = #".(i+1)
  set multiplot layout 1,3 
  #plot 1: fluid density
  set xlabel "x"
  set ylabel "Bx"
  set tmargin at screen TOP-0*DY
  set bmargin at screen TOP-1*DY
  plot "diagcut_B.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:2 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset title
  unset xlabel
  unset ylabel
  #plot 2: fluid velocity
  set tmargin at screen TOP-1*DY-0.05
  set bmargin at screen TOP-2*DY-0.05
  set xlabel "x"
  set ylabel "By"
  plot "diagcut_B.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:3 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  #plot 3: fluid velocity
  set tmargin at screen TOP-1*DY-0.05
  set bmargin at screen TOP-2*DY-0.05
  set xlabel "x"
  set ylabel "Bz"
  plot "diagcut_B.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:4 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel

  unset multiplot
}

unset output
