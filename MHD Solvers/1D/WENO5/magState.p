#Magnetic Field State Plot
#plot from file: output.txt

set terminal gif animate delay 10 size 1280, 680
set output 'MagneticField_Roe.gif'

TOP=0.90
DY = 0.23
nx=200
tSteps=300


set key font ",10"

set autoscale
do for [i=0:tSteps-1] {
  set title "Time Step = #".(i+1)
  set multiplot layout 1,2 
  #plot 1: magnetic field y
  set label 1 'Magnetic Field' at graph 0.92,0.9 font ',8'
  set xlabel "x"
  set ylabel "By"
  plot "output_mag.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:2 with linespoint pointtype 7 lc rgb "blue" notitle  
  #unset xlabel
  unset ylabel
  unset yrange
  #plot 2: magnetic field z
  set label 1 'Magnetic Field' at graph 0.92,0.9 font ',8'
  set xlabel "x"
  set ylabel "Bz"
  #set yrange [-500:100500]
  plot "output_mag.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:3 with linespoint pointtype 7 lc rgb "blue" notitle  
 
  unset xlabel
  unset ylabel
  unset yrange

  unset multiplot
}

unset output
