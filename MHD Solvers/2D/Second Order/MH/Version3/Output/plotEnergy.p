#Energy Plot

set terminal png
set output 'energy.png'


set key font ",10"

set autoscale
  set title "Total Energy" 
  set multiplot layout 1,2 
  #plot 1: fluid density
  set xlabel "t[code units]"
  set ylabel "Energy[code units]"
  plot "energy.txt" using 1:2 with lines lc rgb "blue" title "Kinetic Energy" ,\
  "energy.txt" using 1:4 with lines lc rgb "red" title "Magnetic Energy" 
  unset xlabel
  unset ylabel
  set xlabel "t[code units]"
  set ylabel "Energy[code units]"
  plot "energy.txt" using 1:3 with lines lc rgb "green" title "Internal Energy"
  unset xlabel
  unset ylabel
  unset title

  unset multiplot
unset output
