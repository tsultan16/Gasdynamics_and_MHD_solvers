#Fluid State Plot

set terminal png 
set output 'mass tracers.png'


nx=200

tSteps=100
tSkip=1

a=0.01847

set key font ",10"

set autoscale
    
  
  set multiplot layout 2,1 rowsfirst
  #plot 1: fluid density
  #set label 1 'Fluid Density' at graph 0.92,0.9 font ',8'
  #set xlabel "x"
  set ylabel "Fluid Density[code units]"
  set yrange [0:55]

  j=1
  plot "horcut_fluid.txt" every ::j*nx+1::nx+(j*nx)-1 using 1:2 with lines lc rgb "red" title "t=1"  ,\
"horcut_fluid.txt" every ::(j*200)*nx+1::nx+((j*200)*nx)-1 using 1:2 with lines lc rgb "blue" title "t=200"  ,\
"horcut_fluid.txt" every ::(j*400)*nx+1::nx+((j*400)*nx)-1 using 1:2 with lines lc rgb "green" title "t=400"  
  #unset title
  #unset xlabel
  unset ylabel
  unset yrange
 
  set xlabel "x[code units]"
  set ylabel "Tracer Density[code units]"
  set yrange [0:55]
  plot "horcut_fluid.txt" every 2::j*nx+1::nx+(j*nx)-1 using 1:(a*$4) with lines lc rgb "red" title "t=1" ,\
  "horcut_fluid.txt" every 2::(j*200)*nx+1::nx+((j*200)*nx)-1 using 1:(a*$4) with lines lc rgb "blue" title "t=200" ,\
  "horcut_fluid.txt" every 2::(j*400)*nx+1::nx+((j*400)*nx)-1 using 1:(a*$4) with lines lc rgb "green" title "t=400"
  unset xlabel
  unset ylabel
  unset yrange

  unset multiplot
unset output
