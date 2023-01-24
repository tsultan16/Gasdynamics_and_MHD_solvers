set pm3d map
set palette color
#set palette rgb 33,13,10

set terminal png


set dgrid3d 400,400

TOP=0.90
DY = 0.23
tSteps=20
tSkip=200

set key font ",5"

set autoscale

print "Total frames =".tSteps

do for [i=1:tSteps] { 
  j=i*tSkip

  filename="t=".j.".txt"
  filename2="fluid_t_=".j.".png"
  filename3="velocity_tracer_t_=".j.".png"


  set output filename2
  set title "Time Step = #".i
  
  set xlabel "x"
  set ylabel "y"
  splot filename using 1:2:3  
  unset title
  unset xlabel
  unset ylabel
  unset output

  set output filename3
  set xlabel "x"
  set ylabel "y" 
  splot filename using 1:2:5
  unset xlabel
  unset ylabel
  unset output
 
  print "Frames Completed =".i
  
}

