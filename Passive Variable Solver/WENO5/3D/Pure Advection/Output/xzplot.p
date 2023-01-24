set pm3d map
set palette color

set terminal gif animate delay 5 size 1280, 680
set output 'xz_slice.gif'


set dgrid3d 50,50

tSteps=100
tSkip=1

set key font ",10"

set autoscale

do for [i=0:tSteps-1] { 
  j=i*tSkip
  print "Time Steps Completed =".i
  filename="xz_t=".j.".txt"
  set title "Time Step = #".j
  set multiplot layout 1,2 rowsfirst  
  #plot 1: WENO 

  set cbrange[-1.1:1.1]

  #set label 1 'WENO' at graph 0.92,0.9 font ',8'
  set xlabel "x"
  set ylabel "y"
  splot filename using 1:2:3 
  unset title
  unset xlabel
  unset ylabel
  #plot 2:
  #set label 2 'Exact' at graph 0.92,0.9 font ',8'
  set xlabel "x"
  set ylabel "y" 
  splot filename using 1:2:4
  unset xlabel
  unset ylabel
  unset multiplot

}
unset output
