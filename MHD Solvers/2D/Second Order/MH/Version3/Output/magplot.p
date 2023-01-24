
set terminal gif animate delay  40 size 1280, 680
set output 'magneticFieldVector.gif'


tSteps=120
tSkip=1000

set key font ",10"

set autoscale

do for [i=1:tSteps] { 
  j=i*tSkip
  print "Time Steps Completed =".i
  filename="KH_500x500_120000 timesteps/b_t=".j.".txt"
  set title "Time Step = #".j
  set xlabel "x"
  set ylabel "y"
  set xrange [0:1]
  set yrange [0:1]
  plot filename using 1:2:(0.01*$3):(0.01*$4) with vectors 
  unset title
  unset xlabel
  unset ylabel
}
unset output
