
set terminal gif animate delay 20 size 1280, 680
set output 'magneticFieldLine.gif'


tSteps=40
tSkip=50

set object rectangle from screen 0,0 to screen 1,1 behind fillcolor rgb 'gray' fillstyle solid noborder
set key font ",10"

set autoscale

do for [i=1:tSteps] { 
  j=i*tSkip
  print "Time Steps Completed =".i
  filename="fieldLine1_t=".j.".txt"
  filename2="field2_t=".j.".txt"
  set title "Time Step = #".j
  set xlabel "x"
  set ylabel "y"
  set xrange [0:1]
  set yrange [0:1]
  plot filename using 1:2 with lines lc rgb "blue" notitle #,\
  #filename2 using 1:2 with lines lc rgb "red" notitle
  unset title
  unset xlabel
  unset ylabel
}
unset output
