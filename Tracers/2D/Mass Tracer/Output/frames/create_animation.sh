#!/bin/sh
convert -delay 25 -loop 0 t={500..80000..500}.png output_fluid_slow_t=80000.gif
