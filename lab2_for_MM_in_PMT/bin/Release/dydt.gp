set grid x y
set key top right outside
set title sprintf("%s dy/dt(t)", filename) font ",10"
set xlabel "t"
set ylabel "dy/dt"
plot filename using 1:5  with lines title 'runge2a',\
	 filename using 1:7  with lines title 'runge2b',\
	 filename using 1:9  with lines title 'runge2c',\
	 filename using 1:11 with lines title 'runge3',\
	 filename using 1:13 with lines title 'runge4a',\
	 filename using 1:15 with lines title 'runge4b',\
	 filename using 1:17 with lines title 'adams2',\
	 filename using 1:19 with lines title 'adams3',\
	 filename using 1:21 with lines title 'adams4'
     #filename using 1:3  with lines title 'euler'   
pause mouse keypress	 