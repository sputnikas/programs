set grid x y
set key top right outside
set title sprintf("%s dy/dt(y)", filename) font ",10"
set xlabel "y"
set ylabel "dy/dt"
plot filename using 4:5   with lines title 'runge2a',\
	 filename using 6:7   with lines title 'runge2b',\
	 filename using 8:9   with lines title 'runge2c',\
	 filename using 10:11 with lines title 'runge3',\
	 filename using 12:13 with lines title 'runge4a',\
	 filename using 14:15 with lines title 'runge4b',\
	 filename using 16:17 with lines title 'adams2',\
	 filename using 18:19 with lines title 'adams3',\
	 filename using 20:21	with lines title 'adams4'
	 #filename using 2:3  with lines title 'euler'
pause mouse keypress                                                