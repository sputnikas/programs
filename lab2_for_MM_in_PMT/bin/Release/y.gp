set grid x y
set key top right outside
set title sprintf("%s y(t)", filename) font ",10"
set xlabel "t"
set ylabel "y"
plot filename using 1:4 with lines title 'runge2a',\
	 filename using 1:6 with lines title 'runge2b',\
	 filename using 1:8 with lines title 'runge2c',\
	 filename using 1:10 with lines title 'runge3',\
	 filename using 1:12 with lines title 'runge4a',\
	 filename using 1:14 with lines title 'runge4b',\
	 filename using 1:16 with lines title 'adams2',\
	 filename using 1:18 with lines title 'adams3',\
	 filename using 1:20  with lines title 'adams4'
	 #filename using 1:2 with lines title 'euler'
pause mouse keypress