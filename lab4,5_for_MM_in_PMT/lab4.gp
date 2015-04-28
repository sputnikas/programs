set term png size 800,800 enhanced font 'Verdana,10'
set output 'f/lab4-14-f.png'
set multiplot layout 2,2 rowsfirst
unset key 
set grid 
plot 'lab4.dat' using 1:2 w l lw 2, '' using 1:8 w l lw 2
plot 'lab4.dat' using 3:4 w impulses lw 2, '' using 3:5 w impulses lw 2 
plot 'lab4.dat' using 3:6 w impulses lw 2
plot 'lab4.dat' using 3:7 w impulses lw 2
