set dummy u,v
set samples 51, 51
set isosamples 41, 41
set xrange [ -1.00000 : 1.00000 ] noreverse nowriteback
set yrange [ -1.00000 : 1.00000 ] noreverse nowriteback
set pm3d
set hidden3d
splot [x=-3:3] [y=-3:3] sin(x) * cos(y)