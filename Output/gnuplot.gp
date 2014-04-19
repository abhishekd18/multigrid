set term post enh color solid
set output "./Output/Residual.ps"
set title 'Residual and Error'
set xlabel 'Iterations'
set ylabel 'Relative Residual'
set y2label 'Error'
set y2tics
set log y
set log y2
set grid x y
plot "./Output/ConvergenceHistory.out" u 1:4 w l axes x1y2 lw 1 lt 1 lc rgb "green" title "Error in inf norm","./Output/ConvergenceHistory.out" u 1:3 w l axes x1y1 lw 1 lt 1 lc rgb "red" title "Relative Residual"
reset
set term post enh color solid
set output "./Output/U.ps"
set title 'Grid function U'
set xlabel 'x'
set ylabel 'y'
set zlabel 'U(x,y)'
set zrange [-1.2:1.2]
set view 60,45
splot "./Output/U.out" w l lc rgb "green" title "U(x,y)"
