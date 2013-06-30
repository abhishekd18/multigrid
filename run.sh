echo "set term post enh color solid
set output \"./Output/PRECOND_RES.ps\"
set title 'Residual Comparison'
set xlabel 'Iterations'
set ylabel 'Relative Residual'
set y2label 'Error ||x* - xm||'
set y2tics
set log y
set log y2
set grid x y
plot \"./Output/GMRES_jacobi_true_res\" u 1:2 w l axes x1y2 lt 1 lw 5 lc rgb \"green\" title \"Error for residual based on true absolute residual\",\
\"./Output/GMRES_jacobi_true_res\" u 1:3 w l axes x1y1 lt 1 lc rgb \"red\" title \"Relative Residual based on true absolute residual\",\
\"./Output/GMRES_jacobi_precond_res\" u 1:2 w l axes x1y2 lt 1 lc rgb \"magenta\" title \"Error for residual based on preconditioned absolute residual\",\
\"./Output/GMRES_jacobi_precond_res\" u 1:3 w l axes x1y1 lt 1 lc rgb \"blue\" title \"Relative Residual based on preconditioned absolute residual\"">gnuplot.gp
gnuplot gnuplot.gp
ps2pdf -dEPSCrop -dNOSAFER ./Output/PRECOND_RES.ps ./Output/PRECOND_RES.pdf
rm ./Output/*Error*
rm *.out
./bin/main orsirr_1.mtx GMRES_FULL NULL>run_gmres_full_no.out
echo "set term post enh color solid
set output \"./Output/GMRES_FULL_NO.ps\"
set title 'GMRES FULL No Preconditioning'
set xlabel 'Iterations'
set ylabel 'Residual'
set y2label 'Error ||x* - xm||'
set y2tics
set log y
set log y2
set grid x y
plot \"./Output/Error_GMRES.out\" u 1:2 w l axes x1y2 lt 1 lc rgb \"green\" title \"Error ||x* - xm||\",\
\"./Output/Error_GMRES.out\" u 1:3 w l axes x1y1 lt 1 lc rgb \"red\" title \"Relative Residual\",\
\"./Output/Error_GMRES.out\" u 1:4 w l axes x1y1 lt 1 lc rgb \"blue\" title \"Absolute Residual\"">gnuplot.gp
gnuplot gnuplot.gp
ps2pdf -dEPSCrop -dNOSAFER ./Output/GMRES_FULL_NO.ps ./Output/GMRES_FULL_NO.pdf
mv ./Output/Error_GMRES.out ./Output/Error_GMRES_FULL_NO.out

./bin/main orsirr_1.mtx GMRES_FULL JACOBI>run_gmres_full_jacobi.out
echo "set term post enh color solid
set output \"./Output/GMRES_FULL_JACOBI.ps\"
set title 'GMRES FULL JACOBI'
set xlabel 'Iterations'
set ylabel 'Residual'
set y2label 'Error ||x* - xm||'
set y2tics
set log y
set log y2
set grid x y
plot \"./Output/Error_GMRES.out\" u 1:2 w l axes x1y2 lt 1 lc rgb \"green\" title \"Error ||x* - xm||\",\
\"./Output/Error_GMRES.out\" u 1:3 w l axes x1y1 lt 1 lc rgb \"red\" title \"Relative Residual\",\
\"./Output/Error_GMRES.out\" u 1:4 w l axes x1y1 lt 1 lc rgb \"blue\" title \"Absolute Residual\"">gnuplot.gp
gnuplot gnuplot.gp
ps2pdf -dEPSCrop -dNOSAFER ./Output/GMRES_FULL_JACOBI.ps ./Output/GMRES_FULL_JACOBI.pdf
mv ./Output/Error_GMRES.out ./Output/Error_GMRES_FULL_JACOBI.out

./bin/main orsirr_1.mtx GMRES_FULL GAUSS_SEIDEL>run_gmres_full_gauss_seidel.out
echo "set term post enh color solid
set output \"./Output/GMRES_FULL_GAUSS_SEIDEL.ps\"
set title 'GMRES FULL GAUSS SEIDEL'
set xlabel 'Iterations'
set ylabel 'Residual'
set y2label 'Error ||x* - xm||'
set y2tics
set log y
set log y2
set grid x y
plot \"./Output/Error_GMRES.out\" u 1:2 w l axes x1y2 lt 1 lc rgb \"green\" title \"Error ||x* - xm||\",\
\"./Output/Error_GMRES.out\" u 1:3 w l axes x1y1 lt 1 lc rgb \"red\" title \"Relative Residual\",\
\"./Output/Error_GMRES.out\" u 1:4 w l axes x1y1 lt 1 lc rgb \"blue\" title \"Absolute Residual\"">gnuplot.gp
gnuplot gnuplot.gp
ps2pdf -dEPSCrop -dNOSAFER ./Output/GMRES_FULL_GAUSS_SEIDEL.ps ./Output/GMRES_FULL_GAUSS_SEIDEL.pdf
mv ./Output/Error_GMRES.out ./Output/Error_GMRES_FULL_GAUSS_SEIDEL.out

./bin/main orsirr_1.mtx GMRES_RESTARTED NULL>run_gmres_restarted_no.out
echo "set term post enh color solid
set output \"./Output/GMRES_RESTARTED_NO.ps\"
set title 'GMRES RESTARTED No Preconditioning'
set xlabel 'Iterations'
set ylabel 'Residual'
set y2label 'Error ||x* - xm||'
set y2tics
set log y
set log y2
set grid x y
plot \"./Output/Error_GMRES.out\" u 1:2 w l axes x1y2 lt 1 lc rgb \"green\" title \"Error ||x* - xm||\",\
\"./Output/Error_GMRES.out\" u 1:3 w l axes x1y1 lw 0.5 lt 1 lc rgb \"red\" title \"Relative Residual\",\
\"./Output/Error_GMRES.out\" u 1:4 w l axes x1y1 lw 0.5 lt 1 lc rgb \"blue\" title \"Absolute Residual\"">gnuplot.gp
gnuplot gnuplot.gp
ps2pdf -dEPSCrop -dNOSAFER ./Output/GMRES_RESTARTED_NO.ps ./Output/GMRES_RESTARTED_NO.pdf
mv ./Output/Error_GMRES.out ./Output/Error_GMRES_RESTARTED_NO_internal.out

./bin/main orsirr_1.mtx GMRES_RESTARTED JACOBI>run_gmres_restarted_jacobi.out
echo "set term post enh color solid
set output \"./Output/GMRES_RESTARTED_JACOBI.ps\"
set title 'GMRES RESTARTED JACOBI'
set xlabel 'Iterations'
set ylabel 'Residual'
set y2label 'Error ||x* - xm||'
set y2tics
set log y
set log y2
set grid x y
plot \"./Output/Error_GMRES.out\" u 1:2 w l axes x1y2 lt 1 lc rgb \"green\" title \"Error ||x* - xm||\",\
\"./Output/Error_GMRES.out\" u 1:3 w l axes x1y1 lw 0.5 lt 1 lc rgb \"red\" title \"Relative Residual\",\
\"./Output/Error_GMRES.out\" u 1:4 w l axes x1y1 lw 0.5 lt 1 lc rgb \"blue\" title \"Absolute Residual\"">gnuplot.gp
gnuplot gnuplot.gp
ps2pdf -dEPSCrop -dNOSAFER ./Output/GMRES_RESTARTED_JACOBI.ps ./Output/GMRES_RESTARTED_JACOBI.pdf
mv ./Output/Error_GMRES.out ./Output/Error_GMRES_RESTARTED_JACOBI_internal.out

./bin/main orsirr_1.mtx GMRES_RESTARTED GAUSS_SEIDEL>run_gmres_restarted_gauss_seidel.out
echo "set term post enh color solid
set output \"./Output/GMRES_RESTARTED_GAUSS_SEIDEL.ps\"
set title 'GMRES RESTARTED GAUSS SEIDEL'
set xlabel 'Iterations'
set ylabel 'Residual'
set y2label 'Error ||x* - xm||'
set y2tics
set log y
set log y2
set grid x y
plot \"./Output/Error_GMRES.out\" u 1:2 w l axes x1y2 lt 1 lc rgb \"green\" title \"Error ||x* - xm||\",\
\"./Output/Error_GMRES.out\" u 1:3 w l axes x1y1 lw 0.5 lt 1 lc rgb \"red\" title \"Relative Residual\",\
\"./Output/Error_GMRES.out\" u 1:4 w l axes x1y1 lw 0.5 lt 1 lc rgb \"blue\" title \"Absolute Residual\"">gnuplot.gp
gnuplot gnuplot.gp
ps2pdf -dEPSCrop -dNOSAFER ./Output/GMRES_RESTARTED_GAUSS_SEIDEL.ps ./Output/GMRES_RESTARTED_GAUSS_SEIDEL.pdf
mv ./Output/Error_GMRES.out ./Output/Error_GMRES_RESTARTED_GAUSS_SEIDEL_internal.out

./bin/main s3rmt3m3.mtx CG>run_cg.out
echo "set term post enh color solid
set output \"./Output/CG.ps\"
set title 'Conjugate Gradient'
set xlabel 'Iterations'
set ylabel 'Residual'
set y2label 'Error'
set y2tics
set log y
set log y2
set grid x y
plot \"./Output/Error_CG.out\" u 1:2 w l axes x1y2 lw 2 lt 1 lc rgb \"green\" title \"Error in A norm\",\
\"./Output/Error_CG.out\" u 1:3 w l axes x1y2 lw 2 lt 1 lc rgb \"blue\" title \"Error in 2 norm\",\
\"./Output/Error_CG.out\" u 1:4 w l axes x1y1 lw 0.1 lt 1 lc rgb \"red\" title \"Relative Residual\",\
\"./Output/Error_CG.out\" u 1:5 w l axes x1y1 lw 0.1 lt 1 lc rgb \"magenta\" title \"Absolute Residual\"">gnuplot.gp
gnuplot gnuplot.gp
ps2pdf -dEPSCrop -dNOSAFER ./Output/CG.ps ./Output/CG.pdf
echo "set term post enh color solid
set output \"./Output/CG_1.ps\"
set title 'Conjugate Gradient'
set xlabel 'Iterations'
set ylabel 'Residual'
set y2label 'Error'
set y2tics
set log y
set log y2
set grid x y
plot \"./Output/Error_CG.out\" u 1:2 w l axes x1y2 lw 2 lt 1 lc rgb \"green\" title \"Error in A norm\",\
\"./Output/Error_CG.out\" u 1:4 w l axes x1y1 lw 0.1 lt 1 lc rgb \"red\" title \"Relative Residual\"">gnuplot.gp
gnuplot gnuplot.gp
ps2pdf -dEPSCrop -dNOSAFER ./Output/CG_1.ps ./Output/CG_1.pdf


grep -e "Time" "./run_gmres_full_no.out"|cut -d"=" -f4|cut -d"s" -f1>>timings.out
grep -e "Time" "./run_gmres_full_jacobi.out"|cut -d"=" -f4|cut -d"s" -f1>>timings.out
grep -e "Time" "./run_gmres_full_gauss_seidel.out"|cut -d"=" -f4|cut -d"s" -f1>>timings.out
grep -e "Time" "./run_gmres_full_no.out"|cut -d"=" -f3|cut -d'T' -f1>>final_residual.out
grep -e "Time" "./run_gmres_full_jacobi.out"|cut -d"=" -f3|cut -d'T' -f1>>final_residual.out
grep -e "Time" "./run_gmres_full_gauss_seidel.out"|cut -d"=" -f3|cut -d'T' -f1>>final_residual.out
paste preconditioners timings.out final_residual.out>Times_FULL.out
grep -e "Time" "./run_gmres_restarted_no.out"|cut -d"=" -f4|cut -d"s" -f1>>timings_res.out
grep -e "Time" "./run_gmres_restarted_jacobi.out"|cut -d"=" -f4|cut -d"s" -f1>>timings_res.out
grep -e "Time" "./run_gmres_restarted_gauss_seidel.out"|cut -d"=" -f4|cut -d"s" -f1>>timings_res.out
grep -e "Time" "./run_gmres_restarted_no.out"|cut -d"=" -f3|cut -d'T' -f1>>final_residual_res.out
grep -e "Time" "./run_gmres_restarted_jacobi.out"|cut -d"=" -f3|cut -d'T' -f1>>final_residual_res.out
grep -e "Time" "./run_gmres_restarted_gauss_seidel.out"|cut -d"=" -f3|cut -d'T' -f1>>final_residual_res.out
paste preconditioners timings_res.out final_residual_res.out>Times_Restarted.out
echo "set term post enh color solid
set output \"./Output/Timings_GMRES.ps\"
set title sprintf(\"1. NO Preconditioning\t2. JACOBI Preconditioning\t3. GAUSS SEIDEL Preconditioning\")
set xlabel 'Preconditioning'
set ylabel 'Time (s)'
set grid x y
set xtics 1,1,3
set key right
plot \"./Times_FULL.out\" u 1:3 w linespoints lc rgb \"blue\" pt 4 t 'GMRES Full',\
\"./Times_Restarted.out\" u 1:3 w linespoints lt 1 lc rgb \"red\" pt 4 t 'GMRES Restarted'">gnuplot.gp
gnuplot gnuplot.gp
ps2pdf -dEPSCrop -dNOSAFER ./Output/Timings_GMRES.ps ./Output/Timings_GMRES.pdf

echo "set term post enh color solid
set output \"./Output/Timings_vs_residual_GMRES.ps\"
set title sprintf(\"[1]. NO Preconditioning\t[2]. JACOBI Preconditioning\t[3]. GAUSS SEIDEL Preconditioning\")
set xlabel 'Time (s)'
set ylabel 'Final Residual'
set yrange [2e-9:1.2e-8]
set label \"[1]\" at 1.25,9.56e-9
set label \"[1]\" at 0.29,9.72e-9
set label \"[2]\" at 0.46,9.48e-9
set label \"[2]\" at 0.035,3.72e-9
set label \"[3]\" at 0.18,8.78e-9
set label \"[3]\" at 0.038,2.44e-9
set grid x y
plot \"./Times_FULL.out\" u 3:4 w linespoints lc rgb \"blue\" pt 4 t 'GMRES Full',\
\"./Times_Restarted.out\" u 3:4 w linespoints lt 1 lc rgb \"red\" pt 4 t 'GMRES Restarted'">gnuplot.gp
gnuplot gnuplot.gp
ps2pdf -dEPSCrop -dNOSAFER ./Output/Timings_vs_residual_GMRES.ps ./Output/Timings_vs_residual_GMRES.pdf
