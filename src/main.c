#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mmio.h"
#include "functions.h"

int main(int argc, char* argv[]){

	if (argc < 3){
		fprintf(stderr, "Usage: %s [martix-market-filename] [solver] [preconditioning-for GMRES]\n", argv[0]);
		exit(1);
	}
	
	// Set Input directory 
	char* in = "./Input/";

	// Set Output directory
	char* out = "./Output/";

	// Declare Matrix Structure
	MTX MAT;

	char inMatfile[50];
	strcpy(inMatfile,in);

	// Read Matrix from MatrixMarket format.
	Read_Mat(strcat(inMatfile,argv[1]), &MAT);

	char flnm0[50];
	strcpy(flnm0,out);

	// Convert into CSC format.
	get_CSC(&MAT, strcat(flnm0,"CSC.out"));

	char flnm1[50];
	strcpy(flnm1,out);

	// Convert into CSR format.
	get_CSR(&MAT, strcat(flnm1,"CSR.out"));

	// Exact Solution 
	double* x = (double*) malloc(MAT.ncols*sizeof(double));

	for(int i=0;i<MAT.ncols;i++)
		x[i] = 1.0;
	
	// Right hand side b=A*x
	double *b = (double*) malloc(MAT.ncols*sizeof(double));
	Mat_Vec_Mult(&MAT,x,b);

	// Final Solution vector
	double* xm = (double*) malloc(MAT.ncols*sizeof(double));

	// Initial Guess
	double* x0 = (double*) malloc(MAT.ncols*sizeof(double));
	for(int i=0;i<MAT.ncols;i++)
		x0[i] = 0.0;

	// Final residual, rho and tolerance, tol
	double rho = 1, tol=1e-8;

	if(strcmp(argv[2],"GMRES_FULL")==0){
		/*-------------------------------------------------Full-GMRES-setup-start--------------------------------------------------*/
		int m;					// No. of Krylov vectors. This will be output.
		char res_type[10] = "Relative";		// Monitor residual of this type for convergence
		clock_t start, end; 
		start = clock();
		rho = GMRES(&MAT, x0, xm, b, &m, tol, res_type, argv[3], "full"); 
		end = clock();
		fprintf(stdout,"\n\nIterations = %d\tFinal %s Residual = %e\tTime for computation = %lf s\n"\
				, m, res_type, rho, (end-start)/(double)CLOCKS_PER_SEC);
		/*-------------------------------------------------Full-GMRES-setup-end----------------------------------------------------*/

	}else if(strcmp(argv[2],"GMRES_RESTARTED")==0){
		/*---------------------------------------------Restarted-GMRES-setup-start-------------------------------------------------*/

		// Activate this block(line:81-83) to calculate the best restart parameter "m" and keep the rest commented
		// Initial minimum "m" for restarted GMRES to converge for given tolerence
		/*double t_min = 0.25; // Expected minimum time in seconds
		int m = 15;
		m = Find_Best_m(&MAT, x0, xm, b, m, tol, t_min, argv[3]);*/
		
		// Best found: m = 63 ~0.18 s, m = 40 ~0.26-0.27 s
		// Activate this block(line:87-94) when using GMRES Restarted with the best restart paramter
		int m = 40;				// Best m for minimum run time
		char res_type[10] = "Relative";		// Monitor residual of this type for convergence			
		clock_t start, end; 
		start = clock();
		rho = GMRES(&MAT, x0, xm, b, &m, tol, res_type, argv[3], "restarted"); 
		end = clock();
		fprintf(stdout,"\n\nIterations = %d\tFinal %s Residual = %e\tTime for computation = %lf s\n"\
				, m, res_type, rho, (end-start)/(double)CLOCKS_PER_SEC);

		/*---------------------------------------------Restarted-GMRES-setup-end---------------------------------------------------*/

	}else if(strcmp(argv[2],"CG")==0){
		/*-------------------------------------------Conjugate-Gradient-setup-start------------------------------------------------*/
		//No of iterations	
		int iter = 0;
		clock_t start, end; 
		start = clock();
		rho = CG(&MAT, x0, xm, b, &iter, tol);
		end = clock();
		fprintf(stdout,"\nIterations = %d\tFinal Relative Residual = %e\tTime for computation = %lf s\n"\
				,iter,rho,(end-start)/(double)CLOCKS_PER_SEC);
		/*-------------------------------------------Conjugate-Gradient-setup-end--------------------------------------------------*/
	}else{
		fprintf(stderr, "\nPlease recheck the choice of solver!\nAvailable options are: 1. GMRES_FULL\t2. GMRES_RESTARTED\t3. CG\n");
		exit(1);
	}

	// Print solution vector in file "x.vec"
	FILE *fid;
	char flnm4[50];
	strcpy(flnm4,out);
	fid = fopen(strcat(flnm4,"x.vec"),"w");
	for(int j=0;j<MAT.ncols;j++){
		fprintf(fid,"%1.16e\n",xm[j]);
	}
	fclose(fid);

	free(x);
	free(xm);
	free(x0);
	free(b);
return 0;
}
