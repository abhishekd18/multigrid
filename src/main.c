#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "functions.h"

int main(int argc, char* argv[]){

	if (argc < 2){
		fprintf(stderr, "Usage: %s parameters.inp\n", argv[0]);
		exit(1);
	}
	
	int lf, nu1, nu2, gamma;

	// Default Multigrid Cycle is 'W'
	char cycle_type[5] = "W";	

	// Read parameters from "parameters.inp"
	Read_Inp(argv[1], &lf, cycle_type, &nu1, &nu2); 

	// Default value of gamma corresponding to 'W' cycle
	if(strcmp(cycle_type,"W")==0)	gamma = 2;
	if(strcmp(cycle_type,"V")==0)	gamma = 1;

	// No. of intervals in both dimensions
	int Nx = pow(2,lf), Ny = pow(2,lf); 

	// Interval sizes in both dimensions
	double hx = 1.0/Nx, hy = 1.0/Ny;	
	double hxy = 2.0*(1.0/(hx*hx) + 1.0/(hy*hy));

	// Structure for storing Grid functions on each level
	Level *L = (Level*) malloc(lf*sizeof(Level));

	// Allocate arrays for 'U' and right hand side 'f'
	L[lf-1].U = AllocateDynamicArray(Nx+1,Ny+1);
	L[lf-1].f = AllocateDynamicArray(Nx+1,Ny+1);
	L[lf-1].res = AllocateDynamicArray(Nx+1,Ny+1);
	L[lf-1].err = AllocateDynamicArray(Nx+1,Ny+1);

	// Initial conditions
	for(int k=0;k<Nx+1;k++)
		for(int m=0;m<Ny+1;m++){
			L[lf-1].U[k][m] = 0.0;
			L[lf-1].f[k][m] = 8*PI*PI*sin(2*PI*k*hx)*sin(2*PI*m*hy);
			L[lf-1].res[k][m] = 0.0;
			L[lf-1].err[k][m] = 0.0;
		}
	
	FILE *fd = fopen("../Output/ConvergenceHistory.out","w");
	double Error, Residue, relResidual = 1.0;

	// Initial Residual
	Residual(L[lf-1].U, L[lf-1].f, Nx, Ny, L[lf-1].res);
	double Residue0 = MaxNorm(L[lf-1].res,Nx+1,Ny+1);

	// Allocate array for error
	double **error = AllocateDynamicArray(Nx+1,Ny+1);

	fprintf(stdout,"\n Starting Multigrid Cycle %s .......\n",cycle_type);

	// Iteration loop
	int iter=0;
	while(relResidual>1e-10){
		// Iteration counter
		iter++;

		// Multigrid
		MultiGrid(L, lf-1, Nx, Ny, gamma, nu1, nu2);
	
		// Residual at current iteration
		Residual(L[lf-1].U, L[lf-1].f, Nx, Ny, L[lf-1].res);

		// Calculate the max norm of residual
		Residue = MaxNorm(L[lf-1].res,Nx+1,Ny+1);
	
		// Calculate Relative residual
		relResidual = Residue/Residue0;

		// Find the max norm of error
		for(int k=0;k<Nx+1;k++)
			for(int m=0;m<Ny+1;m++)
				error[k][m] = fabs(L[lf-1].U[k][m]- sin(2*PI*k*hx)*sin(2*PI*m*hy));

		Error = MaxNorm(error, Nx+1, Ny+1);
	
		// Iteration	Absolute_Residual Relative_Residual Error
		fprintf(fd,"%d\t%1.16e\t%1.16e\t%1.16e\n", iter, Residue, relResidual, Error);

		fprintf(stdout, "\nIteration No. = %d\tRelative Residual = %1.16e\tError = %1.16e\n", iter, relResidual, Error);
	}
	
	fclose(fd);
	
	// Print grid function U in a file
	Print_Array("../Output/U.out", L[lf-1].U, Nx+1, Ny+1, "w");

	FreeDynamicArray(L[lf-1].U,Nx+1);
	FreeDynamicArray(L[lf-1].f,Nx+1);
	FreeDynamicArray(L[lf-1].res,Nx+1);
	FreeDynamicArray(L[lf-1].err,Nx+1);	

return 0;
}
