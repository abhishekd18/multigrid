#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"

void MultiGrid(Level *L, int l, int nx_l, int ny_l, int gamma, int nu1, int nu2){
	
	// Presmoothing
	GS_Smoother(L[l].U, L[l].f, nx_l, ny_l, nu1);

	// Residual calculation
	Residual(L[l].U, L[l].f, nx_l, ny_l, L[l].res);

	//printf("\nnx_l = %d ny_l = %d\n",nx_l,ny_l);

	int nx_l_1 = nx_l/2;
	int ny_l_1 = ny_l/2;	

	// Interval sizes in both dimensions
	double hx = 1.0/nx_l_1, hy = 1.0/ny_l_1;	
	double hxy = 2.0*(1.0/(hx*hx) + 1.0/(hy*hy));

	// Allocate 2D arrays for residual, error, U and f at Level 'l-1'
	L[l-1].res = AllocateDynamicArray(nx_l_1+1,ny_l_1+1);
	L[l-1].err = AllocateDynamicArray(nx_l_1+1,ny_l_1+1);
	L[l-1].U = AllocateDynamicArray(nx_l_1+1,ny_l_1+1);
	L[l-1].f = AllocateDynamicArray(nx_l_1+1,ny_l_1+1);

	// Initialize
	for(int i=0;i<nx_l_1+1;i++)
		for(int j=0;j<ny_l_1+1;j++){
			L[l-1].res[i][j] = 0.0;
			L[l-1].err[i][j] = 0.0;
			L[l-1].U[i][j] = 0.0;
			L[l-1].f[i][j] = 0.0;
		}

	// Restriction operator
	Restriction(L[l].res, nx_l_1, ny_l_1, L[l-1].res);

	if(l==1){
		// Exact solve at the coarsest level
		GS_Smoother(L[l-1].res, L[l-1].err, nx_l_1, ny_l_1, 1);
	}else{
		for(int j=0;j<gamma;j++){

			// Assign current error and residual to U and f respectively for next recursive iteration
			for(int i=0;i<nx_l_1+1;i++){
				for(int j=0;j<ny_l_1+1;j++){
					L[l-1].U[i][j] = L[l-1].err[i][j];					
					L[l-1].f[i][j] = L[l-1].res[i][j];
				}
			}
			// Recursive call
			MultiGrid(L, l-1, nx_l_1, ny_l_1, gamma, nu1, nu2);
		}
	}

	// Initialize Error at level l before prolongation
	for(int i=0;i<nx_l+1;i++)
		for(int j=0;j<ny_l+1;j++)
			L[l].err[i][j] = 0.0;

	// Prolongation operator
	Prolongation(L[l-1].err, nx_l_1, ny_l_1, L[l].err);

	//printf("\nnx_l_1 = %d ny_l_1 = %d\n",nx_l_1,ny_l_1);

	// ul = ul + el
	for(int i=0;i<nx_l+1;i++)
		for(int j=0;j<ny_l+1;j++)
			L[l].U[i][j] += L[l].err[i][j];
	
	// Postsmoothing
	GS_Smoother(L[l].U, L[l].f, nx_l, ny_l, nu2);

	// Reassign U to err for the next return of Multigrid
	for(int i=0;i<nx_l+1;i++)
			for(int j=0;j<ny_l+1;j++){
				L[l].err[i][j] = L[l].U[i][j];					
				L[l].res[i][j] = L[l].f[i][j];
			}

	FreeDynamicArray(L[l-1].res,nx_l_1+1);
	FreeDynamicArray(L[l-1].err,nx_l_1+1);
	FreeDynamicArray(L[l-1].U,nx_l_1+1);
	FreeDynamicArray(L[l-1].f,nx_l_1+1);
return;
}
