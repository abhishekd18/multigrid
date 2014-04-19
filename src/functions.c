#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "functions.h"

/*
*	Read input file
*/
void Read_Inp(  char *fileName,				//Input: Input file name
		int *n, 				//Output: Power of 2 (2^n)
	  	char *cycle_type, 			//Output: Type of Cycle
	  	int *nu1,				//Output: Presmoothing factor
	  	int *nu2)		 		//Output: Postsmoothing factor
{
	FILE *ifp;
	ifp = fopen(fileName,"r");

	char buf[100];

	int count=0;
	
	while((fgets(buf,100,ifp))!=NULL) {
		if((buf[0] == '#') || (buf[0] == '\n')) continue;

		count++;

		if(count==1){
			sscanf(buf, "%d\n", n);
		}else if(count==2){
			sscanf(buf, "%s\n",cycle_type);
		}else if(count==3){
			sscanf(buf, "%d\n", nu1);
		}else if(count==4){
			sscanf(buf, "%d\n", nu2);
		}else{
			fprintf(stdout,"\nWarning:Unknown parameter present in the input file! Will be ignored!\n");
			break;
		}
	}

	fclose(ifp);

return;
}

/*
*	Gauss Seidel smoother
*/
void GS_Smoother(double **U, double **f, int nx, int ny, int nu){

	// Interval sizes in both dimensions
	double hx = 1.0/nx, hy = 1.0/ny;	
	double hxy = 2.0*(1.0/(hx*hx) + 1.0/(hy*hy));
	double **error = AllocateDynamicArray(nx+1,ny+1);

	// Gauss Seidel iteration loop
	double Error, Uold;
	int iterations=nu;
	for(int iter=0;iter<nu;iter++)
		for(int k=1;k<nx;k++)
			for(int m=1;m<ny;m++)
				U[k][m] = ((U[k+1][m] + U[k-1][m])/(hx*hx) + (U[k][m+1] + U[k][m-1])/(hy*hy) + f[k][m])/hxy;
//				error[k][m] = fabs(U[k][m]- sin(2*PI*k*hx)*sin(2*PI*m*hy));
//			}
//		Error = MaxNorm(error,nx,ny);
//		printf("\nIteration No. %d Error = %1.16e\n",iter+1,Error);
//	}
return;
}

/*
*	Restriction operator
*/
void Restriction(double **res_l, int nx, int ny, double **res_l_1){
	int kk, mm;
	for(int k=1;k<nx;k++){
		kk = 2*k;
		for(int m=1;m<ny;m++){
			mm = 2*m;
			res_l_1[k][m] = 0.0625*(2*res_l[kk][mm+1] + res_l[kk+1][mm+1] + 2*res_l[kk+1][mm]\
					+ res_l[kk+1][mm-1] + 4*res_l[kk][mm] + 2*res_l[kk][mm-1]\
					+ res_l[kk-1][mm-1] + 2*res_l[kk-1][mm] + res_l[kk-1][mm+1]);
		}
	}

return;
}

/*
*	Prolongation operator
*/
void Prolongation(double **err_l_1, int nx, int ny, double **err_l){
	int kk, mm;		
	for(int k=1;k<nx;k++){
		kk = 2*k;
		for(int m=1;m<ny;m++){
			mm = 2*m;
			err_l[kk][mm+1]   += 0.5*err_l_1[k][m];
			err_l[kk+1][mm+1] += 0.25*err_l_1[k][m];
			err_l[kk+1][mm]   += 0.5*err_l_1[k][m];
			err_l[kk+1][mm-1] += 0.25*err_l_1[k][m];
			err_l[kk][mm-1]   += 0.5*err_l_1[k][m];
			err_l[kk-1][mm-1] += 0.25*err_l_1[k][m];
			err_l[kk-1][mm]   += 0.5*err_l_1[k][m];
			err_l[kk-1][mm+1] += 0.25*err_l_1[k][m];
			err_l[kk][mm] 	  = err_l_1[k][m];
		}
	}
return;
}

/*
*	Exact Solve at the coarsest level
*/
void Solve(double **res, int nx, int ny, double **err){

	// Interval sizes in both dimensions
	double hx = 1.0/nx, hy = 1.0/ny;	
	double hxy = 2.0*(1.0/(hx*hx) + 1.0/(hy*hy));
	
	err[1][1] = res[1][1]/hxy;
	
	//printf("err = %lf",err[1][1]);
return;
}

/*
*	Calculate residual
*/
void Residual(double **U, double **f, int nx, int ny, double **res){

	// Interval sizes in both dimensions
	double hx = 1.0/nx, hy = 1.0/ny;	
	double hxy = 2.0*(1.0/(hx*hx) + 1.0/(hy*hy));

	for(int k=1;k<nx;k++)
		for(int m=1;m<ny;m++)
			res[k][m] = f[k][m] + (U[k+1][m] + U[k-1][m])/(hx*hx) + (U[k][m+1] + U[k][m-1])/(hy*hy) - U[k][m]*hxy;
return;
}

/*
*	Finds the 2-norm of error matrix
*/
double Norm(double **err, int nx, int ny){
	double norm = 0.0;
	for(int i=0;i<nx;i++)
		for(int j=0;j<ny;j++)
			norm = norm + err[i][j]*err[i][j];
	norm = pow(norm,0.5);
return norm;
}

/*
*	Finds the max-norm of error matrix
*/
double MaxNorm(double **err, int nx, int ny){
	double norm = 0.0;
	for(int i=0;i<nx;i++)
		for(int j=0;j<ny;j++)
			if(fabs(err[i][j])>norm) norm = fabs(err[i][j]);
return norm;
}

/* 
*	Allocate 2D array
*/
double** AllocateDynamicArray(int nRows, int nCols){
	double** Array = (double**) malloc(nRows*sizeof(double*));
	for(int j=0;j<nRows;j++)
		Array[j] = (double*) malloc(nCols*sizeof(double));
return Array;
}

/* 
*	Free 2D array
*/
void FreeDynamicArray(double** Array, int nRows){
	for(int j=0;j<nRows;j++)
		free(Array[j]);
	free(Array);
return;
}

/*
*	Print array
*/
void Print_Array(char* fileName, double** Array, int nRows, int nCols, char *mode){
	FILE *fid;
	char M[10];
	strcpy(M,mode);
	fid = fopen(fileName,M);
	double hx = 1.0/(nCols-1), hy = 1.0/(nRows-1);
	for(int i=0;i<nRows;i++){
		for(int j=0;j<nCols;j++)
			fprintf(fid,"%1.16e\t%1.16e\t%1.16e\n",i*hx,j*hy,Array[i][j]);
		fprintf(fid,"\n");
	}
	fclose(fid);		
return;
}
