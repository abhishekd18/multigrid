/* 
*   Declaration of functions 
*/

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#define PI 3.14159265359

typedef struct{
		double **err, **res, **U, **f;
}Level;

void Read_Inp(  char *fileName,				//Input: Input file name
		int *n, 				//Output: Power of 2 (2^n) / Level
	  	char *cycle_type, 			//Output: Type of Cycle
	  	int *nu1,				//Output: Presmoothing factor
	  	int *nu2);		 		//Output: Postsmoothing factor

void GS_Smoother(double **U, double **f, int nx, int ny, int nu);

void Restriction(double **res_l, int nx, int ny, double **res_l_1);

void Prolongation(double **err_l_1, int nx, int ny, double **err_l);

void Solve(double **res, int nx, int ny, double **err);

void Residual(double **U, double **f, int nx, int ny, double **res);

double Norm(double **err, int nx, int ny);

double MaxNorm(double **err, int nx, int ny);

void MultiGrid(Level *L, int l, int nx, int ny, int gamma, int nu1, int nu2);

double** AllocateDynamicArray(int nRows, int nCols);

void FreeDynamicArray(double** Array, int nRows);

void Print_Array(char* fileName, double** Array, int nRows, int nCols, char *mode);

#endif
