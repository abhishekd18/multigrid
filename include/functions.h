/* 
*   Declaration of functions 
*/

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

typedef struct {
	int *IA, *JA, *col_ptr, *row_ptr;
	double* val;
	int nrows, ncols, nz, n_col_ptr, n_row_ptr;
}MTX;

void Print_CSR(MTX *MAT, char* fileName);

void Print_CSC(MTX *MAT, char* fileName);

void Read_Mat(char* file, MTX *MAT);

void Read_Vec(MTX *MAT, char* fileName, double* b);

void Sort(MTX *MAT);

void get_CSR(MTX *MAT, char* fileName);

void get_CSC(MTX *MAT, char* fileName);

void Mat_Vec_Mult(MTX *MAT, double* inVec, double* outVec);

void CSR_Mat_Vec_Mult(MTX *MAT, double* inVec, double* outVec);

void CSC_Mat_Vec_Mult(MTX *MAT, double* inVec, double* outVec);

void SYM_CSR_Mat_Vec_Mult(MTX *MAT, double* inVec, double* outVec);

void Back_Substitute(double** R, int nRows, int nCols, double* invec, double* outvec);

void get_diag(MTX *MAT, double* diag);

void check_orthonormality(double **V, int m, int nrows, double eps);

void LeftPreconditioning(MTX *MAT, double* vec, char* preconditioner);

double GMRES(MTX *MAT, double* x0, double* xm, double* b, int *m, double tol, char* res, char* preconditioner, char* mode);

void get_solution(MTX *MAT, double** R, double** V, double* g, int j, double* x0, double* xm);

int Find_Best_m(MTX* MAT, double* x0, double* xm, double* b, int m, double tol, double t_min, char* preconditioner);

void get_Krylov(MTX *MAT, double** V, double** H, int j, char* preconditioner);

double vecnorm(double* vec, int n);

double scalarProd(double* vec1, double* vec2, int n);

double CG(MTX *MAT, double* x0, double* xm, double* b, int *m, double tol);

double** AllocateDynamicArray(int nRows, int nCols);

void FreeDynamicArray(double** Array, int nRows);

void Print_Array(char* fileName, double** Array, int nRows, int nCols);

#endif
