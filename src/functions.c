#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mmio.h"
#include "functions.h"


/*
*	Print a CSR matrix in a file
*/
void Print_CSR(MTX *MAT, char* fileName){
	FILE *fid;
	fid = fopen(fileName,"w");
	mm_write_banner(fid,matcode);
    	mm_write_mtx_crd_size(fid, MAT->nrows, MAT->ncols, MAT->nz);
	for(int j=0;j<MAT->nz;j++){
		if(j<MAT->n_row_ptr)
			fprintf(fid,"%lg\t%d\t%d\n",MAT->val[j],MAT->JA[j],MAT->row_ptr[j]);
		else
			fprintf(fid,"%lg\t%d\n",MAT->val[j],MAT->JA[j]);
	}
return;
}

/* 
*	Print a CSC matrix in a file 
*/
void Print_CSC(MTX *MAT, char* fileName){
	FILE *fid;
	fid = fopen(fileName,"w");
	mm_write_banner(fid,matcode);
    	mm_write_mtx_crd_size(fid, MAT->n_col_ptr-1, MAT->n_col_ptr-1, MAT->nz);
	for(int j=0;j<MAT->nz;j++){
		if(j<MAT->n_col_ptr)
			fprintf(fid,"%lg\t%d\t%d\n",MAT->val[j],MAT->IA[j],MAT->col_ptr[j]);
		else
			fprintf(fid,"%lg\t%d\n",MAT->val[j],MAT->IA[j]);
	}
return;
}

/* 
*	Sort by rows for CSR format 
*/
void Sort(MTX *MAT){
    int swap, temp1,temp2;
    double temp3;
    do{
	swap = 0;
	for(int i=0;i<MAT->nz-1;i++){
		if((MAT->IA[i+1])<(MAT->IA[i])){
			temp1 = MAT->IA[i];
			MAT->IA[i] = MAT->IA[i+1];
			MAT->IA[i+1] = temp1;

			temp2 = MAT->JA[i];
			MAT->JA[i] = MAT->JA[i+1];
			MAT->JA[i+1] = temp2;

			temp3 = MAT->val[i];
			MAT->val[i] = MAT->val[i+1];
			MAT->val[i+1] = temp3;
			swap = 1;
		}
	}
    }while(swap!=0);

return;
}


/* 
*   Adapted from Matrix Market I/O example program read.c
*
*   NOTES:
*
*   1) Matrix Market files are always 1-based, i.e. the index of the first
*      element of a matrix is (1,1), not (0,0) as in C.  ADJUST THESE
*      OFFSETS ACCORDINGLY offsets accordingly when reading and writing 
*      to files.
*
*   2) ANSI C requires one to use the "l" format modifier when reading
*      double precision floating point numbers in scanf() and
*      its variants.  For example, use "%lf", "%lg", or "%le"
*      when reading doubles, otherwise errors will occur.
*/
void Read_Mat(char* file, MTX *MAT){
    int ret_code;

    FILE *f;   

    if ((f = fopen(file, "r")) == NULL) 
        exit(1);

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }


    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
            mm_is_sparse(matcode) )
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &MAT->nrows, &MAT->ncols, &MAT->nz)) !=0)
        exit(1);


    /* reseve memory for matrices */

    MAT->IA = (int *) malloc(MAT->nz * sizeof(int));
    MAT->JA = (int *) malloc(MAT->nz * sizeof(int));
    MAT->val = (double *) malloc(MAT->nz * sizeof(double));


    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    for(int i=0; i<MAT->nz; i++)
    {
        if(fscanf(f, "%d %d %lg\n", &MAT->IA[i], &MAT->JA[i], &MAT->val[i])==0)
		fprintf(stderr,"\nFailed to read input Matrix\n");
        MAT->IA[i]--;  /* adjust from 1-based to 0-based */
        MAT->JA[i]--;
    }

return;
}

/*
*	Read vector from a file
*/
void Read_Vec(MTX *MAT, char* fileName, double* b){
	FILE *fid;
	fid = fopen(fileName,"r");
	for(int i=0;i<MAT->ncols;i++){
        	if(fscanf(fid, "%lg\n", &b[i])==0)
			fprintf(stderr,"\nFailed to read input vector\n");
	}
return;
}

/*
*	Convert a matrix from Coordinate (COO) format to Compressed Sparse Row format
*/

void get_CSR(MTX *MAT, char* fileName){

	Sort(MAT);

	MAT->n_row_ptr = MAT->nrows+1;

	/* reseve memory for row_ptr */
	MAT->row_ptr = (int *) malloc(MAT->n_row_ptr * sizeof(int));

	MAT->row_ptr[0] = 0;

	int start=0;
	/* Find the value of row pointer */
	for(int j=1;j<MAT->n_row_ptr-1;j++){
		for(int i=start;i<MAT->nz;i++){
			if(MAT->IA[i+1] != MAT->IA[i]){
				MAT->row_ptr[j] = i+1;
				start = i+1;
				break;
			}
		}
	}

	MAT->row_ptr[MAT->n_row_ptr-1] = MAT->nz;
	
	matcode[1] = 'X';
	strcpy(crd,"CSR");
	Print_CSR(MAT,fileName);	
return;
}

/*
*	Convert a matrix from Coordinate (COO) format to Compressed Sparse Column format
*/

void get_CSC(MTX *MAT, char* fileName){

	MAT->n_col_ptr = MAT->ncols+1;

	/* reseve memory for col_ptr */
	MAT->col_ptr = (int *) malloc(MAT->n_col_ptr * sizeof(int));

	MAT->col_ptr[0] = 0;

	int start = 0;
	/* Find the value of column pointer */
	for(int j=1;j<MAT->n_col_ptr-1;j++){
		for(int i=start;i<MAT->nz;i++){
			if(MAT->JA[i+1] != MAT->JA[i]){
				MAT->col_ptr[j] = i+1;
				start = i+1;
				break;
			}
		}
	}
	
	MAT->col_ptr[MAT->n_col_ptr-1] = MAT->nz;

	matcode[1] = 'Y';
	strcpy(crd,"CSC");
	//Print_CSC(fileName);	
return;
}


/*
*	Matrix Vector product selection based on type of matrix
*/
void Mat_Vec_Mult(MTX *MAT, double* inVec, double* outVec){
	if((matcode[1] == 'X')&&(matcode[3] == 'S'))
		SYM_CSR_Mat_Vec_Mult(MAT,inVec, outVec);
	else if((matcode[1] == 'X')&&(matcode[3] == 'G'))
		CSR_Mat_Vec_Mult(MAT,inVec,outVec);
	else if((matcode[1] == 'Y')&&(matcode[3] == 'G'))
		CSC_Mat_Vec_Mult(MAT,inVec,outVec);
	else{
		fprintf(stderr,"\nDesired multiplication with the desired matrix format is not implemented.\n");
		exit(1);
	}
return;
}

/*
*	CSR:Matrix Vector product for general matrix
*/
void CSR_Mat_Vec_Mult(MTX *MAT, double* inVec, double* outVec){

	//Initialize outVec to zeros
	for(int j=0;j<MAT->n_row_ptr-1;j++)
		outVec[j] = 0.0;

	for(int j=0;j<MAT->n_row_ptr-1;j++)
		for(int i=MAT->row_ptr[j];i<MAT->row_ptr[j+1];i++)
			outVec[j] = outVec[j] + MAT->val[i] * inVec[MAT->JA[i]];

return;
}

/*
*	CSC:Matrix Vector product for general matrix
*/
void CSC_Mat_Vec_Mult(MTX *MAT, double* inVec, double* outVec){

	for(int j=0;j<MAT->n_col_ptr-1;j++)
		outVec[j] = 0.0;

	for(int j=0;j<MAT->n_col_ptr-1;j++)
		for(int i=MAT->col_ptr[j];i<MAT->col_ptr[j+1];i++)
			outVec[MAT->IA[i]] = outVec[MAT->IA[i]] + MAT->val[i] * inVec[j];

return;
}

/*
*	CSR:Matrix Vector product for symmetric matrix
*/
void SYM_CSR_Mat_Vec_Mult(MTX *MAT, double* inVec, double* outVec){

	//Initialize outVec to zeros
	for(int j=0;j<MAT->nrows;j++)
		outVec[j] = 0.0;

	for(int j=0;j<MAT->nrows;j++)
		for(int i=MAT->row_ptr[j];i<MAT->row_ptr[j+1];i++)
			outVec[j] = outVec[j] + MAT->val[i] * inVec[MAT->JA[i]];

	for(int j=0;j<MAT->nrows;j++)
		for(int i=MAT->row_ptr[j];i<MAT->row_ptr[j+1];i++){
			if(MAT->JA[i]!=j)
				outVec[MAT->JA[i]] = outVec[MAT->JA[i]] + MAT->val[i] * inVec[j];
		}
return;
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
*	Print simple array
*/
void Print_Array(char* fileName, double** Array, int nRows, int nCols){
	FILE *fid;
	fid = fopen(fileName,"w");
	for(int i=0;i<nRows;i++){
		for(int j=0;j<nCols;j++)
			fprintf(fid,"%lg\t",Array[i][j]);
		fprintf(fid,"\n");
	}
	fclose(fid);		
return;
}

/*
*	Check orthogonality of Krylov vectors
*/
void check_orthonormality(double **V, int m, int nrows, double eps){
	int Orthogonality=1, Normality=1;
	double prod;	

	for(int i=0;i<m;i++)
		for(int j=0;j<m;j++){
			prod = scalarProd(V[i], V[j], nrows);

			if((i==j)&&((pow(fabs(prod),0.5)>1+eps)||(pow(fabs(prod),0.5)<1-eps)))
				Normality = 0;
			if((i!=j)&&((prod>eps)||(prod<-eps)))	
				Orthogonality = 0;
		}

	if((Orthogonality==1)&&(Normality==1))
		fprintf(stdout,"\nKrylov Vectors are Orthonormal!\n");
	if(Orthogonality==0)
		fprintf(stdout,"\nWarning! Krylov Vectors are not Orthogonal! Error = %1.16E\n",eps);
	if(Normality==0)
		fprintf(stdout,"\nWarning! Krylov Vectors are not Normal! Error = %1.16E\n",eps);

return;
}

/*
*	Find best m corresponding to given expected minimum time
*/
int Find_Best_m(MTX* MAT, double* x0, double* xm, double* b, int m, double tol, double t_min, char* preconditioner){

	double t = 1000.0, rho = 1;
	int iter = m;
	char res_type[10] = "Relative";

	clock_t start, end;

	while(t>t_min){

		/* Initial Guess */
		for(int i=0;i<MAT->n_row_ptr-1;i++)
			x0[i] = 0.0;

		start = clock();	

		rho = GMRES(MAT, x0, xm, b, &iter, tol, res_type, preconditioner, "restarted"); 

		end = clock();

		m++;		
	
		t = (end-start)/(double)CLOCKS_PER_SEC;

		if(t>t_min)	iter = m;

		/* Stop if cannot find best m for given expected time within (no. of rows in Matrix)/8 sweeps*/
		if(m > (MAT->nrows)/10 ){
			fprintf(stderr,"\nSorry! Coudn't find the best \"m\" for given expected minimum time\nCurrent time = %lf\n",t);
			t_min += 0.01;
			m = 15;
		}
	}
	
	fprintf(stdout,"\nMinimum time for computation = %lf s for the best restart m = %d Iterations = %d\n", t, m-1, iter);

return m;
}
