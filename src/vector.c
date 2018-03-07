#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <R_ext/Utils.h>
#include <R.h>
#include "vector.h"
#include <complex.h>
#include <R_ext/Complex.h>

Rcomplex* compArray(int num) {
  
  Rcomplex *cArray = (Rcomplex *)malloc(num * sizeof(Rcomplex));
  if (!cArray)
    error("Out of memory error in compArray\n");
  return cArray;
}


Rcomplex** compMatrix(int row, int col) {

  int i;
  Rcomplex **cMatrix = (Rcomplex **)malloc(row * sizeof(Rcomplex *));
  if (!cMatrix)
    error("Out of memory error in compMatrix\n");
  for (i = 0; i < row; i++) {
     cMatrix[i] = (Rcomplex *)calloc(col, sizeof(Rcomplex));
    if (!cMatrix[i])
      error("Out of memory error in compMatrix\n");
  }
  return cMatrix;
}

void PcompMatrix(Rcomplex **cmatrix, int row, int col) {
   int i, j;
   for (i = 0; i < row; i++) {
      for (j = 0; j < col; j++) {
	 if (cmatrix[i][j].i >= 0) {
	    Rprintf("%6.4f+%6.4fi ", cmatrix[i][j].r, cmatrix[i][j].i);
	 } 
	 if (cmatrix[i][j].i < 0) {
	    Rprintf("%6.4f%6.4fi ", cmatrix[i][j].r, cmatrix[i][j].i);
	 } 
      }
      Rprintf("\n");      
   }
}

int* intArray(int num) {
  int *iArray = (int *)malloc(num * sizeof(int));
  if (!iArray)
    error("Out of memory error in intArray\n");
  return iArray;
}

void PintArray(int *ivector, int length) {
  int i;
  for (i = 0; i < length; i++)
    Rprintf("%5d\n", ivector[i]);
}

int** intMatrix(int row, int col) {
  int i;
  int **iMatrix = (int **)malloc(row * sizeof(int *));
  if (!iMatrix) 
    error("Out of memory error in intMatrix\n");
  for (i = 0; i < row; i++) {
    iMatrix[i] = (int *)malloc(col *  sizeof(int));
    if (!iMatrix[i]) 
      error("Out of memory error in intMatrix\n");
  }
  return iMatrix;
}

void PintMatrix(int **imatrix, int row, int col) {
  int i, j;
  for (i = 0; i < row; i++) {
    for (j = 0; j < col; j++)
      Rprintf("%5d", imatrix[i][j]);
    Rprintf("\n");
  }
}


double* doubleArray(int num) {
  double *dArray = (double *)malloc(num * sizeof(double));
  if (!dArray)
    error("Out of memory error in doubleArray\n");
  return dArray;
}

void PdoubleArray(double *dvector, int length) {
  int i;
  for (i = 0; i < length; i++)
    Rprintf("%14g\n", dvector[i]);
}

double** doubleMatrix(int row, int col) {
  int i;
  double **dMatrix = (double **)malloc((size_t)(row * sizeof(double *)));
  if (!dMatrix) 
    error("Out of memory error in doubleMatrix\n");
  for (i = 0; i < row; i++) {
    dMatrix[i] = (double *)malloc((size_t)(col * sizeof(double)));
    if (!dMatrix[i])
      error("Out of memory error in doubleMatrix\n");
  }
  return dMatrix;
}

void PdoubleMatrix(double **dmatrix, int row, int col) {
  int i, j;
  for (i = 0; i < row; i++) {
    for (j = 0; j < col; j++)
      Rprintf("%14g", dmatrix[i][j]);
    Rprintf("\n");
  }
}

double*** doubleMatrix3D(int x, int y, int z) {
  int i;
  double ***dM3 = (double ***)malloc(x * sizeof(double **));
  if (!dM3) 
    error("Out of memory error in doubleMatrix3D\n");
  for (i = 0; i < x; i++) 
    dM3[i] = doubleMatrix(y, z);
  return dM3;
}

void PdoubleMatrix3D(double ***dmatrix3D, int x, int y, int z) {
  int i, j, k;
  for (i = 0; i < x; i++) {
    Rprintf("First dimension = %5d\n", i);
    for (j = 0; j < y; j++) {
      for (k = 0; k < z; k++)
	Rprintf("%14g", dmatrix3D[i][j][k]);
      Rprintf("\n");
    }
  }
}

double**** doubleMatrix4D(int x, int y, int z, int zz) {
  int i;
  double ****dM4 = (double ****)malloc(x * sizeof(double ***));
  if (!dM4) 
    error("Out of memory error in doubleMatrix4D\n");
  for (i = 0; i < x; i++) 
    dM4[i] = doubleMatrix3D(y, z, zz);
  return dM4;
}

long* longArray(int num) {
  long *lArray = (long *)malloc(num * sizeof(long));
  if (!lArray)
    error("Out of memory error in longArray\n");
  return lArray;
}

void FreeMatrix(double **Matrix, int row) {
  int i;
  for (i = 0; i < row; i++)
    free(Matrix[i]);
  free(Matrix);
}

void FreeintMatrix(int **Matrix, int row) {
  int i;
  for (i = 0; i < row; i++)
    free(Matrix[i]);
  free(Matrix);
}

void FreecompMatrix(Rcomplex **Matrix, int row) {
  int i; 
  for (i = 0; i < row; i++)
    free(Matrix[i]);
  free(Matrix);
}


void Free3DMatrix(double ***Matrix, int index, int row) {
  int i;
  for (i = 0; i < index; i++)
    FreeMatrix(Matrix[i], row);
  free(Matrix);
}

void Free4DMatrix(double ****Matrix, int index, int index1, int row) {
  int i;
  for (i = 0; i < index; i++)
    Free3DMatrix(Matrix[i], index1, row);
  free(Matrix);
}
		
	
		
			
