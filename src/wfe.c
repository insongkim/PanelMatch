#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <R.h>
#include <R_ext/Complex.h>
#include <complex.h>
#include "vector.h"
#include "wfe.h"


/* get absolute value for float */
double absolute(double value) {
  if (value < 0) {
    return -value;
  }
  else {
    return value;  
  }
}


/* multiply two complex numbers */
Rcomplex compMultiply(Rcomplex a, Rcomplex b) {

  Rcomplex c;
  c.r = a.r * b.r - a.i * b.i;
  c.i = a.i * b.r + a.r * b.i;

  return c;
}

/* add two complex numbers */
Rcomplex compAdd(Rcomplex a, Rcomplex b) {
  
  Rcomplex c;
  c.r = a.r + b.r;
  c.i = a.i + b.i;

  return c;
}


/* multiply two complex matrices: Matrix = a %*% b */
Rcomplex** compMultiplyMatrix(Rcomplex **a, Rcomplex **b, 
			      int rowa, int cola, int colb) {
  
  int i, j, k;
  Rcomplex temp;
  Rcomplex **Matrix = compMatrix(rowa, colb);
  
  for (i = 0; i < rowa; i++) {
    for (j = 0; j < colb; j++) {
      for (k = 0; k < cola; k++) {
	temp = compMultiply(a[i][k], b[k][j]);
	Matrix[i][j].r += temp.r;
	Matrix[i][j].i += temp.i;
      }
    }
  }
  return Matrix;
} 


/* crossproduct of two complex matrices: Matrix = a %*% t(b) */
Rcomplex** compcrossprod(Rcomplex **a, Rcomplex **b, int rowa, int cola, int rowb) {

  int i, j, k;
  Rcomplex temp1, temp2;
  Rcomplex **Matrix = compMatrix(rowa, rowb);
 
  for (i = 0; i < rowa; i++) {
    for (j = 0; j < rowb; j++) {
      Matrix[i][j].r = 0;
      Matrix[i][j].i = 0;
      for (k = 0; k < cola; k++) {
	temp1.r = b[j][k].r;
	temp1.i = b[j][k].i; /* excluded - */
	temp2 = compMultiply(a[i][k], temp1);
	Matrix[i][j].r += temp2.r;
	Matrix[i][j].i += temp2.i;
      }
    }
  } 
  return Matrix;
}




/* tcrossproduct of two complex matrices: Matrix = t(a) %*% b */
Rcomplex** comptcrossprod(Rcomplex **a, Rcomplex **b, int rowa, int cola, int colb) {

  int i, j, k;
  Rcomplex temp1, temp2;
  Rcomplex **Matrix = compMatrix(cola, colb);
 
  for (i = 0; i < cola; i++) {
    for (j = 0; j < colb; j++) {
      Matrix[i][j].r = 0;
      Matrix[i][j].i = 0;
      for (k = 0; k < rowa; k++) {
	temp1.r = a[k][i].r;
	temp1.i = a[k][i].i;  /* excluded - */
	temp2 = compMultiply(temp1, b[k][j]);
	Matrix[i][j].r += temp2.r;
	Matrix[i][j].i += temp2.i;
      }
    }
  }
  return Matrix;
}




/* Projection onto complex plane */
void ProjectionM(Rcomplex *Q_QQinv, Rcomplex *Q, 
		 Rcomplex *P1_first, Rcomplex *P1_second, 
		 Rcomplex *P1_third, Rcomplex *P1_fourth, Rcomplex *P1_fifth, 
		 Rcomplex *Y, Rcomplex *X, 
		 int *len_data, int *n_col, 
		 int *n_var, int *n_p1, int *n_p2, int *n_p3, int *n_p4,
		 Rcomplex *Projection) {
   
   
  /* temporary storages */
  int i, j, itemp;
      
  Rcomplex **QQQ_M = compMatrix(*len_data, *n_col);
  Rcomplex **Q_M = compMatrix(*len_data, *n_col);
  Rcomplex **QQQQ_M = compMatrix(*len_data, *len_data);
  Rcomplex **P1_M = compMatrix(*len_data, *len_data);
  Rcomplex **P1_M1 = compMatrix(*len_data, *n_p1);
  Rcomplex **P1_M2 = compMatrix(*len_data, (*n_p2-*n_p1));
  Rcomplex **P1_M3 = compMatrix(*len_data, (*n_p3-*n_p2));
  Rcomplex **P1_M4 = compMatrix(*len_data, (*n_p4-*n_p3));
  Rcomplex **P1_M5 = compMatrix(*len_data, (*len_data-*n_p4));
  Rcomplex **Y_M = compMatrix(*len_data,1);
  Rcomplex **X1_M = compMatrix(*len_data, (*n_var-1));
  Rcomplex **DataM = compMatrix(*len_data, *n_var);
  Rcomplex **PMatrix = compMatrix(*len_data, *len_data);
  Rcomplex **Result = compMatrix(*len_data, *n_var);
   
   

  /* packing QQQ_M, Q_M */
  itemp = 0;
  for (j = 0; j < *n_col; j++) {
    for (i = 0; i < *len_data; i++) {
      QQQ_M[i][j].r = Q_QQinv[itemp].r;
      QQQ_M[i][j].i = Q_QQinv[itemp].i;
      Q_M[i][j].r = Q[itemp].r;
      Q_M[i][j].i = Q[itemp].i;
      itemp++;
    }
  }

  /* packing QQQQ_M */
   
  QQQQ_M = compcrossprod(QQQ_M, Q_M, *len_data, *n_col, *len_data);
  /* PcompMatrix(QQQQ_M, 6, 6); */

  /* packing P1_M1, P1_M2 P1_M3*/
  itemp = 0;
  for (j = 0; j < *n_p1; j++) {
    for (i = 0; i < *len_data; i++) {
      P1_M1[i][j].r = P1_first[itemp].r;
      P1_M1[i][j].i = P1_first[itemp].i;
      itemp++;
    }
  }

  itemp = 0;
  for (j = 0; j < (*n_p2 - *n_p1); j++) {
    for (i = 0; i < *len_data; i++) {
      P1_M2[i][j].r = P1_second[itemp].r;
      P1_M2[i][j].i = P1_second[itemp].i;
      itemp++;
    }
  }
    

  itemp = 0;
  for (j = 0; j < (*n_p3 - *n_p2); j++) {
    for (i = 0; i < *len_data; i++) {
      P1_M3[i][j].r = P1_third[itemp].r;
      P1_M3[i][j].i = P1_third[itemp].i;
      itemp++;
    }
  }
    
  itemp = 0;
  for (j = 0; j < (*n_p4 - *n_p3); j++) {
    for (i = 0; i < *len_data; i++) {
      P1_M4[i][j].r = P1_fourth[itemp].r;
      P1_M4[i][j].i = P1_fourth[itemp].i;
      itemp++;
    }
  }
    
  itemp = 0;
  for (j = 0; j < (*len_data - *n_p4); j++) {
    for (i = 0; i < *len_data; i++) {
      P1_M5[i][j].r = P1_fifth[itemp].r;
      P1_M5[i][j].i = P1_fifth[itemp].i;
      itemp++;
    }
  }
    

    

  /*  /\* packing P1_M *\/ */
  /*  itemp = 0; */
  /*  for (j = 0; j < *len_data; j++) { */
  /*    for (i = 0; i < *len_data; i++) { */
  /* 	P1_M[i][j].r = P1[itemp].r; */
  /* 	P1_M[i][j].i = P1[itemp].i; */
  /* 	itemp++; */
  /*    } */
  /* } */

  /* packing P1_M */
  itemp = 0;
  for (j = 0 ; j < *len_data; j++) {
    for (i = 0; i < *len_data; i++) {
      if(i < *n_p1){
	P1_M[j][i].r = P1_M1[j][i].r;
	P1_M[j][i].i = P1_M1[j][i].i;
      } else if ( (i >= *n_p1) && (i < *n_p2) ) {
	P1_M[j][i].r = P1_M2[j][i- *n_p1].r;
	P1_M[j][i].i = P1_M2[j][i- *n_p1].i;
      } else if ( (i >= *n_p2) && (i < *n_p3) ) {
	P1_M[j][i].r = P1_M3[j][i- *n_p2].r;
	P1_M[j][i].i = P1_M3[j][i- *n_p2].i;
      } else if ( (i >= *n_p3) && (i < *n_p4) ) {
	P1_M[j][i].r = P1_M4[j][i- *n_p3].r;
	P1_M[j][i].i = P1_M4[j][i- *n_p3].i;
      } else {
	P1_M[j][i].r = P1_M5[j][i- *n_p4].r;
	P1_M[j][i].i = P1_M5[j][i- *n_p4].i;
      }
    }
  }






   
  /* packing X1_M */
  itemp = 0;
  for (j = 0; j < (*n_var-1); j++) {
    for (i = 0; i < *len_data; i++) {
      X1_M[i][j].r = X[itemp].r;
      X1_M[i][j].i = X[itemp].i;
      Y_M[i][0].r = Y[i].r;
      Y_M[i][0].i = Y[i].i;
      itemp++;
    }
  }
   
  /* Y_M is correct, should check X1_M */
  
  /* packing DataM */

  itemp = 0;
  for (j = 0 ; j < *len_data; j++) {
    for (i = 0; i < *n_var; i++) {
      if(i == 0){
	DataM[j][i].r = Y_M[j][i].r;
	DataM[j][i].i = Y_M[j][i].i;
      } else {
	DataM[j][i].r = X1_M[j][i-1].r;
	DataM[j][i].i = X1_M[j][i-1].i;
      }
    }
  }
  
  /* PcompMatrix(DataM, 10, *n_var); */
   
  /* Rprintf("done up to here\n"); */
   

  
  /* initializing PMatrix */
  for (j = 0 ; j < *len_data; j++) {
    for (i = 0; i < *len_data; i++) {
      if (i == j){
	PMatrix[j][i].r = 1 -(P1_M[j][i].r + QQQQ_M[j][i].r);
	PMatrix[j][i].i = 0-(P1_M[j][i].i + QQQQ_M[j][i].i);
      } else {
	PMatrix[j][i].r = 0-(P1_M[j][i].r + QQQQ_M[j][i].r);
	PMatrix[j][i].i = 0-(P1_M[j][i].i + QQQQ_M[j][i].i);
      }
	
    }
  }
  
  /* PcompMatrix(PMatrix, 6, 6); */
  
  /* packing Result */
  Result = compMultiplyMatrix(PMatrix, DataM, *len_data, *len_data, *n_var);
  /* PcompMatrix(Result, 6, *n_var); */
   


  /* storing results */

  itemp = 0;
  for (i = 0; i < *n_var; i++) {
    for (j = 0; j < *len_data; j++) {
      Projection[itemp].r = Result[j][i].r;
      Projection[itemp].i = Result[j][i].i;
      itemp++;
    }
  }



  FreecompMatrix(QQQ_M, *len_data);
  FreecompMatrix(Q_M, *len_data);
  FreecompMatrix(QQQQ_M, *len_data);
  FreecompMatrix(P1_M, *len_data);
  FreecompMatrix(P1_M1, *len_data);
  FreecompMatrix(P1_M2, *len_data);
  FreecompMatrix(P1_M3, *len_data);
  FreecompMatrix(P1_M4, *len_data);
  FreecompMatrix(P1_M5, *len_data);
  FreecompMatrix(X1_M, *len_data);
  FreecompMatrix(Y_M, *len_data);
  FreecompMatrix(DataM, *len_data);
  FreecompMatrix(PMatrix, *len_data);
  FreecompMatrix(Result, *len_data);
   
}





/* Omega.hat for DID model with complex numbers Eq.12*/
/* Calculates /sum_i t(X1.i)%*%u1.i%*%t(u2.i)%*%X2.i  */

void comp_OmegaHAC(Rcomplex *X_1, Rcomplex *u_1, 
		   Rcomplex *X_2, Rcomplex *u_2, int *len_data,
		   int *n_cov, int *unit_index, int *len_uniq_u_index,
		   Rcomplex *OmegaHAC) {


  /* temporary storages */
  int i, j, k, itemp;
  int n;
  int count;

  Rcomplex **X1 = compMatrix(*len_data, *n_cov);
  Rcomplex **X2 = compMatrix(*len_data, *n_cov);
  Rcomplex *u1 = compArray(*len_data);
  Rcomplex *u2 = compArray(*len_data);
  Rcomplex **OmegaMatrix = compMatrix(*n_cov, *n_cov);

  /* packing X1, X2, u1, u2 */
  itemp = 0;
  for (j = 0; j < *n_cov; j++) {
    for (i = 0; i < *len_data; i++) {
      X1[i][j].r = X_1[itemp].r;
      X1[i][j].i = X_1[itemp].i;
      X2[i][j].r = X_2[itemp].r;
      X2[i][j].i = X_2[itemp].i;
      u1[i].r = u_1[i].r;
      u1[i].i = u_1[i].i;
      u2[i].r = u_2[i].r;
      u2[i].i = u_2[i].i;

      itemp++;
    }
  }

  /* packing OmegaMatrix */

  for (j = 0; j < *n_cov; j++) {
    for (i = 0; i < *n_cov; i++) {
      OmegaMatrix[i][j].r = 0;
      OmegaMatrix[i][j].i = 0;
    }
  }
 
   
  /* calculating t(X1.i) %*% u1.i %*% t(u2.i) %*% X2.i */
   
  for (n = 0; n < *len_uniq_u_index; n++) {
    count = 0;
    /* Rprintf("unit is  %d\n", (n+1)); */
    
    for (k = 0; k < *len_data; k++) {
      if (unit_index[k] == (n+1)) {
	count ++;
      }
    }
    /* Rprintf("T_%d is  %d\n", (n+1), count); */
      
    /* Rcomplex **tX1_i = compMatrix(*n_cov, count); */
    Rcomplex **X1_i = compMatrix(count, *n_cov);
    Rcomplex **X2_i = compMatrix(count, *n_cov);
    Rcomplex **u1_i = compMatrix(count, 1);
    Rcomplex **u2_i = compMatrix(count, 1);
    Rcomplex **uu_i = compMatrix(count, count);
    Rcomplex **Xuu_i = compMatrix(*n_cov, count);
    Rcomplex **XuuX_i = compMatrix(*n_cov, *n_cov);


      
    /* packing X1_i X2_i, u1.i, u2.i*/
    itemp = 0;
    for (k = 0; k < *len_data; k++) {
      if (unit_index[k] == (n+1)) {
	for (j = 0; j < *n_cov; j++) {
	  X1_i[itemp][j].r = X1[k][j].r;
	  X1_i[itemp][j].i = X1[k][j].i;
	  X2_i[itemp][j].r = X2[k][j].r;
	  X2_i[itemp][j].i = X2[k][j].i;
	  u1_i[itemp][0].r = u1[k].r;
	  u1_i[itemp][0].i = u1[k].i;
	  u2_i[itemp][0].r = u2[k].r;
	  u2_i[itemp][0].i = u2[k].i;
	}
	itemp++;
      }
    }

    /* PcompMatrix(X1_i, count, *n_cov); */
    
    /* /\* packing tX_i *\/ */

    /* for (j = 0; j < *n_cov; j++) { */
    /*   for (k = 0; k < count; k++) { */
    /* 	tX1_i[j][k] = X1_i[k][j]; */
    /*   } */
    /* } */

    
    /*   /\* Rprintf("number of observations for unit %d is %d\n", (n+1), count ); *\/ */
   
    /* packing uu_i */
    uu_i = compcrossprod(u1_i, u2_i, count, 1, count);
    /* PcompMatrix(uu_i, count, count); */

    /* packing Xuu_i */
    Xuu_i = comptcrossprod(X1_i, uu_i, count, *n_cov, count);
    /* PcompMatrix(Xuu_i, *n_cov, count); */

    /* packing XuuX_i */
    XuuX_i = compMultiplyMatrix(Xuu_i, X2_i, *n_cov, count, *n_cov);
    /* PcompMatrix(XuuX_i, *n_cov, *n_cov); */


    /* Adding results in OmegaMatrix  */
    for(i=0; i < *n_cov; i++) {
      for(j=0; j < *n_cov; j++){
	OmegaMatrix[i][j].r +=  XuuX_i[i][j].r;
	OmegaMatrix[i][j].i +=  XuuX_i[i][j].i;
      }
    }

    /* FreecompMatrix(tX1_i, *n_cov); */
    FreecompMatrix(X1_i, count);
    FreecompMatrix(X2_i, count);
    FreecompMatrix(u1_i, count);
    FreecompMatrix(u2_i, count);
    FreecompMatrix(uu_i, count);
    FreecompMatrix(Xuu_i, *n_cov);
    FreecompMatrix(XuuX_i, *n_cov);
      
  }
   
  /* storing results */

  itemp = 0;
  for (i = 0; i < *n_cov; i++) {
    for (j = 0; j < *n_cov; j++) {
      OmegaHAC[itemp].r = OmegaMatrix[i][j].r;
      OmegaHAC[itemp].i = OmegaMatrix[i][j].i;
      itemp++;
    }
  }


  FreecompMatrix(X1, *len_data);
  FreecompMatrix(X2, *len_data);
  free(u1);
  free(u2);
  FreecompMatrix(OmegaMatrix, *n_cov);


}




/* Omega.hat for DID model with complex numbers */
/* Calculates /sum_i t(X1.i)%*%diag(u1.i%*%t(u2.i))%*%X2.i  */

void comp_OmegaHC(Rcomplex *X_1, Rcomplex *u_1, 
		  Rcomplex *X_2, Rcomplex *u_2, int *len_data,
		  int *n_cov, int *unit_index, int *len_uniq_u_index,
		  Rcomplex *OmegaHC) {
   
   
  /* temporary storages */
  int i, j, k, itemp;
  int n;
  int count;
   
  Rcomplex **X1 = compMatrix(*len_data, *n_cov);
  Rcomplex **X2 = compMatrix(*len_data, *n_cov);
  Rcomplex *u1 = compArray(*len_data);
  Rcomplex *u2 = compArray(*len_data);
  Rcomplex **OmegaMatrix = compMatrix(*n_cov, *n_cov);

  /* packing X1, X2, u1, u2 */
  itemp = 0;
  for (j = 0; j < *n_cov; j++) {
    for (i = 0; i < *len_data; i++) {
      X1[i][j].r = X_1[itemp].r;
      X1[i][j].i = X_1[itemp].i;
      X2[i][j].r = X_2[itemp].r;
      X2[i][j].i = X_2[itemp].i;
      u1[i].r = u_1[i].r;
      u1[i].i = u_1[i].i;
      u2[i].r = u_2[i].r;
      u2[i].i = u_2[i].i;

      itemp++;
    }
  }

  /* packing OmegaMatrix */

  for (j = 0; j < *n_cov; j++) {
    for (i = 0; i < *n_cov; i++) {
      OmegaMatrix[i][j].r = 0;
      OmegaMatrix[i][j].i = 0;
    }
  }
 
   
  /* calculating t(X1.i) %*% u1.i %*% t(u2.i) %*% X2.i */
   
  for (n = 0; n < *len_uniq_u_index; n++) {
    count = 0;
    /* Rprintf("unit is  %d\n", (n+1)); */
    
    for (k = 0; k < *len_data; k++) {
      if (unit_index[k] == (n+1)) {
	count ++;
      }
    }
    /* Rprintf("T_%d is  %d\n", (n+1), count); */
      
    /* Rcomplex **tX1_i = compMatrix(*n_cov, count); */
    Rcomplex **X1_i = compMatrix(count, *n_cov);
    Rcomplex **X2_i = compMatrix(count, *n_cov);
    Rcomplex **u1_i = compMatrix(count, 1);
    Rcomplex **u2_i = compMatrix(count, 1);
    Rcomplex **uu_i = compMatrix(count, count);
    Rcomplex **diag_uu_i = compMatrix(count, count);
    Rcomplex **Xuu_i = compMatrix(*n_cov, count);
    Rcomplex **XuuX_i = compMatrix(*n_cov, *n_cov);


      
    /* packing X1_i X2_i, u1.i, u2.i*/
    itemp = 0;
    for (k = 0; k < *len_data; k++) {
      if (unit_index[k] == (n+1)) {
	for (j = 0; j < *n_cov; j++) {
	  X1_i[itemp][j].r = X1[k][j].r;
	  X1_i[itemp][j].i = X1[k][j].i;
	  X2_i[itemp][j].r = X2[k][j].r;
	  X2_i[itemp][j].i = X2[k][j].i;
	  u1_i[itemp][0].r = u1[k].r;
	  u1_i[itemp][0].i = u1[k].i;
	  u2_i[itemp][0].r = u2[k].r;
	  u2_i[itemp][0].i = u2[k].i;
	}
	itemp++;
      }
    }

    /* PcompMatrix(X1_i, count, *n_cov); */
    
    /* /\* packing tX_i *\/ */

    /* for (j = 0; j < *n_cov; j++) { */
    /*   for (k = 0; k < count; k++) { */
    /* 	tX1_i[j][k] = X1_i[k][j]; */
    /*   } */
    /* } */

    
    /*   /\* Rprintf("number of observations for unit %d is %d\n", (n+1), count ); *\/ */
   
    /* packing uu_i */
    uu_i = compcrossprod(u1_i, u2_i, count, 1, count);
    /* PcompMatrix(uu_i, count, count); */


    /* packing diag_uu_i */

    for(i=0; i < count; i++) {
      for(j=0; j < count; j++){
	if (i == j){
	  diag_uu_i[i][j] = uu_i[i][j];
	} else {
	  diag_uu_i[i][j].r = 0;
	  diag_uu_i[i][j].i = 0;
	}
      }
    }
      

    /* packing Xuu_i */
    Xuu_i = comptcrossprod(X1_i, diag_uu_i, count, *n_cov, count);
    /* PcompMatrix(Xuu_i, *n_cov, count); */

    /* packing XuuX_i */
    XuuX_i = compMultiplyMatrix(Xuu_i, X2_i, *n_cov, count, *n_cov);
    /* PcompMatrix(XuuX_i, *n_cov, *n_cov); */


    /* Adding results in OmegaMatrix  */
    for(i=0; i < *n_cov; i++) {
      for(j=0; j < *n_cov; j++){
	OmegaMatrix[i][j].r +=  XuuX_i[i][j].r;
	OmegaMatrix[i][j].i +=  XuuX_i[i][j].i;
      }
    }

    /* FreecompMatrix(tX1_i, *n_cov); */
    FreecompMatrix(X1_i, count);
    FreecompMatrix(X2_i, count);
    FreecompMatrix(u1_i, count);
    FreecompMatrix(u2_i, count);
    FreecompMatrix(uu_i, count);
    FreecompMatrix(diag_uu_i, count);
    FreecompMatrix(Xuu_i, *n_cov);
    FreecompMatrix(XuuX_i, *n_cov);
      
  }
   
  /* storing results */

  itemp = 0;
  for (i = 0; i < *n_cov; i++) {
    for (j = 0; j < *n_cov; j++) {
      OmegaHC[itemp].r = OmegaMatrix[i][j].r;
      OmegaHC[itemp].i = OmegaMatrix[i][j].i;
      itemp++;
    }
  }

  FreecompMatrix(X1, *len_data);
  FreecompMatrix(X2, *len_data);
  free(u1);
  free(u2);
  FreecompMatrix(OmegaMatrix, *n_cov);


}








// To do: change C.it into double

void Index(int *index, int *uniq_index, int *len_u_index, int *len_data, int *result){
  
  // initializing index_result
  memset(result, 0, sizeof(int)*(*len_data));

  int i;
  for (i = 0 ; i < *len_u_index; i++) {
    
    int k;
    for (k = 0 ; k < *len_data; k++) {
      if ( index[k] == uniq_index[i] ) {
	result[k] = i+1;
	//Rprintf("allocated index: %d\n", result[k]);
      }
    }

  }

}



void VectorizeC(double *Wvec, int *nrow, int *ncol, int *time_index, int *dyad_index, int *n_obs, double *results) {

  int i, j, itemp;
  
  double **W = doubleMatrix(*nrow, *ncol);

  itemp = 0;
  for (j = 0; j < *ncol; j++)
    for (i = 0; i < *nrow; i++)
      W[i][j] = Wvec[itemp++];
  
  for (i = 0; i < *n_obs; i++) {
    results[i] = W[time_index[i]-1][dyad_index[i]-1];
  }
  
  FreeMatrix(W, *nrow);
}



void Transform(double *y, int *n, int *treat, double *pscore, double *ytrans) {

  int i, sumTreat;
  double psDenom0, psDenom1;

  sumTreat = 0; psDenom1 = 0; psDenom0 = 0;
  for (i = 0; i < *n; i++) {
    sumTreat += treat[i];
    if (treat[i] == 1) {
      psDenom1 += (1/pscore[i]);
    } else { 
      psDenom0 += (1/(1-pscore[i]));
    }
  }

  for (i = 0; i < *n; i++) {
    if (treat[i] == 1) {
      ytrans[i] = y[i]*sumTreat / (pscore[i] * psDenom1);
    } else {
      ytrans[i] = y[i]*(*n - sumTreat) / ((1 - pscore[i]) * psDenom0);      
    }
  }

}


// Generating time index when time index is missing
void GenTime(int* unit_index, 
	     int* len_data, /* total number of rows in the data */
	     int* len_u_index,
	     double* time_index) {
  
  int i;
  int k;
  int count;
  for (i = 0 ; i < *len_u_index; i++) {
    count = 0;
    for (k = 0 ; k < *len_data ; k++) {
      if ( unit_index[k] == (i+1) ) {
	count ++;
	time_index[k] = count;
      }
    }
  } 
}




// Demeaning numerical vector by index
void Demean(double *var,
	    int *index, 
	    int *len_index, /* unique number of unit index, e.g. number of countries*/
	    int *len_data, /* total number of rows of data */
	    double *demean) {


  double count;
  double sum_var;
  double *mean = doubleArray(*len_index);
   
  int i;
  int k;
  
  for (i = 0 ; i < *len_index; i++) {
    sum_var = 0;
    count = 0;
    for (k = 0 ; k < *len_data ; k++) {
      if ( index[k] == (i+1) ) {
	count ++;
	sum_var += var[k];
      }
    }
    mean[i] = sum_var / count;
    // Rprintf("1] Wmean_i: %f\n", Wmean[i]);
  }
  
  for (i = 0 ; i < *len_index; i++) {
    for (k = 0 ; k < *len_data ; k++) {
      if ( index[k] == (i+1) ) {
	demean[k] = var[k] - mean[i];
	//Rprintf("2] dmean_k: %f\n", demean[k]);
      }
    }
  } 

  free(mean);

}
    





// Weighted Demeaning numerical vector by index


void WDemean(double *var, double *weight,
	     int *index, 
	     int *len_index, /* unique number of unit index, e.g. number of countries*/
	     int *len_data, /* total number of rows of data */
	     double *Wdemean) {


  double sum_weight;
  double Wsum_var;
  double *Wmean = doubleArray(*len_index);
   
  int i;
  int k;
  
  for (i = 0 ; i < *len_index; i++) {
    Wsum_var = 0;
    sum_weight = 0;
    for (k = 0 ; k < *len_data ; k++) {
      if ( index[k] == (i+1) ) {
	sum_weight += weight[k];
	Wsum_var += weight[k] * var[k];
      }
    }
    Wmean[i] = Wsum_var / sum_weight;
    // Rprintf("1] Wmean_i: %f\n", Wmean[i]);
  }
  
  for (i = 0 ; i < *len_index; i++) {
    for (k = 0 ; k < *len_data ; k++) {
      if ( index[k] == (i+1) ) {
	Wdemean[k] = var[k] - Wmean[i];
	// Rprintf("2] Wmean_k: %f\n", Wdemean[k]);
      }
    }
  } 

  free(Wmean);
}
    


// Sqrt(W_it) * Weighted Demean_it


void WWDemean(double *var, double *weight,
	      int *index, 
	      int *len_index, /* unique number of unit index, e.g. number of countries*/
	      int *len_data, /* total number of rows of data */
	      double *WWdemean) {


  double sum_weight;
  double Wsum_var;
  double *Wmean = doubleArray(*len_index);
   
  int i;
  int k;
  
  for (i = 0 ; i < *len_index; i++) {
    Wsum_var = 0;
    sum_weight = 0;
    for (k = 0 ; k < *len_data ; k++) {
      if ( index[k] == (i+1) ) {
	sum_weight += weight[k];
	Wsum_var += weight[k] * var[k];
      }
    }
    if (sum_weight !=0) {
      Wmean[i] = Wsum_var / sum_weight;
    } else {
      Wmean[i] =0;
    }
    /* Rprintf("1] Wmean_i: %f\n", Wmean[i]); */
  }
  
  for (i = 0 ; i < *len_index; i++) {
    for (k = 0 ; k < *len_data ; k++) {
      if ( index[k] == (i+1) ) {
	WWdemean[k] = sqrt(weight[k]) * (var[k] - Wmean[i]);
	/* Rprintf("3] WWmean_k: %f\n", WWdemean[k]); */
      }
    }
  } 

  free(Wmean);

}



/* Two-way demean: x_it - bar{x_i} - bar{x_t} + bar{x_it} */




void TwayDemean(double *var,
		int *unit_index, int *time_index, 
		int *len_u_index, /* unique number of unit index, e.g. number of countries*/
		int *len_t_index, /* unique number of time index, e.g. number of years*/
		int *len_data, /* total number of rows of data */
		double *TwayDemean) {
  
  double count;
  double sum_var;
  double *umean = doubleArray(*len_u_index);
  double *tmean = doubleArray(*len_t_index);
  double omean;
  int i;
  int j;
  int k;
  
  /* calculating unit mean */
  for (i = 0 ; i < *len_u_index; i++) {
    sum_var = 0;
    count = 0;
    for (k = 0 ; k < *len_data ; k++) {
      if ( unit_index[k] == (i+1) ) {
	count ++;
	sum_var += var[k];
      }
    }
    umean[i] = sum_var / count;
    /* Rprintf("1] Unit mean_i: %f\n", umean[i]); */
  }


  /* calculating time mean */
  for (j = 0 ; j < *len_t_index; j++) {
    sum_var = 0;
    count = 0;
    for (k = 0 ; k < *len_data ; k++) {
      if ( time_index[k] == (j+1) ) {
	count ++;
	sum_var += var[k];
      }
    }
    tmean[j] = sum_var / count;
    /* Rprintf("2] Time mean_j: %f\n", tmean[j]); */
  }
  
  count = 0;
  sum_var = 0;
  /* calculating overall mean */
  for (k = 0 ; k < *len_data ; k++) {
    count ++;
    sum_var += var[k];
  }
  omean = sum_var / count;
  /* Rprintf("3] overall mean: %f\n", omean); */

  /* #pragma omp parallel for   */
  for (i =0 ; i < *len_u_index ; i++) {
    for (j =0 ; j < *len_t_index ; j++) {
      for (k = 0 ; k < *len_data ; k++) {
	if ( unit_index[k] == (i+1) && time_index[k] == (j+1) ) {
	  TwayDemean[k] = var[k] - umean[i] - tmean[j] + omean;	
	  
	}
      }
    }
  }

  free(umean);
  free(tmean);
}



/* A function checking index exist  */

int is_index_exist(int* unit_index, int* time_index,
		   int* len_u_index,
		   int* len_t_index, 
		   int* len_data, 
		   int** exist) {
  int i, j, k;
  int empty = 0;
  /* packing exist: 1 if exist 0 otherwise */
  
  /* initializing exist */
  for (j = 0 ; j < *len_t_index; j++) {
    for (i = 0; i < *len_u_index ; i++) {
      exist[j][i] = 0;
    }
  }
 
  /* #pragma omp parallel for   */
  for (j = 0 ; j < *len_t_index; j++) {
    for (i=0; i < *len_u_index; i++) {
      for (k = 0 ; k < *len_data; k++) {
	if( (unit_index[k] == (i+1)) && (time_index[k] == (j+1))) {
	  exist[j][i] = 1;
	  break;
	}
      }
    }
  }
  return empty;
}




int is_time_index_exist(int* u_i, int* t_i, int i, int j, int size) {
  int iter = 0;
  int exist = 0;
  for (iter = 0 ; iter < size ; iter++) {
    if (u_i[iter] == i && t_i[iter] == j) {
      exist = 1;
      break;
    }
  }
  return exist;
}



/* A function generating weights for unit weighted fixed effects */


void GenWeightsUnit(int *unit_index, int *time_index, int *tr, int *C_it,
		    int *len_data, /* total number of rows in the data */
		    int *len_u_index,
		    int *len_t_index,
		    int *ate, int *att,
		    int *verbose,  /* 1 if extra print is needed */ 
		    double *weight) {
  int i, j;
  /* Rprintf("u_index: %d, t_index: %d, len_data: %d\n", *len_u_index, *len_t_index, *len_data); */

  int** exist = intMatrix(*len_t_index, *len_u_index);
  is_index_exist(unit_index, time_index, len_u_index, len_t_index, len_data, exist);

  /* #pragma omp parallel for   */
  for (i = 0 ; i < *len_u_index ; i++) {


    if(*verbose) {
      int percent = *len_u_index/10;
      if (i % percent == 0){
	/*  Rprintf("- %d th year calculated:\n", (j+1)); */
	Rprintf(".");
	R_FlushConsole();
      }
    }


    for (j = 0 ; j < *len_t_index ; j++) {
      double w_it[*len_t_index];

      // initialize all elements in w_it to 0
      memset(w_it, 0, sizeof(double)*(*len_t_index));
      

      int c_it = 0;
      int t_it = 0;

      // initialize c_it and t_it
      int k;
      for (k = 0 ; k < *len_data ; k++) {
	if ( unit_index[k] == (i+1) && time_index[k] == (j+1) ) {
	  c_it = C_it[k];
	  t_it = tr[k];
	  /* Rprintf(" t_it: %d\n", t_it); */
	  /* Rprintf("0] unit: %d\n", i+1); */
	}
      }

  
      /* Rprintf("time index exist for unit %d and time %d: %d\n",
	 i+1, j+1, is_time_index_exist(unit_index, time_index, i+1,
	 j+1, *len_data )); */

      /* if (is_time_index_exist(unit_index, time_index, i+1, j+1, *len_data)) */
      if (exist[j][i]) {
	if (t_it == 1) { /* ifTRUE(sub[,treat]...) == 1 */
	  // count control
	  double count_control = 0;
	  int k = 0;
	  for (k = 0 ; k < *len_data ; k++) {
	    if (unit_index[k] == (i+1) && tr[k] == 0)
	      count_control++;
	  }
          
	  if (count_control > 0) { 
	    double v_it = 1 / count_control;
	    /* Rprintf("1] time: %d\n", j+1); */
	    /* Rprintf("5.1] number of control: %f\n", count_control); */
	    /* Rprintf("6] v_it: %f\n", v_it); */
	    w_it[j] = 1;
            
	    for (k = 0 ; k < *len_data ; k++) {
	      if (unit_index[k] == (i+1) && tr[k] == 0) {
		int t_index = time_index[k] - 1;
		// Rprintf("2] opposite treatment index: %d\n", t_index+1);
		w_it[t_index] = v_it;
                
	      }
	    }

	  }
	  /* PdoubleArray(w_it, *len_t_index); */
	  if (*ate == 1) {
	    int k = 0;
	    for (k = 0 ; k < *len_t_index ; k++)
	      weight[k*(*len_u_index)+i] = weight[k*(*len_u_index)+i] + (w_it[k] * c_it);
	  }
	  else if (*att == 1) {
	    int k = 0;
	    for (k = 0 ; k < *len_t_index ; k++)
	      weight[k*(*len_u_index)+i] = weight[k*(*len_u_index)+i] + (w_it[k] * c_it * t_it);
	  }
	}        
	else if (t_it == 0) { /* ifTRUE(sub[,treat]...) == 0 */
	  // count treate
	  double count_treat = 0;
	  int k = 0;
	  for (k = 0 ; k < *len_data ; k++) {
	    if (unit_index[k] == (i+1) && tr[k] == 1)
	      count_treat++;
	  }

	  if (count_treat > 0) {
	    double v_it = 1 / count_treat;
	    /* Rprintf("3] j: %d\n", j+1); */
	    /* Rprintf("5.2] number of treated: %f\n", count_treat); */
	    /* Rprintf("6] v_it: %f\n", v_it); */
	    w_it[j] = 1;
           
	    for (k = 0 ; k < *len_data ; k++) {
	      if (unit_index[k] == (i+1) && tr[k] == 1) {
		int t_index = time_index[k] - 1;
		/* Rprintf("4] opposite treatment index: %d\n", t_index+1); */
		w_it[t_index] = v_it;
	      }
	    }
	  }
	  /* PdoubleArray(w_it, *len_t_index); */
	  if (*ate == 1) {
	    int k = 0;
	    for (k = 0 ; k < *len_t_index ; k++)
	      weight[k*(*len_u_index)+i] = weight[k*(*len_u_index)+i] + (w_it[k] * c_it);
	  }
	  else if (*att == 1) {
	    int k = 0;
	    for (k = 0 ; k < *len_t_index ; k++)
	      weight[k*(*len_u_index)+i] = weight[k*(*len_u_index)+i] + (w_it[k] * c_it * t_it);
	  }

	}

      }
    }
  }
 
  FreeintMatrix(exist, *len_t_index); 
}




/* A function generating weights for time weighted fixed effects */

void GenWeightsTime(int *time_index, int *unit_index, int *tr, int *C_it,
		    int *len_data, /* total number of rows in the data */
		    int *len_t_index,
		    int *len_u_index,
		    int *ate, int *att,
		    int *verbose,  /* 1 if extra print is needed */ 
		    double *weight) {
  int i, j;
  int** exist = intMatrix(*len_t_index, *len_u_index);
  is_index_exist(unit_index, time_index, len_u_index, len_t_index, len_data, exist);
  
  /* #pragma omp parallel for   */
  for (i = 0 ; i < *len_t_index ; i++) {

    if(*verbose) {
      int percent = *len_t_index/10;
      if (i % percent == 0){
	/*  Rprintf("- %d th year calculated:\n", (j+1)); */
	Rprintf(".");
	R_FlushConsole();
      }
    }

    
    for (j = 0 ; j < *len_u_index ; j++) {
    
      double w_it[*len_u_index];

      memset(w_it, 0, sizeof(double)*(*len_u_index));
      

      int c_it = 0;
      int t_it = 0;
	 
      // initialize c_it and t_it
      int k;
      for (k = 0 ; k < *len_data ; k++) {
	if ( time_index[k] == (i+1) && unit_index[k] == (j+1) ) {
	  c_it = C_it[k];
	  t_it = tr[k];
	}
      }

      if (exist[i][j]) {
	if (t_it == 1) { /* ifTRUE(sub[,treat]...) == 1 */
	  // count control
	  double count_control = 0;
	  int k = 0;
	  for (k = 0 ; k < *len_data ; k++) {
	    if (time_index[k] == (i+1) && tr[k] == 0)
	      count_control++;
	  }
          
	  if (count_control > 0) { 
	    double v_it = 1 / count_control;
	    w_it[j] = 1;
	    for (k = 0 ; k < *len_data ; k++) {
	      if (time_index[k] == (i+1) && tr[k] == 0) {
		int u_index = unit_index[k] - 1;
		w_it[u_index] = v_it;
                
	      }
	    }

	  }
	    
	  if (*ate == 1) {
	    int k = 0;
	    for (k = 0 ; k < *len_u_index ; k++)
	      weight[i*(*len_u_index)+k] = weight[i*(*len_u_index)+k] + (w_it[k] * c_it);
	  }
	  else if (*att == 1) {
	    int k = 0;
	    for (k = 0 ; k < *len_u_index ; k++)
	      weight[i*(*len_u_index)+k] = weight[i*(*len_u_index)+k] + (w_it[k] * c_it * t_it);
	  }
	}        
	else if (t_it == 0) { /* ifTRUE(sub[,treat]...) == 0 */
	  // count treate
	  double count_treat = 0;
	  int k = 0;
	  for (k = 0 ; k < *len_data ; k++) {
	    if (time_index[k] == (i+1) && tr[k] == 1)
	      count_treat++;
	  }

	  if (count_treat > 0) {
	    double v_it = 1 / count_treat;
	    w_it[j] = 1;
	    for (k = 0 ; k < *len_data ; k++) {
	      if (time_index[k] == (i+1) && tr[k] == 1) {
		int u_index = unit_index[k] - 1;
		w_it[u_index] = v_it;
	      }
	    }
	  }
	 
	  if (*ate == 1) {
	    int k = 0;
	    for (k = 0 ; k < *len_u_index ; k++)
	      weight[i*(*len_u_index)+k] = weight[i*(*len_u_index)+k] + (w_it[k] * c_it);
	  }
	  else if (*att == 1) {
	    int k = 0;
	    for (k = 0 ; k < *len_u_index ; k++)
	      weight[i*(*len_u_index)+k] = weight[i*(*len_u_index)+k] + (w_it[k] * c_it * t_it);
	  }

	}

      }
    }
  }
  FreeintMatrix(exist, *len_t_index);  
}



// Generate Weights for first-difference(FD)


void GenWeightsFD(int* unit_index, int* time_index, int* tr, int* C_it,
		  int* len_data, /* total number of rows in the data */
		  int* len_u_index,
		  int* len_t_index,
		  int* ate, int* att,
		  int* verbose,  /* 1 if extra print is needed */ 
		  double* weightfd) {
  int i, j, k;
  int** exist = intMatrix(*len_t_index, *len_u_index);
  is_index_exist(unit_index, time_index, len_u_index, len_t_index, len_data, exist);
  
  /* #pragma omp parallel for   */
  for (j = 0 ; j < *len_t_index ; j++) {

    if(*verbose) {
      int percent = *len_t_index/10;
      if (j % percent == 0){
	/*  Rprintf("- %d th year calculated:\n", (j+1)); */
	Rprintf(".");
	R_FlushConsole();
      }
    }
      
      
    for (i = 0 ; i < *len_u_index ; i++) {

      if (j == 0) {
	double w_it[*len_t_index];
	// initialize all elements in w_it to 0
	memset(w_it, 0, sizeof(double)*(*len_t_index));
      } else {

	double w_it[*len_t_index];
	// initialize all elements in w_it to 0
	memset(w_it, 0, sizeof(double)*(*len_t_index));
	 

	int c_it;
	int t_it;
	// initialize c_it and t_it

	for (k = 0 ; k < *len_data ; k++) {
	  if ( unit_index[k] == (i+1) && time_index[k] == (j+1) ) {
	    c_it = C_it[k];          
	    t_it = tr[k];
	    break;
	    /* Rprintf("2] t_it: %d\n", t_it); */
	  }
	}
	if (exist[j][i] && exist[j-1][i]) {
	  // checking whether previous time exit, i.e, j>1
	  /* if (is_time_index_exist(unit_index, time_index, i+1, j+1, *len_data) && is_time_index_exist(unit_index, time_index, i+1, j, *len_data)) { */

	  int t_it_1;
	  // initialize t_it_1 (treatment for t-1)
	  for (k = 0 ; k < *len_data ; k++) {
	    if ( unit_index[k] == (i+1) && time_index[k] == (j) ) {
	      t_it_1 = tr[k];
	      /* Rprintf("3] t_i(t-1): %d\n", t_it_1); */
	    }
	  }


	  if ( t_it != t_it_1 ) { /* check whether treatment status changes  */
	    w_it[j] = 1;
	    w_it[j-1] = 1;            
	  }
	  /* PdoubleArray(w_it, *len_t_index); */

	  if (*ate == 1) {
	    for (k = 0 ; k < *len_t_index ; k++)
	      weightfd[k*(*len_u_index)+i] = weightfd[k*(*len_u_index)+i] + (w_it[k] * c_it);
	  }
	  else if (*att == 1) {
	    for (k = 0 ; k < *len_t_index ; k++)
	      weightfd[k*(*len_u_index)+i] = weightfd[k*(*len_u_index)+i] + (w_it[k] * c_it * t_it);
	  }
	}
      }
    } 
  }
  FreeintMatrix(exist, *len_t_index);  
}




// checking whether treatment status remains same from t-1 to t

int is_t_t1_same(int* u_i, int* t_i, int i, int j, int* tr, int size) {
  int iter = 0;
  int same = 0;
  int t_it;
  int t_it_1;

  /* #pragma omp parallel for    */
  for (iter = 0 ; iter < size ; iter++) {
    /* initialize t_it and t_it_1*/   
    if ( u_i[iter] == i && t_i[iter] == j )
      t_it = tr[iter];
    if ( u_i[iter] == i && t_i[iter] == (j-1) )
      t_it_1 = tr[iter];
  }
  if (t_it == t_it_1) {
    same = 1;
  }

  return same;
}



/* A function checking treatment is same as one period before  */

int t_t1_same(int* unit_index, int* time_index,
	      int* len_u_index,
	      int* len_t_index,
	      int* len_data,
	      int* tr,
	      int** same) {
  int i, j, k;
  int t_ij;
  int t_ij_1;
  int empty = 0;

  /* packing exist: 1 if exist 0 otherwise */
  
  /* initializing exist */
  for (j = 0 ; j < *len_t_index; j++) {
    for (i = 0; i < *len_u_index ; i++) {
      same[j][i] = 0;
    }
  }
 
  /* #pragma omp parallel for   */
  for (j = 1 ; j < *len_t_index; j++) {
    for (i=0; i < *len_u_index; i++) {
      for (k = 0 ; k < *len_data; k++) {
	if( (unit_index[k] == (i+1)) && (time_index[k] == (j))) {
	  t_ij_1 = tr[k];
	}
	if( (unit_index[k] == (i+1)) && (time_index[k] == (j+1))) {
	  t_ij = tr[k];
	}
      }
      if (t_ij == t_ij_1) {
	same[j][i] = 1;
      }
    }
  }
  return empty;
}



/* Calculate did estimate */

void CalDID(int *unit_index, int *time_index, int *tr, int *C_it, double *y,
            int *len_data, /* total number of rows in the data */
            int *len_u_index,
            int *len_t_index,
            int *ate, int *att,
            int *verbose,  /* 1 if extra print is needed */ 
	    double *did) {

  int i, j, k;
  int iter;
  double tcount = 0;
   
  // initializing did 
  memset(did, 0, sizeof(double));


  int** exist = intMatrix(*len_t_index, *len_u_index);
  is_index_exist(unit_index, time_index, len_u_index, len_t_index, len_data, exist);
  
  int** same = intMatrix(*len_t_index, *len_u_index);
  t_t1_same(unit_index, time_index, len_u_index, len_t_index, len_data, tr, same);
  
  /* PintMatrix(exist, *len_t_index, *len_u_index); */
  /* PintMatrix(same,  *len_t_index, *len_u_index); */


  /* #pragma omp parallel for */
  for (j = 1 ; j < *len_t_index ; j++) {
    /* Rprintf("- %d th year calculated:\n", (j+1)); */
    
    /* if(*verbose) { */
    /* 	 if (j % 10 == 0){ */
    /* 	    Rprintf("- %d th year calculated:\n", (j+1)); */
    /* 	    /\* Rprintf(".");  *\/ */
    /* 	    R_FlushConsole(); */
    /* 	 }  */
    /* } */
  
    for (i = 0 ; i < *len_u_index ; i++) {
       
      /* proceeds if treatment status changes  */
      if ( (exist[j][i]) && (exist[j-1][i]) && (same[j][i]!=1) ) { 
	/* Rprintf("treatment status change: unit %d time %d\n", i+1, j+1);  */
	double c_it;
	double t_it;
	double y_it; 
	double y_it1;
	double did_it = 0;
	double diff = 0;
	double count_control = 0;
	// initialize c_it, t_it, y_it
	for (k = 0 ; k < *len_data ; k++) {
	  if ( unit_index[k] == i+1 && time_index[k] == j+1 ) {
	    /* Rprintf("unit index[k]: %d\n", unit_index[k]); */
	    c_it = C_it[k];          
	    t_it = tr[k];
	    y_it = y[k];
	    break;
	    //Rprintf("1] t_it: %f\n", t_it); 
	  }
	}
	// initialize y_it1
	for (k =0 ; k < *len_data ; k++) {
	  if ( unit_index[k] == i+1 && time_index[k] == j ) {
	    y_it1 = y[k];
	  }
	}
            
	if (*ate == 1) { 
	      
	  if (t_it == 1) { /* ifTRUE(sub[,treat]...) == 1 */
	    // count control

	    int iprime;
	    int jprime;
	    double yprime;
	    double diffprime = 0;
	    double diff0;

	    diff = y_it - y_it1;
	      
	    for (k = 0 ; k < *len_data ; k++) {
	      /* Rprintf("data index: %d\n", k); */
	      if (unit_index[k] != i+1 && time_index[k] == j+1 && tr[k] == 0) {
		  
		  
		iprime = unit_index[k];
		jprime = time_index[k];
		yprime = y[k];
               
		if (exist[jprime-2][iprime-1]) {
		  if (same[jprime-1][iprime-1]) {

		    count_control++;
			  
		    double yprime1;
		    /* double diff0; */
		    for (iter = 0 ; iter < *len_data ; iter++) {
		      if (unit_index[iter] == iprime && time_index[iter] == j) {
			yprime1 = y[iter];
			diff0 = yprime - yprime1; 
				   
			/* Rprintf("yprime: %f\n", yprime); */
			/* Rprintf("yprime1: %f\n", yprime1); */
			/* Rprintf("y_it-y_it1: %f\n", diff); */
			/* Rprintf("diff0: %f\n", diff0); */
                                 
			break;
		      }
		    }
		    diffprime += diff0;			
		    /* Rprintf("diffprime: %f\n", diffprime);  */
                        
		  }

		}

		     
	      }
               
	      /* increase relevant unit counts by 1 */
		
	    }


	       
	      
	    /* Rprintf("unit %d time %d's num of counter-factual unit: %f\n", i+1, j+1, count_control); */
	
	  

	    if (count_control > 0) {
	      tcount++;
	      /* Rprintf("print of total count: %f\n", tcount); */
   
	      did_it += diff - (diffprime * 1/count_control);
	      *did += did_it;
	      /* Rprintf("did is %f\n", did[0]); */

	    }

	  } else if (t_it == 0) { /* ifTRUE(sub[,treat]...) == 1 */
	    // count control

	    int iprime;
	    int jprime;
	    double yprime;
	    double diffprime = 0;
	    double diff0;
	      
	    diff = y_it1 - y_it;
	  
	    for (k = 0 ; k < *len_data ; k++) {
	      if (unit_index[k] != (i+1) && time_index[k] == (j+1) && tr[k] == 1) {
		  
		  
		iprime = unit_index[k];
		jprime = time_index[k];
		yprime = y[k];
               
		if (exist[jprime-2][iprime-1]) {
		  if (same[jprime-1][iprime-1]) {

		    count_control++;
			  
		    double yprime1;
		    /* double diff0; */
		    for (iter = 0 ; iter < *len_data ; iter++) {
		      if (unit_index[iter] == iprime && time_index[iter] == j) {
			yprime1 = y[iter];
			diff0 = yprime1 - yprime; 
				   
			/* Rprintf("yprime: %f\n", yprime); */
			/* Rprintf("yprime1: %f\n", yprime1); */
			/* Rprintf("y_it1-y_it: %f\n", diff); */
			/* Rprintf("diff0: %f\n", diff0); */
                                 
			break;
		      }
		    }
		    diffprime += diff0;			
		    /* Rprintf("diffprime: %f\n", diffprime);  */
                        
		  }

		}

		     
	      }
               
	      /* increase relevant unit counts by 1 */
		
	    }


	       
	      
	    /* Rprintf("unit %d time %d's num of counter-factual unit: %f\n", i+1, j+1, count_control); */
	

	    if (count_control > 0) {
	      tcount++;
	      /* Rprintf("print of total count: %f\n", tcount); */
   
	      did_it += diff - (diffprime * 1/count_control);
	      *did += did_it;
	      /* Rprintf("did is %f\n", did[0]); */

	    }

	  }

	} else if (*att == 1) { 


	  if (t_it == 1) { /* ifTRUE(sub[,treat]...) == 1 */
	    // count control

	    int iprime;
	    int jprime;
	    double yprime;
	    double diffprime = 0;
	    double diff0;

	    diff = y_it - y_it1;
	      
	    for (k = 0 ; k < *len_data ; k++) {
	      if (unit_index[k] != i+1 && time_index[k] == j+1 && tr[k] == 0) {
		  
		  
		iprime = unit_index[k];
		jprime = time_index[k];
		yprime = y[k];
               
		if (exist[jprime-2][iprime-1]) {
		  if (same[jprime-1][iprime-1]) {

		    count_control++;
			  
		    double yprime1;
		    /* double diff0; */
		    for (iter = 0 ; iter < *len_data ; iter++) {
		      if (unit_index[iter] == iprime && time_index[iter] == j) {
			yprime1 = y[iter];
			diff0 = yprime - yprime1; 
				   
			/* Rprintf("yprime: %f\n", yprime); */
			/* Rprintf("yprime1: %f\n", yprime1); */
			/* Rprintf("y_it-y_it1: %f\n", diff); */
			/* Rprintf("diff0: %f\n", diff0); */
                                 
			break;
		      }
		    }
		    diffprime += diff0;			
		    /* Rprintf("diffprime: %f\n", diffprime);  */
                        
		  }

		}

		     
	      }
               
	      /* increase relevant unit counts by 1 */
		
	    }


	       
	      
	    /* Rprintf("unit %d time %d's num of counter-factual unit: %f\n", i+1, j+1, count_control); */
	
	  

	    if (count_control > 0) {
	      tcount++;
	      /* Rprintf("print of total count: %f\n", tcount); */
   
	      did_it += diff - (diffprime * 1/count_control);
	      *did += did_it;
	      /* Rprintf("did is %f\n", did[0]); */

	    }

	  }

	}
	    
      }
    }
  }
   
  /* Rprintf("total number of relevant units are %f\n", tcount); */
   
  if(tcount > 0) { 
    *did = *did * (1/tcount);
    /* Rprintf("DID estimate is %f\n", *did); */
  }
   
  FreeintMatrix(exist, *len_t_index);
  FreeintMatrix(same, *len_t_index);

}




void DemeanDID(double *var, double *weights,
	       int *unit_index, int *time_index, 
	       int *len_u_index, /* unique number of unit index, e.g. number of countries*/
	       int *len_t_index, /* unique number of time index, e.g. number of years*/
	       int *len_data, /* total number of rows of data */
	       double *DemeanDID) {
   
  double sum_wvar;
  double sum_weights;
  double *umean = doubleArray(*len_u_index);
  double *tmean = doubleArray(*len_t_index);
  double omean;
  int i;
  int j;
  int k;
  
  /* calculating unit mean */
  for (i = 0 ; i < *len_u_index; i++) {
    sum_wvar = 0;
    sum_weights = 0;
    for (k = 0 ; k < *len_data ; k++) {
      if ( unit_index[k] == (i+1) ) {
	sum_wvar += weights[k] * var[k];
	sum_weights += weights[k];
      }
    }
    /* sum_weights = floor(sum_weights * 10000000000.0) / 10000000000.0; */
    /* Rprintf("sum of weights of unit %d is zero: %d\n", (i+1), sum_weights == 0); */
    /* if (fprec(sum_weights, 10) !=0) { */
    if (sum_weights < 0.00000000001 && sum_weights > -0.00000000001) {
      umean[i] = 0;
    } else { 
      umean[i] = sum_wvar / sum_weights;      
    }
    Rprintf("1] Unit mean %d: %f\n",(i+1), umean[i]);

  }
  /* PdoubleArray(umean, *len_u_index); */

  /* calculating time mean */
  for (j = 0 ; j < *len_t_index; j++) {
    sum_wvar = 0;
    sum_weights = 0;
    for (k = 0 ; k < *len_data ; k++) {
      if ( time_index[k] == (j+1) ) {
	sum_wvar += weights[k] * var[k];
	sum_weights += weights[k];
      }
    }
      
    /* sum_weights = floor(sum_weights * 10000000000.0) / 10000000000.0; */
    /* Rprintf("sum of weights of time %d is zero: %d\n", (j+1), sum_weights == 0); */

    if (sum_weights < 0.00000000001 && sum_weights > -0.00000000001) {
      tmean[j] = 0;
    } else { 
      tmean[j] = sum_wvar / sum_weights;
    }

    Rprintf("2] Time mean %d: %f\n",(j+1), tmean[j]);
      
  }
  /* PdoubleArray(tmean, *len_t_index); */
   
   
  /* calculating overall mean */
  sum_wvar  = 0;
  sum_weights = 0;
  for (k = 0 ; k < *len_data ; k++) {
    sum_wvar += weights[k] * var[k];
    sum_weights += weights[k];
  }
  omean = sum_wvar / sum_weights;
  /* Rprintf("3] overall mean: %f\n", omean); */

  /* #pragma omp parallel for   */
  for (i =0 ; i < *len_u_index ; i++) {
    for (j =0 ; j < *len_t_index ; j++) {
      for (k = 0 ; k < *len_data ; k++) {
	if ( unit_index[k] == (i+1) && time_index[k] == (j+1) && umean[i]!= 0 && tmean[j] == 0) {
	  DemeanDID[k] = var[k] - umean[i];
	} 
	if ( unit_index[k] == (i+1) && time_index[k] == (j+1) && umean[i] == 0 && tmean[j]!=0 ) {
	  DemeanDID[k] = var[k] - tmean[j];
	} 
	if ( (unit_index[k] == (i+1) && time_index[k] == (j+1)) && umean[i] == 0 && tmean[j] == 0) {
	  DemeanDID[k] = var[k] - omean;
	} 
	if ((unit_index[k] == (i+1) && time_index[k] == (j+1)) && umean[i] != 0 && tmean[j] != 0) {
	  DemeanDID[k] = var[k] - umean[i] - tmean[j] + omean;  
	}
      }
    }
  }


  free(umean);
  free(tmean);
}


// Generate Weights for difference-in-differences (DID)


void GenWeightsDID(int* unit_index, int* time_index, int* tr, int* C_it,
		   int* len_data, /* total number of rows in the data */
		   int* len_u_index,
		   int* len_t_index,
		   int* ate, int* att,
		   int *verbose,  /* 1 if extra print is needed */ 
		   double* weightdid) {

  double **Wdid = doubleMatrix(*len_t_index, *len_u_index);   
  double **W = doubleMatrix(*len_t_index, *len_u_index);   
  int i, j, k;
  int m, n;
  int iter;

  int** exist = intMatrix(*len_t_index, *len_u_index);
  is_index_exist(unit_index, time_index, len_u_index, len_t_index, len_data, exist);
  
  int** same = intMatrix(*len_t_index, *len_u_index);
  t_t1_same(unit_index, time_index, len_u_index, len_t_index, len_data, tr, same);
  
 

  for (n = 0; n < *len_u_index; n++) {
    for (m = 0; m < *len_t_index; m++) {
      Wdid[m][n] = 0;
    }
  }


  /* #pragma omp parallel for   */
  for (j = 1 ; j < *len_t_index ; j++) {

    if(*verbose) {
      int percent = *len_t_index/10;
      if (j % percent == 0){
	/*  Rprintf("- %d th year calculated:\n", (j+1)); */
	Rprintf(".");
	R_FlushConsole();
      }
    }
   

  
    for (i = 0 ; i < *len_u_index ; i++) {

      /* initialize all elements of W to 0 */

      for (n = 0; n < *len_u_index; n++) {
	for (m = 0; m < *len_t_index; m++) {
	  W[m][n] = 0;
	}
      }

      double c_it = 0;
      double t_it = 0;
      // initialize c_it and t_it
      for (k = 0 ; k < *len_data ; k++) {
	if ( unit_index[k] == i+1 && time_index[k] == j+1 ) {
	  c_it = C_it[k];          
	  t_it = tr[k];
	  break;
	  //Rprintf("1] t_it: %f\n", t_it); 
	}
      }

      if ( (exist[j][i]==1) && (exist[j-1][i]==1) && (same[j][i]!=1) ) { 
	/* is_t_t1_same(unit_index, time_index, i+1, j+1, tr, *len_data) !=1) ) { */

	if (t_it == 1) { /* ifTRUE(sub[,treat]...) == 1 */
	  // count control
	  double count_control = 0;
	  int iprime;
	  int jprime;
	  for (k = 0 ; k < *len_data ; k++) {
	    if (unit_index[k] != i+1 && time_index[k] == j+1 && tr[k] == 0) {
	      iprime = unit_index[k];
	      jprime = time_index[k];

	      if (exist[jprime-2][iprime-1]==1){
		/* if (is_time_index_exist(unit_index, time_index, iprime, jprime-1, *len_data)) { */
		if (same[jprime-1][iprime-1]==1){
		  /* if ( (is_t_t1_same( unit_index, time_index, iprime, jprime, tr, *len_data)) == 1) { */
		  count_control++;
		} 
	      }
	    }
	  }
	  /* Rprintf("unit %d time %d's number of counter-factual units is %f\n", i+1, j+1, count_control); */
	  if (count_control > 0) {
	    double v_it = 1 / count_control;
	    W[j][i] = 1;
	    W[j-1][i] = 1;

	    for (k = 0 ; k < *len_data ; k++) {
	      if (unit_index[k] != i+1 && time_index[k] == j+1 && tr[k] == 0) {
		iprime = unit_index[k];
		jprime = time_index[k];

		if (exist[jprime-2][iprime-1]==1){
		  /* if (is_time_index_exist(unit_index, time_index, iprime, jprime-1, *len_data)) { */
		  if (same[jprime-1][iprime-1]==1){	
		    /* if (is_t_t1_same( unit_index, time_index, iprime, jprime, tr, *len_data)) { */
		    W[jprime-1][iprime-1] = v_it;
		    W[jprime-2][iprime-1] = - v_it;
  	    	  
		  }
		}
	      }
	    }
	    /* PdoubleMatrix(W, *len_t_index, *len_u_index); */
	    /* Rprintf("printed\n"); */
	  }
	}

	if (t_it == 0) { /* ifTRUE(sub[,treat]...) == 1 */
	  // count control
	  double count_treat = 0;
	  int iprime;
	  int jprime;
	  for (k = 0 ; k < *len_data ; k++) {
	    if (unit_index[k] != (i+1) && time_index[k] == (j+1) && tr[k] == 1) {
	      iprime = unit_index[k];
	      jprime = time_index[k];
	      if (exist[jprime-2][iprime-1]==1){
		/* if (is_time_index_exist(unit_index, time_index, iprime, jprime-1, *len_data)) { */
		if (same[jprime-1][iprime-1]==1){
		  /* if ( (is_t_t1_same( unit_index, time_index, iprime, jprime, tr, *len_data)) == 1) { */
		  count_treat++;
		} 
	      }
	    }
	  }
	  /* Rprintf("unit %d time %d's number of counter-factual units is %f\n", i+1, j+1, count_treat); */
	  if (count_treat > 0) {
	    double v_it = 1 / count_treat;
	    /* Rprintf("1] v_it: %f\n", v_it); */

	    W[j][i] = 1;
	    W[j-1][i] = 1;

	    for (k = 0 ; k < *len_data ; k++) {
	      if (unit_index[k] != (i+1) && time_index[k] == (j+1) && tr[k] == 1) {
		iprime = unit_index[k];
		jprime = time_index[k];

		if (exist[jprime-2][iprime-1]==1){
		  /* if (is_time_index_exist(unit_index, time_index, iprime, jprime-1, *len_data)) { */
		  if (same[jprime-1][iprime-1]==1){
		    /* if (is_t_t1_same( unit_index, time_index, iprime, jprime, tr, *len_data)) { */
		    W[jprime-1][iprime-1] = v_it;
		    W[jprime-2][iprime-1] = - v_it;
  	  
		  }
		}
	      }
	    }
	    /* PdoubleMatrix(W, *len_t_index, *len_u_index); */
	    /* Rprintf("printed\n"); */
	  }
	}

  	
	if (*ate == 1) {
	  for (n = 0 ; n < *len_u_index ; n++) {
	    for (m = 0 ; m < *len_t_index ; m++) {
	      Wdid[m][n] = Wdid[m][n] + (c_it*W[m][n]);
	    }
	  }
	  //Rprintf("Print Wdid %d okay?\n", 1);
	  //PdoubleMatrix(Wdid, *len_t_index, *len_u_index);     
	}
	    
	else if (*att == 1) {
	  for (n = 0 ; n < *len_u_index ; n++) {
	    for (m = 0 ; m < *len_t_index ; m++) {
	      Wdid[m][n] = Wdid[m][n] + (c_it*t_it*W[m][n]);
	    }
	  }
	  //Rprintf("Print Wdid %d okay?\n", 1);
	  //PdoubleMatrix(Wdid, *len_t_index, *len_u_index);    
	}

      }
	 
    }
      
  } 
   
  
  /* PdoubleMatrix(Wdid, *len_t_index, *len_u_index);     */
  
  iter = 0;
  for (m = 0; m < *len_t_index; m++) {
    for (n = 0; n < *len_u_index; n++) {
      weightdid[iter] = Wdid[m][n];
      iter++;
    }
  }  

  FreeMatrix(Wdid, *len_t_index);
  FreeMatrix(W, *len_t_index);
  FreeintMatrix(exist, *len_t_index);  
  FreeintMatrix(same, *len_t_index);


}




/* // Generate Weights for "matched" difference-in-differences (DID): */
/* // Nearest Neighbor Matching on the pre-treatment outcome */
/* void GenWeightsMDID(int* unit_index, int* time_index, int* tr, int* C_it, */
/* 		    double *y, /\* needed for checking pretreatment outcome *\/ */
/* 		    double *maxdev, /\* user provided maximum deviation */
/* 				       in past outcome NULL if */
/* 				       negative *\/ */
/* 		    int* len_data, /\* total number of rows in the data *\/ */
/* 		    int* len_u_index, */
/* 		    int* len_t_index, */
/* 		    int* ate, int* att, */
/* 		    int *verbose,  /\* 1 if extra print is needed *\/  */
/* 		    double* weightdid) { */

/*   double **Wdid = doubleMatrix(*len_t_index, *len_u_index);    */
/*   double **W = doubleMatrix(*len_t_index, *len_u_index);    */
/*   int i, j, k; */
/*   int m, n; */
/*   int iter; */

/*   int** exist = intMatrix(*len_t_index, *len_u_index); */
/*   is_index_exist(unit_index, time_index, len_u_index, len_t_index, len_data, exist); */
  
/*   int** same = intMatrix(*len_t_index, *len_u_index); */
/*   t_t1_same(unit_index, time_index, len_u_index, len_t_index, len_data, tr, same); */
  
/*   /\* initialize the Weight Matrix *\/ */
/*   for (n = 0; n < *len_u_index; n++) { */
/*     for (m = 0; m < *len_t_index; m++) { */
/*       Wdid[m][n] = 0; */
/*     } */
/*   } */

/*   /\* #pragma omp parallel for   *\/ */
/*   for (j = 1 ; j < *len_t_index ; j++) { */

/*     if(*verbose) { */
/*       int percent = *len_t_index/10; */
/*       if (j % percent == 0){ */
/* 	/\*  Rprintf("- %d th year calculated:\n", (j+1)); *\/ */
/* 	Rprintf("."); */
/* 	R_FlushConsole(); */
/*       } */
/*     } */
  
/*     for (i = 0 ; i < *len_u_index ; i++) { */

/*       /\* initialize all elements of W to 0 *\/ */
/*       for (n = 0; n < *len_u_index; n++) { */
/* 	for (m = 0; m < *len_t_index; m++) { */
/* 	  W[m][n] = 0; */
/* 	} */
/*       } */

/*       double c_it = 0; */
/*       double t_it = 0; */
/*       double y_it;  */
/*       double y_it1; */


/*       /\* IF observation (i,t) and (i,t-1) exists and their treatment status differs *\/ */
/*       if ( (exist[j][i]==1) && (exist[j-1][i]==1) && (same[j][i]!=1) ) {  */
	 
/* 	// initialize c_it, t_it, y_it */
/* 	for (k = 0 ; k < *len_data ; k++) { */
/* 	  if ( unit_index[k] == i+1 && time_index[k] == j+1 ) { */
/* 	    c_it = C_it[k];           */
/* 	    t_it = tr[k]; */
/* 	    y_it = y[k]; */
/* 	    break; */
/* 	    //Rprintf("1] t_it: %f\n", t_it);  */
/* 	  } */
/* 	} */
/* 	// initialize y_it1 */
/* 	for (k =0 ; k < *len_data ; k++) { */
/* 	  if ( unit_index[k] == i+1 && time_index[k] == j ) { */
/* 	    y_it1 = y[k]; */
/* 	    break; */
/* 	  } */
/* 	} */


/* 	if (t_it == 1) { /\* ifTRUE(sub[,treat]...) == 1 *\/ */
/* 	  // count control */
/* 	  double count_control = 0; */
/* 	  int iprime; */
/* 	  int jprime; */

/* 	  /\* /\\* place holder for minium difference between y_{i,t-1} vs y_{i^\prime,t-1} *\\/ *\/ */
/* 	  double tmpydiff; */
/* 	  /\* double minydiff; *\/ */

/* 	  /\* /\\* finding nearest neighbor  *\\/ *\/ */
/* 	  /\* if (*maxdev < 0) { *\/ */
/* 	  /\*   /\\* nearest neighbor *\\/ *\/ */
/* 	  /\*   int matchiprime; *\/ */
/* 	  /\*   int matchjprime; *\/ */

/* 	  /\*   for (k = 0 ; k < *len_data ; k++) { *\/ */
/* 	  /\*     if (unit_index[k] != i+1 && time_index[k] == j && tr[k] == 0) { *\/ */
/* 	  /\* 	/\\* k corresponds to the index for the t-1 observation *\/ */
/* 	  /\* 	   of the matched observation *\\/ *\/ */
/* 	  /\* 	iprime = unit_index[k]; *\/ */
/* 	  /\* 	jprime = time_index[k]; *\/ */

/* 	  /\* 	/\\* observation for time t exists and the treatment status is control *\\/ *\/ */
/* 	  /\* 	if (exist[jprime][iprime-1]==1 && (same[jprime][iprime-1]==1)){ *\/ */
/* 	  /\* 	  tmpydiff = absolute(y_it1 - y[k]); /\\* absolute diff *\\/ *\/ */

/* 	  /\* 	  if(count_control == 0){ /\\* first matched observation *\\/ *\/ */
/* 	  /\* 	    minydiff = tmpydiff; *\/ */
/* 	  /\* 	    matchiprime = iprime; *\/ */
/* 	  /\* 	    matchjprime = jprime; *\/ */
/* 	  /\* 	    count_control++; *\/ */
/* 	  /\* 	  } else if (count_control > 0 && tmpydiff<= minydiff) { *\/ */
/* 	  /\* 	    minydiff = tmpydiff; *\/ */
/* 	  /\* 	    matchiprime = iprime; *\/ */
/* 	  /\* 	    matchjprime = jprime; *\/ */
/* 	  /\* 	    count_control++; *\/ */
/* 	  /\* 	  } else if (count_control > 0 && tmpydiff > minydiff) { *\/ */
/* 	  /\* 	    ; *\/ */
/* 	  /\* 	  } *\/ */
		     
/* 	  /\* 	} *\/ */
/* 	  /\*     } *\/ */
/* 	  /\*   } *\/ */

/* 	  /\*   /\\* if match is found *\\/ *\/ */
/* 	  /\*   if (count_control > 0 && exist[matchjprime][matchiprime-1]==1 && (same[matchjprime][matchiprime-1]==1)) { *\/ */
/* 	  /\*     W[j][i] = 1.0; *\/ */
/* 	  /\*     W[j-1][i] = 1.0; *\/ */
/* 	  /\*     W[matchjprime][matchiprime-1] = 1.0;		  *\/ */
/* 	  /\*     W[matchjprime-1][matchiprime-1] = -1.0; *\/ */

/* 	  /\*     /\\* Rprintf("Unit %d Time %d outcome %f: Nearest Neighbor is Unit %d Time %d difference %f\n", *\\/ *\/ */
/* 	  /\*     /\\* 	       (i+1), (j+1), y_it, (matchiprime), (matchjprime+1), minydiff); *\\/ *\/ */

/* 	  /\*   } *\/ */
/* 	  /\* } *\/ */

/* 	  /\* finding nearest neighbor  *\/ */
/* 	  if (*maxdev < 0) { */
/* 	    int iprime_idx; */
/* 	    double *dist_to_yit1 = doubleArray(*len_u_index); */

/* 	    /\* packing dist_to_yit1 *\/ */
/* 	    for (k = 0; k < *len_u_index; k++) { */
/* 	      dist_to_yit1[k] = -1; /\* negative value for irrelevant matched units *\/ */
/* 	    } */

/* 	    /\* filling in dist_to_yit1 with distance when treatment */
/* 	       status stays the same *\/ */
/* 	    for (k = 0 ; k < *len_data ; k++) { */
/* 	      if (unit_index[k] != i+1 && time_index[k] == j && tr[k] == 0) { */

/* 		iprime = unit_index[k]; */
/* 		jprime = time_index[k]; */
		
/* 		if (exist[jprime][iprime-1]==1 && (same[jprime][iprime-1]==1)){  */
/* 		  dist_to_yit1[iprime-1] = absolute(y_it1 - y[k]); */
/* 		} */

/* 	      } */
/* 	    } */
/* 	    /\* Now finding the index with minumum distance *\/ */

/* 	    iprime_idx = 0; */
/* 	    double min = dist_to_yit1[iprime_idx]; */
/* 	    for (k = 0; k < *len_u_index; k++) { */
/* 	      if (dist_to_yit1[k] >=0 && dist_to_yit1[k] <= min){ */
/* 		min = (double)dist_to_yit1[k]; */
/* 		iprime_idx = k; */
/* 		count_control++; */
/* 		Rprintf("-- Min: Unit %d Time %d outcome %f\n", (iprime_idx+1), (j), (min)); */
/* 	      } */
/* 	    } */
	    
/* 	    /\* if match is found *\/ */
/* 	    if (count_control > 0) { */
/* 	      W[j][i] = 1.0; */
/* 	      W[j-1][i] = 1.0; */
/* 	      W[j][iprime_idx] = 1.0; */
/* 	      W[j-1][iprime_idx] = -1.0; */

/* 	      Rprintf("Unit %d Time %d outcome %f: Nearest Neighbor is Unit %d Time %d\n", */
/* 		      (i+1), (j), y_it1, (iprime_idx+1), (j)); */

/* 	    } */
/* 	  } */
	  
/* 	  else if (*maxdev >=0) { */

/* 	  for (k = 0 ; k < *len_data ; k++) { */
/* 	    if (unit_index[k] != i+1 && time_index[k] == j && tr[k] == 0) { */
/* 	      /\* k corresponds to the index for the t-1 observation */
/* 		 of the matched observation *\/ */
/* 	      iprime = unit_index[k]; */
/* 	      jprime = time_index[k]; */
		     
/* 	      if (exist[jprime][iprime-1]==1 && (same[jprime][iprime-1]==1)){ */
/* 		tmpydiff = absolute(y_it1 - y[k]); /\* absolute diff *\/ */
/* 		if(tmpydiff <= *maxdev) { */
/* 		  count_control++; */
/* 		} */
/* 	      } */
/* 	    } */
/* 	  } */

/* 	  /\* if match is found *\/ */
/* 	  if (count_control > 0) { */
/* 	    double v_it = 1/count_control; */
/* 	    W[j][i] = 1.0; */
/* 	    W[j-1][i] = 1.0; */

/* 	    for (k = 0 ; k < *len_data ; k++) { */
/* 	      if (unit_index[k] != i+1 && time_index[k] == j && tr[k] == 0) { */
/* 		iprime = unit_index[k]; */
/* 		jprime = time_index[k]; */

/* 		if (exist[jprime][iprime-1]==1 && (same[jprime][iprime-1]==1)){ */
/* 		  tmpydiff = absolute(y_it1 - y[k]); /\* absolute diff *\/ */
/* 		  if(tmpydiff <= *maxdev) { */
/* 		    W[jprime][iprime-1] =  v_it; */
/* 		    W[jprime-1][iprime-1] = - v_it; */
/* 		    /\* Rprintf("Unit %d Time %d outcome %f: Matched to Unit %d Time %d: the difference is within %f\n", *\/ */
/* 		    /\* 	   (i+1), (j), y_it1, (iprime), (jprime), *maxdev); *\/ */
/* 		  } */
/* 		} */
/* 	      } */
/* 	    } */
/* 	  } */
	       
/* 	} */
	     
	     
/*       } */

/*       if (t_it == 0) {  */
/* 	// count treat */
/* 	double count_treat = 0; */
/* 	int iprime; */
/* 	int jprime; */

/* 	/\* place holder for minium difference between y_{i,t-1} vs y_{i^\prime,t-1} *\/ */
/* 	double tmpydiff; */
	
/* 	/\* double minydiff; *\/ */

/* 	/\* /\\* finding nearest neighbor  *\\/ *\/ */
/* 	/\* if (*maxdev < 0) { *\/ */
/* 	/\*   /\\* nearest neighbor *\\/ *\/ */
/* 	/\*   int matchiprime; *\/ */
/* 	/\*   int matchjprime; *\/ */
	       
/* 	/\*   for (k = 0 ; k < *len_data ; k++) { *\/ */
/* 	/\*     if (unit_index[k] != i+1 && time_index[k] == j && tr[k] == 1) { *\/ */
/* 	/\*       iprime = unit_index[k]; *\/ */
/* 	/\*       jprime = time_index[k]; *\/ */

/* 	/\*       if (exist[jprime][iprime-1]==1 && (same[jprime][iprime-1]==1)){ *\/ */
/* 	/\* 	tmpydiff = absolute(y_it1 - y[k]); /\\* absolute diff *\\/ *\/ */

/* 	/\* 	if(count_treat == 0){ /\\* first matched observation *\\/ *\/ */
/* 	/\* 	  minydiff = tmpydiff; *\/ */
/* 	/\* 	  matchiprime = iprime; *\/ */
/* 	/\* 	  matchjprime = jprime; *\/ */
/* 	/\* 	  count_treat++; *\/ */
/* 	/\* 	} else if (count_treat > 0 && tmpydiff<= minydiff) { *\/ */
/* 	/\* 	  minydiff = tmpydiff; *\/ */
/* 	/\* 	  matchiprime = iprime; *\/ */
/* 	/\* 	  matchjprime = jprime; *\/ */
/* 	/\* 	  count_treat++; *\/ */
/* 	/\* 	} else if (count_treat > 0 && tmpydiff > minydiff) { *\/ */
/* 	/\* 	  ; *\/ */
/* 	/\* 	} *\/ */

/* 	/\*       } *\/ */
/* 	/\*     } *\/ */
/* 	/\*   } *\/ */

/* 	/\*   /\\* if match is found *\\/ *\/ */
/* 	/\*   if (count_treat > 0 && exist[matchjprime][matchiprime-1]==1 && (same[matchjprime][matchiprime-1]==1)) { *\/ */
/* 	/\*     W[j][i] = 1.0; *\/ */
/* 	/\*     W[j-1][i] = 1.0; *\/ */
/* 	/\*     W[matchjprime][matchiprime-1] = 1.0; *\/ */
/* 	/\*     W[matchjprime-1][matchiprime-1] = -1.0; *\/ */
/* 	/\*     /\\* Rprintf("Unit %d Time %d outcome %f: Nearest Neighbor is Unit %d Time %d difference %f\n", *\\/ *\/ */
/* 	/\*     /\\* 	       (i+1), (j+1), y_it, (matchiprime-1), (matchjprime+1), minydiff); *\\/ *\/ */
	       
/* 	/\*   } *\/ */

/* 	/\* } *\/ */


/* 	/\* finding nearest neighbor  *\/ */
/* 	if (*maxdev < 0) { */
/* 	  int iprime_idx; */
/* 	  double *dist_to_yit1 = doubleArray(*len_u_index); */

/* 	  /\* packing dist_to_yit1 *\/ */
/* 	  for (k = 0; k < *len_u_index; k++) { */
/* 	    dist_to_yit1[k] = -1; /\* negative value for irrelevant matched units *\/ */
/* 	  } */

/* 	  /\* filling in dist_to_yit1 with distance when treatment */
/* 	     status stays the same *\/ */
/* 	  for (k = 0 ; k < *len_data ; k++) { */
/* 	    if (unit_index[k] != i+1 && time_index[k] == j && tr[k] == 1) { */

/* 	      iprime = unit_index[k]; */
/* 	      jprime = time_index[k]; */
		
/* 	      if (exist[jprime][iprime-1]==1 && (same[jprime][iprime-1]==1)){  */
/* 		dist_to_yit1[iprime-1] = absolute(y_it1 - y[k]); */
/* 	      } */

/* 	    } */
/* 	  } */
/* 	  /\* Now finding the index with minumum distance *\/ */

/* 	  iprime_idx = 0; */
/* 	  double min = dist_to_yit1[iprime_idx]; */
/* 	  for (k = 0; k < *len_u_index; k++) { */
/* 	    if (dist_to_yit1[k] >=0 && dist_to_yit1[k] <= min){ */
/* 	      min = (double)dist_to_yit1[k]; */
/* 	      iprime_idx = k; */
/* 	      count_treat++; */
/* 	      Rprintf("-- Min: Unit %d Time %d outcome %f\n", (iprime_idx+1), (j), (min)); */
	      
/* 	    } */
/* 	  } */
	    
/* 	  /\* if match is found *\/ */
/* 	  if (count_treat > 0) { */
/* 	    W[j][i] = 1.0; */
/* 	    W[j-1][i] = 1.0; */
/* 	    W[j][iprime_idx] = 1.0; */
/* 	    W[j-1][iprime_idx] = -1.0; */

/* 	    Rprintf("Unit %d Time %d outcome %f: Nearest Neighbor is Unit %d Time %d\n", */
/* 		    (i+1), (j), y_it1, (iprime_idx+1), (j)); */

/* 	  } */
/* 	} */


	
/* 	else if (*maxdev >=0) { */
	       
/* 	  for (k = 0 ; k < *len_data ; k++) { */
/* 	    if (unit_index[k] != i+1 && time_index[k] == j && tr[k] == 1) { */
/* 	      /\* k corresponds to the index for the t-1 observation */
/* 		 of the matched observation *\/ */
/* 	      iprime = unit_index[k]; */
/* 	      jprime = time_index[k]; */
		     
/* 	      if (exist[jprime][iprime-1]==1 && (same[jprime][iprime-1]==1)){ */
/* 		tmpydiff = absolute(y_it1 - y[k]); /\* absolute diff *\/ */
/* 		if(tmpydiff<= *maxdev) { */
/* 		  count_treat++; */
/* 		} */
/* 	      } */
/* 	    } */
/* 	  } */

/* 	  /\* if match is found *\/ */
/* 	  if (count_treat > 0) { */
/* 	    double v_it = 1/count_treat; */
/* 	    W[j][i] = 1.0; */
/* 	    W[j-1][i] = 1.0; */

/* 	    for (k = 0 ; k < *len_data ; k++) { */
/* 	      if (unit_index[k] != i+1 && time_index[k] == j && tr[k] == 1) { */
/* 		iprime = unit_index[k]; */
/* 		jprime = time_index[k]; */

/* 		if (exist[jprime][iprime-1]==1 && (same[jprime][iprime-1]==1)){ */
/* 		  tmpydiff = absolute(y_it1 - y[k]); /\* absolute diff *\/ */
/* 		  if(tmpydiff <= *maxdev) { */
/* 		    W[jprime][iprime-1] =  v_it; */
/* 		    W[jprime-1][iprime-1] = - v_it; */
/* 		    /\* Rprintf("Unit %d Time %d outcome %f: Matched to Unit %d Time %d: the difference is within %f\n", *\/ */
/* 		    /\* 	   (i+1), (j), y_it1, (iprime), (jprime), *maxdev); *\/ */
/* 		  } */
/* 		} */
/* 	      } */
/* 	    } */
/* 	  } */
/* 	} */
/*       } */


	     	
/*       if (*ate == 1) { */
/* 	for (n = 0 ; n < *len_u_index ; n++) { */
/* 	  for (m = 0 ; m < *len_t_index ; m++) { */
/* 	    Wdid[m][n] = Wdid[m][n] + (c_it*W[m][n]); */
/* 	  } */
/* 	} */
/* 	//Rprintf("Print Wdid %d okay?\n", 1); */
/* 	//PdoubleMatrix(Wdid, *len_t_index, *len_u_index);      */
/*       }  */
/*       else if (*att == 1) { */
/* 	for (n = 0 ; n < *len_u_index ; n++) { */
/* 	  for (m = 0 ; m < *len_t_index ; m++) { */
/* 	    Wdid[m][n] = Wdid[m][n] + (c_it*t_it*W[m][n]); */
/* 	  } */
/* 	} */
/* 	//Rprintf("Print Wdid %d okay?\n", 1); */
/* 	//PdoubleMatrix(Wdid, *len_t_index, *len_u_index);     */
/*       } */
	   
/*     } */
	 
/*   } */
      
/* }  */
   
  
/* /\* PdoubleMatrix(Wdid, *len_t_index, *len_u_index); *\/ */
   
/* iter = 0; */
/* for (m = 0; m < *len_t_index; m++) { */
/*   for (n = 0; n < *len_u_index; n++) { */
/*     weightdid[iter] = Wdid[m][n]; */
/*     iter++; */
/*   } */
/*  }   */
   
/* FreeMatrix(Wdid, *len_t_index); */
/* FreeMatrix(W, *len_t_index); */
/* FreeintMatrix(exist, *len_t_index);   */
/* FreeintMatrix(same, *len_t_index); */


/* } */




// Generate Weights for "matched" difference-in-differences (DID):
// Nearest Neighbor Matching on the pre-treatment outcome
void GenWeightsMDID(int* unit_index, int* time_index, int* tr, int* C_it,
		    double *y, /* needed for checking pretreatment outcome */
		    double *maxdev, /* user provided maximum deviation
				       in past outcome NULL if
				       negative */
		    int* len_data, /* total number of rows in the data */
		    int* len_u_index,
		    int* len_t_index,
		    int* ate, int* att,
		    int *verbose,  /* 1 if extra print is needed */ 
		    double* weightdid) {

  double **Wdid = doubleMatrix(*len_t_index, *len_u_index);   
  double **W = doubleMatrix(*len_t_index, *len_u_index);   
  int i, j, k;
  int m, n;
  int iter;

  int** exist = intMatrix(*len_t_index, *len_u_index);
  is_index_exist(unit_index, time_index, len_u_index, len_t_index, len_data, exist);
  
  int** same = intMatrix(*len_t_index, *len_u_index);
  t_t1_same(unit_index, time_index, len_u_index, len_t_index, len_data, tr, same);
  
  /* initialize the Weight Matrix */
  for (n = 0; n < *len_u_index; n++) {
    for (m = 0; m < *len_t_index; m++) {
      Wdid[m][n] = 0;
    }
  }

  /* #pragma omp parallel for   */
  for (j = 1 ; j < *len_t_index ; j++) {

    if(*verbose) {
      int percent = *len_t_index/10;
      if (j % percent == 0){
	/*  Rprintf("- %d th year calculated:\n", (j+1)); */
	Rprintf(".");
	R_FlushConsole();
      }
    }
  
    for (i = 0 ; i < *len_u_index ; i++) {

      /* initialize all elements of W to 0 */
      for (n = 0; n < *len_u_index; n++) {
	for (m = 0; m < *len_t_index; m++) {
	  W[m][n] = 0;
	}
      }

      double c_it = 0;
      double t_it = 0;
      double y_it; 
      double y_it1;


      /* IF observation (i,t) and (i,t-1) exists and their treatment status differs */
      if ( (exist[j][i]==1) && (exist[j-1][i]==1) && (same[j][i]!=1) ) { 
	 
	// initialize c_it, t_it, y_it
	for (k = 0 ; k < *len_data ; k++) {
	  if ( unit_index[k] == i+1 && time_index[k] == j+1 ) {
	    c_it = C_it[k];          
	    t_it = tr[k];
	    y_it = y[k];
	    break;
	    //Rprintf("1] t_it: %f\n", t_it); 
	  }
	}
	// initialize y_it1
	for (k =0 ; k < *len_data ; k++) {
	  if ( unit_index[k] == i+1 && time_index[k] == j ) {
	    y_it1 = y[k];
	    break;
	  }
	}


	if (t_it == 1) { /* ifTRUE(sub[,treat]...) == 1 */
	  // count control
	  double count_control = 0;
	  int iprime;
	  int jprime;

	  /* /\* place holder for minium difference between y_{i,t-1} vs y_{i^\prime,t-1} *\/ */
	  double tmpydiff;
	  double minydiff;

	  /* finding nearest neighbor  */
	  if (*maxdev < 0) {
	    /* nearest neighbor */
	    int matchiprime;
	    int matchjprime;

	    for (k = 0 ; k < *len_data ; k++) {
	      if (unit_index[k] != i+1 && time_index[k] == j && tr[k] == 0) {
	  	/* k corresponds to the index for the t-1 observation
	  	   of the matched observation */
	  	iprime = unit_index[k];
	  	jprime = time_index[k];

	  	/* observation for time t exists and the treatment status is control */
	  	if (exist[jprime][iprime-1]==1 && (same[jprime][iprime-1]==1)){
	  	  tmpydiff = absolute(y_it1 - y[k]); /* absolute diff */

	  	  if(count_control == 0){ /* first matched observation */
	  	    minydiff = tmpydiff;
	  	    matchiprime = iprime;
	  	    matchjprime = jprime;
	  	    count_control++;
	  	  } else if (count_control > 0 && tmpydiff<= minydiff) {
	  	    minydiff = tmpydiff;
	  	    matchiprime = iprime;
	  	    matchjprime = jprime;
	  	    count_control++;
	  	  } else if (count_control > 0 && tmpydiff > minydiff) {
	  	    ;
	  	  }
		     
	  	}
	      }
	    }

	    /* if match is found */
	    if (count_control > 0 && exist[matchjprime][matchiprime-1]==1 && (same[matchjprime][matchiprime-1]==1)) {
	      W[j][i] = 1.0;
	      W[j-1][i] = 1.0;
	      W[matchjprime][matchiprime-1] = 1.0;
	      W[matchjprime-1][matchiprime-1] = -1.0;

	      /* Rprintf("Unit %d Time %d outcome %f: Nearest Neighbor is Unit %d Time %d difference %f\n", */
	      /* 	       (i+1), (j+1), y_it, (matchiprime), (matchjprime+1), minydiff); */

	    }
	  }

	  else if (*maxdev >=0) {

	  for (k = 0 ; k < *len_data ; k++) {
	    if (unit_index[k] != i+1 && time_index[k] == j && tr[k] == 0) {
	      /* k corresponds to the index for the t-1 observation
		 of the matched observation */
	      iprime = unit_index[k];
	      jprime = time_index[k];
		     
	      if (exist[jprime][iprime-1]==1 && (same[jprime][iprime-1]==1)){
		tmpydiff = absolute(y_it1 - y[k]); /* absolute diff */
		if(tmpydiff <= *maxdev) {
		  count_control++;
		}
	      }
	    }
	  }

	  /* if match is found */
	  if (count_control > 0) {
	    double v_it = 1/count_control;
	    W[j][i] = 1.0;
	    W[j-1][i] = 1.0;

	    for (k = 0 ; k < *len_data ; k++) {
	      if (unit_index[k] != i+1 && time_index[k] == j && tr[k] == 0) {
		iprime = unit_index[k];
		jprime = time_index[k];

		if (exist[jprime][iprime-1]==1 && (same[jprime][iprime-1]==1)){
		  tmpydiff = absolute(y_it1 - y[k]); /* absolute diff */
		  if(tmpydiff <= *maxdev) {
		    W[jprime][iprime-1] =  v_it;
		    W[jprime-1][iprime-1] = - v_it;
		    /* Rprintf("Unit %d Time %d outcome %f: Matched to Unit %d Time %d: the difference is within %f\n", */
		    /* 	   (i+1), (j), y_it1, (iprime), (jprime), *maxdev); */
		  }
		}
	      }
	    }
	  }
	       
	}
	     
	     
      }

      if (t_it == 0) { 
	// count treat
	double count_treat = 0;
	int iprime;
	int jprime;

	/* place holder for minium difference between y_{i,t-1} vs y_{i^\prime,t-1} */
	double tmpydiff;
	double minydiff;

	/* finding nearest neighbor  */
	if (*maxdev < 0) {
	  /* nearest neighbor */
	  int matchiprime;
	  int matchjprime;
	       
	  for (k = 0 ; k < *len_data ; k++) {
	    if (unit_index[k] != i+1 && time_index[k] == j && tr[k] == 1) {
	      iprime = unit_index[k];
	      jprime = time_index[k];

	      if (exist[jprime][iprime-1]==1 && (same[jprime][iprime-1]==1)){
		tmpydiff = absolute(y_it1 - y[k]); /* absolute diff */

		if(count_treat == 0){ /* first matched observation */
		  minydiff = tmpydiff;
		  matchiprime = iprime;
		  matchjprime = jprime;
		  count_treat++;
		} else if (count_treat > 0 && tmpydiff<= minydiff) {
		  minydiff = tmpydiff;
		  matchiprime = iprime;
		  matchjprime = jprime;
		  count_treat++;
		} else if (count_treat > 0 && tmpydiff > minydiff) {
		  ;
		}

	      }
	    }
	  }

	  /* if match is found */
	  if (count_treat > 0 && exist[matchjprime][matchiprime-1]==1 && (same[matchjprime][matchiprime-1]==1)) {
	    W[j][i] = 1.0;
	    W[j-1][i] = 1.0;
	    W[matchjprime][matchiprime-1] = 1.0;
	    W[matchjprime-1][matchiprime-1] = -1.0;
	    /* Rprintf("Unit %d Time %d outcome %f: Nearest Neighbor is Unit %d Time %d difference %f\n", */
	    /* 	       (i+1), (j+1), y_it, (matchiprime-1), (matchjprime+1), minydiff); */
	       
	  }

	}

	else if (*maxdev >=0) {
	       
	  for (k = 0 ; k < *len_data ; k++) {
	    if (unit_index[k] != i+1 && time_index[k] == j && tr[k] == 1) {
	      /* k corresponds to the index for the t-1 observation
		 of the matched observation */
	      iprime = unit_index[k];
	      jprime = time_index[k];
		     
	      if (exist[jprime][iprime-1]==1 && (same[jprime][iprime-1]==1)){
		tmpydiff = absolute(y_it1 - y[k]); /* absolute diff */
		if(tmpydiff<= *maxdev) {
		  count_treat++;
		}
	      }
	    }
	  }

	  /* if match is found */
	  if (count_treat > 0) {
	    double v_it = 1/count_treat;
	    W[j][i] = 1.0;
	    W[j-1][i] = 1.0;

	    for (k = 0 ; k < *len_data ; k++) {
	      if (unit_index[k] != i+1 && time_index[k] == j && tr[k] == 1) {
		iprime = unit_index[k];
		jprime = time_index[k];

		if (exist[jprime][iprime-1]==1 && (same[jprime][iprime-1]==1)){
		  tmpydiff = absolute(y_it1 - y[k]); /* absolute diff */
		  if(tmpydiff <= *maxdev) {
		    W[jprime][iprime-1] =  v_it;
		    W[jprime-1][iprime-1] = - v_it;
		    /* Rprintf("Unit %d Time %d outcome %f: Matched to Unit %d Time %d: the difference is within %f\n", */
		    /* 	   (i+1), (j), y_it1, (iprime), (jprime), *maxdev); */
		  }
		}
	      }
	    }
	  }
	}
      }


	     	
      if (*ate == 1) {
	for (n = 0 ; n < *len_u_index ; n++) {
	  for (m = 0 ; m < *len_t_index ; m++) {
	    Wdid[m][n] = Wdid[m][n] + (c_it*W[m][n]);
	  }
	}
	//Rprintf("Print Wdid %d okay?\n", 1);
	//PdoubleMatrix(Wdid, *len_t_index, *len_u_index);     
      } 
      else if (*att == 1) {
	for (n = 0 ; n < *len_u_index ; n++) {
	  for (m = 0 ; m < *len_t_index ; m++) {
	    Wdid[m][n] = Wdid[m][n] + (c_it*t_it*W[m][n]);
	  }
	}
	//Rprintf("Print Wdid %d okay?\n", 1);
	//PdoubleMatrix(Wdid, *len_t_index, *len_u_index);    
      }
	   
    }
	 
  }
      
} 
   
  
/* PdoubleMatrix(Wdid, *len_t_index, *len_u_index); */
   
iter = 0;
for (m = 0; m < *len_t_index; m++) {
  for (n = 0; n < *len_u_index; n++) {
    weightdid[iter] = Wdid[m][n];
    iter++;
  }
 }  
   
FreeMatrix(Wdid, *len_t_index);
FreeMatrix(W, *len_t_index);
FreeintMatrix(exist, *len_t_index);  
FreeintMatrix(same, *len_t_index);


}





/* Generate Dummy Variables */
 
void MDummy(int* index, int* len_index, int* len_data, int* dummy) { 

  int **DMatrix = intMatrix(*len_data, *len_index);
  int i;
  int k;
  int iter = 0;
  
  for (i = 0 ; i < *len_index ; i++) {
    for (k = 0 ; k < *len_data ; k++) {
      if ( index[k] == (i+1) ) {
	DMatrix[k][i] = 1;
      }
      else {
	DMatrix[k][i] = 0;
      }
    }
  }
  
  for (i = 0; i < *len_index; i++) {
    for (k = 0; k < *len_data; k++) {
      dummy[iter] = DMatrix[k][i];
      iter++;
    }
  }
  
  FreeintMatrix(DMatrix, *len_data);
}

   
/* Function for summing over i t(X_i)%*% (X_i) */

void XXiSum(int *len_data, int *n_cov,
	    int *unit_index, int *len_uniq_u_index,
	    double *Xtilde,
	    double *result) {


  /* temporary storages */
  int i, j, k, itemp;
  int n;
  
  int count;
  double **X;
  double **XX;

  /* defining vectors and matricies */
  X = doubleMatrix(*len_data, *n_cov);
  Rprintf("x allocated\n");
  XX = doubleMatrix(*n_cov, *n_cov);
  Rprintf("XX allocated\n");

  /* packing X */
  itemp = 0;
  for (j = 0; j < *n_cov; j++) {
    for (i = 0; i < *len_data; i++) {
      X[i][j] = Xtilde[itemp++];
    }
  }

  /* packing XX */
  for (i = 0; i < *n_cov; i++) {
    for (j = 0; j < *n_cov; j++) {
      XX[i][j] = 0;
    }
  }

  /* calculating t(X.tilda) %*% X.tilde */


  for (n = 0; n < *len_uniq_u_index; n++) {
    count = 0;
    /* Rprintf("unit is  %d\n", (n+1)); */
    /* calculating t(X.tilda.i) %*% X.tilde.i */

    for (k = 0; k < *len_data; k++) {
      if (unit_index[k] == (n+1)) {
	count ++;
      }
    }
   

        
    double **X_i;
    double **tX_i;
    double **XX_i;

    X_i = doubleMatrix(count, *n_cov);
    Rprintf("x_i allocated\n");
    tX_i = doubleMatrix(*n_cov, count);
    Rprintf("tx_i allocated\n");
    XX_i = doubleMatrix(*n_cov, *n_cov);
    Rprintf("xx_i allocated\n");
  
    /* initializing XX_i */
    for(i=0; i < *n_cov; i++) {
      for(j=0; j < *n_cov; j++){
	XX_i[i][j]=0;
      }
    }
      
    /* packing X_i */
    itemp = 0;
    for (k = 0; k < *len_data; k++) {
      if (unit_index[k] == (n+1)) {
	for (j = 0; j < *n_cov; j++) {
	  X_i[itemp][j] = X[k][j];
	}
	itemp++;
      }
    }
   

    /* packing tX_i */

    for (j = 0; j < *n_cov; j++) {
      for (k = 0; k < count; k++) {
	tX_i[j][k] = X_i[k][j];
      }
    }


    /* Matrix multiplication */
    for(i=0; i < *n_cov; i++) {
      for(j=0; j < *n_cov; j++){
	for (k = 0; k < count; k++) {
	  XX_i[i][j] += tX_i[i][k] * X_i[k][j];
	}
      }
    }

    /* Adding results in XX */
    for(i=0; i < *n_cov; i++) {
      for(j=0; j < *n_cov; j++){
	XX[i][j] +=  XX_i[i][j];
      }
    }
      
      
    FreeMatrix(X_i, count);
    FreeMatrix(tX_i, *n_cov);
    FreeMatrix(XX_i, *n_cov);

      
  }

  /* Rprintf("XX matrix in C \n"); */
  /* PdoubleMatrix(XX, *n_cov, *n_cov); */


  itemp = 0;
  for (i = 0; i < *n_cov; i++) {
    for (j = 0; j < *n_cov; j++) {
      result[itemp++] = XX[i][j];
    }
  }

  FreeMatrix(X, *len_data);
  FreeMatrix(XX, *n_cov);

   
}





/* Function for summing over i W_i t(X_i)%*% (X_i) */

void XWXiSum(int *len_data, int *n_cov,
	     int *unit_index, int *len_uniq_u_index,
	     double *Xtilde,
	     double *weights,
	     double *result) {


  /* temporary storages */
  int i, j, k, itemp;
  int n;
  int count;
   
  double **X;
  double *w;
  double **XWX;

  /* defining vectors and matricies */
  X = doubleMatrix(*len_data, *n_cov);
  w = doubleArray(*len_data);
  XWX = doubleMatrix(*n_cov, *n_cov);

  /* packing X */
  itemp = 0;
  for (j = 0; j < *n_cov; j++) {
    for (i = 0; i < *len_data; i++) {
      X[i][j] = Xtilde[itemp++];
    }
  }
   
  /* PdoubleMatrix(X, *len_data, *n_cov); */
      
   
  /* packing w */
  for (i = 0; i < *len_data; i++) {
    w[i] =  weights[i];
  }

  /* packing XWX */
  for (j = 0; j < *n_cov; j++) {
    for (i = 0; i < *n_cov; i++) {
      XWX[i][j] = 0;
    }
  }
   


  for (n = 0; n < *len_uniq_u_index; n++) {
      
    /* Rprintf("unit %d\n", n+1); */

    count = 0;
   
    for (k = 0; k < *len_data; k++) {
      if (unit_index[k] == (n+1)) {
	count ++;
      }
    }
   

        
    double **X_i;
    double **tX_i;
    double **XW_i;
    double **XWX_i;
    double *w_i;
    double **W_i;

    X_i = doubleMatrix(count, *n_cov);
    tX_i = doubleMatrix(*n_cov, count);
    XW_i = doubleMatrix(*n_cov, count);
    XWX_i = doubleMatrix(*n_cov, *n_cov);
    w_i = doubleArray(count);
    W_i = doubleMatrix(count, count);

      
    /* packing X_i */
    itemp = 0;
    for (k = 0; k < *len_data; k++) {
      if (unit_index[k] == (n+1)) {
	for (j = 0; j < *n_cov; j++) {
	  X_i[itemp][j] = X[k][j];
	}
	itemp++;
      }
    }
   
    /* PdoubleMatrix(X_i, count, *n_cov); */



    /* packing tX_i */

    for (j = 0; j < *n_cov; j++) {
      for (k = 0; k < count; k++) {
	tX_i[j][k] = X_i[k][j];
      }
    }


 
    /* packing w_i */
    itemp = 0;
    for (k = 0; k < *len_data; k++) {
      if (unit_index[k] == (n+1)) {
	w_i[itemp] = w[k];
	itemp++;
      }
    }
     

       
    /* packing W_i */
    for(i=0; i < count; i++) {
      for(j=0; j < count; j++){
	if (i == j){
	  W_i[i][j] = w_i[i];
	} else {
	  W_i[i][j] = 0;
	}
      }
    }

        

    /* packing XW_i */
    for(i=0; i < *n_cov; i++) {
      for(j=0; j < count; j++){
	XW_i[i][j] = 0;
	for (k = 0; k < count; k++) {
	  XW_i[i][j] += tX_i[i][k] * W_i[k][j];
	}
      }
    }
  

    /* packing XWX_i */
    for(i=0; i < *n_cov; i++) {
      for(j=0; j < *n_cov; j++){
	XWX_i[i][j] = 0;
	for (k = 0; k < count; k++) {
	  XWX_i[i][j] += XW_i[i][k] * X_i[k][j];
	}
      }
    }
  
    /* PdoubleMatrix(XWX_i, *n_cov, *n_cov); */



    /* Adding results in XWX */
    for(i=0; i < *n_cov; i++) {
      for(j=0; j < *n_cov; j++){
	XWX[i][j] +=  XWX_i[i][j];
      }
    }
      
      
    FreeMatrix(X_i, count);
    FreeMatrix(tX_i, *n_cov);
    FreeMatrix(XW_i, *n_cov);
    FreeMatrix(XWX_i, count);
    free(w_i);
    FreeMatrix(W_i, count);
    
  


  
  }



  itemp = 0;
  for (i = 0; i < *n_cov; i++) {
    for (j = 0; j < *n_cov; j++) {
      result[itemp++] = XWX[i][j];
    }
  }

  FreeMatrix(X, *len_data);
  free(w);
  FreeMatrix(XWX, *n_cov);
      
}






void OmegaHatHAC(int* len_data, int* n_cov,
		 int* unit_index, int* len_uniq_u_index,
		 double* Xtilde,
		 double* utilde,
		 double* Omega_hat_HAC) {


  /* temporary storages */
  int i, j, k, itemp;
  int n;
  int count;

  double **X;
  double **OmegaMatrix;
  double *u;

  /* defining vectors and matricies */
  X = doubleMatrix(*len_data, *n_cov);
  OmegaMatrix = doubleMatrix(*n_cov, *n_cov);
  u = doubleArray(*len_data);


  /* packing X */
  itemp = 0;
  for (j = 0; j < *n_cov; j++) {
    for (i = 0; i < *len_data; i++) {
      X[i][j] = Xtilde[itemp++];
    }
  }
 
  /* packing OmegaMatrix */
  itemp = 0;
  for (j = 0; j < *n_cov; j++) {
    for (i = 0; i < *n_cov; i++) {
      OmegaMatrix[i][j] = 0;
    }
  }
 
  /* packing u */
  for (i = 0; i < *len_data; i++) {
    u[i] = utilde[i];
  }

   
  /* PdoubleArray(u, *len_data); */


  /* calculating t(X.tilda.i) %*% u.tilde.i %*% t(u.tilde.i) %*% X.tilde.i */
   
  /* Rprintf("unit length is  %d\n", *len_uniq_u_index); */

  for (n = 0; n < *len_uniq_u_index; n++) {
    count = 0;
    /* Rprintf("unit is  %d\n", (n+1)); */
    
    /* calculating t(X.tilda.i) %*% X.tilde.i */


    for (k = 0; k < *len_data; k++) {
      if (unit_index[k] == (n+1)) {
	count ++;
      }
    }
   

    /* Rprintf("count is  %d\n", count);  */
      
    double **X_i;
    double **tX_i;
    double *u_i;
    double **uu_i;
    double **Xuu_i;
    double **XuuX_i;


    X_i = doubleMatrix(count, *n_cov);
    tX_i = doubleMatrix(*n_cov, count);
    u_i = doubleArray(count);
    uu_i = doubleMatrix(count, count);
    Xuu_i = doubleMatrix(*n_cov, count);
    XuuX_i = doubleMatrix(*n_cov, *n_cov);

      
    /* packing X_i */
    itemp = 0;
    for (k = 0; k < *len_data; k++) {
      if (unit_index[k] == (n+1)) {
	for (j = 0; j < *n_cov; j++) {
	  X_i[itemp][j] = X[k][j];
	}
	itemp++;
      }
    }
   
    /* Rprintf("X_i for unit %d is \n", (n+1)); */
    /* PdoubleMatrix(X_i, count, *n_cov);  */
      
    /* Rprintf("number of observations for unit %d is %d\n", (n+1), count ); */

    /* packing tX_i */

    for (j = 0; j < *n_cov; j++) {
      for (k = 0; k < count; k++) {
	tX_i[j][k] = X_i[k][j];
      }
    }

    /* Rprintf("t(X) for unit %d is \n", (n+1)); */
    /* PdoubleMatrix(tX_i, *n_cov, count); */
 
    /* packing u_i */
    itemp = 0;
    for (k = 0; k < *len_data; k++) {
      if (unit_index[k] == (n+1)) {
	u_i[itemp] = u[k];
	itemp++;
      }
    }
   
    /* Rprintf("u_i for unit %d is \n", (n+1)); */
    /* PdoubleArray(u_i, count); */
   
    /* packing uu_i */

    for(i=0; i < count; i++) {
      for(j=0; j < count; j++){
	for (k = 0; k < count; k++) {
	  uu_i[i][j] = u_i[i] * u_i[j];
	}
      }
    }

    /* Rprintf("uu_i for unit %d is \n", (n+1)); */
    /* PdoubleMatrix(uu_i, count, count); */
         

    /* Matrix multiplication */


    for(i=0; i < *n_cov; i++) {
      for(j=0; j < count; j++){
	Xuu_i[i][j] = 0;
	for (k = 0; k < count; k++) {
	  Xuu_i[i][j] += tX_i[i][k] * uu_i[k][j];
	}
      }
    }
  
    /* Rprintf("Xuu_i for unit %d is \n", (n+1)); */
    /* PdoubleMatrix(Xuu_i, *n_cov, count); */
    


    for(i=0; i < *n_cov; i++) {
      for(j=0; j < *n_cov; j++){
	XuuX_i[i][j] = 0;
	for (k = 0; k < count; k++) {
	  XuuX_i[i][j] += Xuu_i[i][k] * X_i[k][j];
	}
      }
    }
  
    /* Rprintf("XuuX_i for unit %d is \n", (n+1)); */
    /* PdoubleMatrix(XuuX_i, *n_cov, *n_cov); */
    

    /* Adding results in OmegaMatrix  */
    for(i=0; i < *n_cov; i++) {
      for(j=0; j < *n_cov; j++){
	OmegaMatrix[i][j] +=  XuuX_i[i][j];
      }
    }



    FreeMatrix(X_i, count);
    FreeMatrix(tX_i, *n_cov);
    free(u_i);
    FreeMatrix(uu_i, count);
    FreeMatrix(Xuu_i, *n_cov);
    FreeMatrix(XuuX_i, *n_cov);


      
  }
   
  /* storing results */

  itemp = 0;
  for (i = 0; i < *n_cov; i++) {
    for (j = 0; j < *n_cov; j++) {
      Omega_hat_HAC[itemp++] = OmegaMatrix[i][j];
    }
  }

  FreeMatrix(X, *len_data);
  FreeMatrix(OmegaMatrix, *n_cov);
  free(u);

}





/* Heteroskedasticity constitent Omega */


void OmegaHatHC(int* len_data, int* n_cov,
		int* unit_index, int* len_uniq_u_index,
		double* Xtilde,
		double* utilde,
		double* Omega_hat_HC) {


  /* temporary storages */
  int i, j, k, itemp;
  int n;
  int count;

  double **X;
  double **OmegaMatrix;
  double *u;

  /* defining vectors and matricies */
  X = doubleMatrix(*len_data, *n_cov);
  OmegaMatrix = doubleMatrix(*n_cov, *n_cov);
  u = doubleArray(*len_data);


  /* packing X */
  itemp = 0;
  for (j = 0; j < *n_cov; j++) {
    for (i = 0; i < *len_data; i++) {
      X[i][j] = Xtilde[itemp++];
    }
  }
 
  /* packing OmegaMatrix */
  itemp = 0;
  for (j = 0; j < *n_cov; j++) {
    for (i = 0; i < *n_cov; i++) {
      OmegaMatrix[i][j] = 0;
    }
  }
 
  /* packing u */
  for (i = 0; i < *len_data; i++) {
    u[i] = utilde[i];
  }

   
  /* PdoubleArray(u, *len_data); */


  /* calculating t(X.tilda.i) %*% u.tilde.i %*% t(u.tilde.i) %*% X.tilde.i */
   
  /* Rprintf("unit length is  %d\n", *len_uniq_u_index); */

  for (n = 0; n < *len_uniq_u_index; n++) {
    count = 0;
    /* Rprintf("unit is  %d\n", (n+1)); */
    
    /* calculating t(X.tilda.i) %*% X.tilde.i */

    for (k = 0; k < *len_data; k++) {
      if (unit_index[k] == (n+1)) {
	count ++;
      }
    }
   

    /* Rprintf("count is  %d\n", count);  */
      
    double **X_i;
    double **tX_i;
    double *u_i;
    double **uu_i;
    double **Xuu_i;
    double **XuuX_i;
    double *u_i2;


    X_i = doubleMatrix(count, *n_cov);
    tX_i = doubleMatrix(*n_cov, count);
    u_i = doubleArray(count);
    uu_i = doubleMatrix(count, count);
    Xuu_i = doubleMatrix(*n_cov, count);
    XuuX_i = doubleMatrix(*n_cov, *n_cov);
    u_i2 = doubleArray(count);
      
    /* packing X_i */
    itemp = 0;
    for (k = 0; k < *len_data; k++) {
      if (unit_index[k] == (n+1)) {
	for (j = 0; j < *n_cov; j++) {
	  X_i[itemp][j] = X[k][j];
	}
	itemp++;
      }
    }
   
    /* Rprintf("X_i for unit %d is \n", (n+1)); */
    /* PdoubleMatrix(X_i, count, *n_cov);  */
      
    /* Rprintf("number of observations for unit %d is %d\n", (n+1), count ); */

    /* packing tX_i */

    for (j = 0; j < *n_cov; j++) {
      for (k = 0; k < count; k++) {
	tX_i[j][k] = X_i[k][j];
      }
    }

    /* Rprintf("t(X) for unit %d is \n", (n+1)); */
    /* PdoubleMatrix(tX_i, *n_cov, count); */
 
    /* packing u_i */
    itemp = 0;
    for (k = 0; k < *len_data; k++) {
      if (unit_index[k] == (n+1)) {
	u_i[itemp] = u[k];
	itemp++;
      }
    }
   
    /* Rprintf("u_i for unit %d is \n", (n+1)); */
    /* PdoubleArray(u_i, count); */
   

    /* packing u_i2 */

    for (k = 0; k < count; k++) {
      u_i2[k] = u_i[k] * u_i[k];
    }
   



    /* packing uu_i */

    for(i=0; i < count; i++) {
      for(j=0; j < count; j++){
	if (i == j){
	  uu_i[i][j] = u_i2[i];
	} else {
	  uu_i[i][j] = 0;
	}
      }
    }

    /* Rprintf("uu_i for unit %d is \n", (n+1)); */
    /* PdoubleMatrix(uu_i, count, count); */
         

    /* Matrix multiplication */


    for(i=0; i < *n_cov; i++) {
      for(j=0; j < count; j++){
	Xuu_i[i][j] = 0;
	for (k = 0; k < count; k++) {
	  Xuu_i[i][j] += tX_i[i][k] * uu_i[k][j];
	}
      }
    }
  
    /* Rprintf("Xuu_i for unit %d is \n", (n+1)); */
    /* PdoubleMatrix(Xuu_i, *n_cov, count); */
    


    for(i=0; i < *n_cov; i++) {
      for(j=0; j < *n_cov; j++){
	XuuX_i[i][j] = 0;
	for (k = 0; k < count; k++) {
	  XuuX_i[i][j] += Xuu_i[i][k] * X_i[k][j];
	}
      }
    }
  
    /* Rprintf("XuuX_i for unit %d is \n", (n+1)); */
    /* PdoubleMatrix(XuuX_i, *n_cov, *n_cov); */
    

    /* Adding results in OmegaMatrix  */
    for(i=0; i < *n_cov; i++) {
      for(j=0; j < *n_cov; j++){
	OmegaMatrix[i][j] +=  XuuX_i[i][j];
      }
    }



    FreeMatrix(X_i, count);
    FreeMatrix(tX_i, *n_cov);
    free(u_i);
    FreeMatrix(uu_i, count);
    FreeMatrix(Xuu_i, *n_cov);
    FreeMatrix(XuuX_i, *n_cov);
    free(u_i2);


      
  }
   
  /* storing results */

  itemp = 0;
  for (i = 0; i < *n_cov; i++) {
    for (j = 0; j < *n_cov; j++) {
      Omega_hat_HC[itemp++] = OmegaMatrix[i][j];
    }
  }

  FreeMatrix(X, *len_data);
  FreeMatrix(OmegaMatrix, *n_cov);
  free(u);

}





void OmegaDiDHAC(int* len_data, int* n_cov,
		 int* unit_index, int* len_uniq_u_index,
		 double* Xtilde,
		 double* utilde,
		 double* W,
		 double* Omega_DiD_HAC) {


  /* temporary storages */
  int i, j, k, itemp;
  int n;
  int count;

  double **X;
  double **OmegaMatrix;
  double *u;
  double *w;

  /* defining vectors and matricies */
  X = doubleMatrix(*len_data, *n_cov);
  OmegaMatrix = doubleMatrix(*n_cov, *n_cov);
  u = doubleArray(*len_data);
  w = doubleArray(*len_data);

  /* packing X */
  itemp = 0;
  for (j = 0; j < *n_cov; j++) {
    for (i = 0; i < *len_data; i++) {
      X[i][j] = Xtilde[itemp++];
    }
  }
 
  /* packing OmegaMatrix */
  itemp = 0;
  for (j = 0; j < *n_cov; j++) {
    for (i = 0; i < *n_cov; i++) {
      OmegaMatrix[i][j] = 0;
    }
  }
 
  /* packing u */
  for (i = 0; i < *len_data; i++) {
    u[i] = utilde[i];
  }

  /* packing w */
  for (i = 0; i < *len_data; i++) {
    w[i] = W[i];
  }
   

  /* PdoubleArray(u, *len_data); */


  /* calculating t(X.tilda.i) %*% u.tilde.i %*% t(u.tilde.i) %*% X.tilde.i */
   
  /* Rprintf("unit length is  %d\n", *len_uniq_u_index); */

  for (n = 0; n < *len_uniq_u_index; n++) {
    count = 0;
    /* Rprintf("unit is  %d\n", (n+1)); */
    
    /* calculating t(X.tilda.i) %*% X.tilde.i */

    for (k = 0; k < *len_data; k++) {
      if (unit_index[k] == (n+1)) {
	count ++;
      }
    }
   

    /* Rprintf("count is  %d\n", count);  */
      
    double **X_i;
    double **tX_i;
    double *u_i;
    double **uu_i;
    double **XWuu_i;
    double **XWuuWX_i;
    double *w_i;
    double **W_i;
    double **XW_i;
    double **WX_i;


    X_i = doubleMatrix(count, *n_cov);
    tX_i = doubleMatrix(*n_cov, count);
    u_i = doubleArray(count);
    uu_i = doubleMatrix(count, count);
    XWuu_i = doubleMatrix(*n_cov, count);
    XWuuWX_i = doubleMatrix(*n_cov, *n_cov);
    w_i = doubleArray(count);
    W_i = doubleMatrix(count, count);
    XW_i = doubleMatrix(*n_cov, count);
    WX_i = doubleMatrix(count, *n_cov);
      
      
    /* packing X_i */
    itemp = 0;
    for (k = 0; k < *len_data; k++) {
      if (unit_index[k] == (n+1)) {
	for (j = 0; j < *n_cov; j++) {
	  X_i[itemp][j] = X[k][j];
	}
	itemp++;
      }
    }
   
    /* Rprintf("X_i for unit %d is \n", (n+1)); */
    /* PdoubleMatrix(X_i, count, *n_cov);  */
      
    /* Rprintf("number of observations for unit %d is %d\n", (n+1), count ); */

    /* packing tX_i */

    for (j = 0; j < *n_cov; j++) {
      for (k = 0; k < count; k++) {
	tX_i[j][k] = X_i[k][j];
      }
    }

    /* Rprintf("t(X) for unit %d is \n", (n+1)); */
    /* PdoubleMatrix(tX_i, *n_cov, count); */
 
    /* packing u_i */
    itemp = 0;
    for (k = 0; k < *len_data; k++) {
      if (unit_index[k] == (n+1)) {
	u_i[itemp] = u[k];
	itemp++;
      }
    }
 
    /* packing w_i */
    itemp = 0;
    for (k = 0; k < *len_data; k++) {
      if (unit_index[k] == (n+1)) {
	w_i[itemp] = w[k];
	itemp++;
      }
    }
     
    /* Rprintf("u_i for unit %d is \n", (n+1)); */
    /* PdoubleArray(u_i, count); */
       
    /* packing W_i */
    for(i=0; i < count; i++) {
      for(j=0; j < count; j++){
	if (i == j){
	  W_i[i][j] = w_i[i];
	} else {
	  W_i[i][j] = 0;
	}
      }
    }


    /* packing XW_i */
    for(i=0; i < *n_cov; i++) {
      for(j=0; j < count; j++){
	XW_i[i][j] = 0;
	for (k = 0; k < count; k++) {
	  XW_i[i][j] += tX_i[i][k] * W_i[k][j];
	}
      }
    }
  
    /* packing WX_i */
    for(i=0; i < count; i++) {
      for(j=0; j < *n_cov; j++){
	WX_i[i][j] = 0;
	for (k = 0; k < count; k++) {
	  WX_i[i][j] += W_i[i][k] * X_i[k][j]; 
	}
      }
    }
  


    /* packing uu_i */

    for(i=0; i < count; i++) {
      for(j=0; j < count; j++){
	uu_i[i][j] = u_i[i] * u_i[j];
      }
    }
      

 
    /* Rprintf("uu_i for unit %d is \n", (n+1)); */
    /* PdoubleMatrix(uu_i, count, count); */
         

    /* Matrix multiplication */


    for(i=0; i < *n_cov; i++) {
      for(j=0; j < count; j++){
	XWuu_i[i][j] = 0;
	for (k = 0; k < count; k++) {
	  XWuu_i[i][j] += XW_i[i][k] * uu_i[k][j];
	}
      }
    }
  
    /* Rprintf("Xuu_i for unit %d is \n", (n+1)); */
    /* PdoubleMatrix(Xuu_i, *n_cov, count); */
    


    for(i=0; i < *n_cov; i++) {
      for(j=0; j < *n_cov; j++){
	XWuuWX_i[i][j] = 0;
	for (k = 0; k < count; k++) {
	  XWuuWX_i[i][j] += XWuu_i[i][k] * WX_i[k][j];
	}
      }
    }
  
    /* Rprintf("XWuuWX_i for unit %d is \n", (n+1)); */
    /* PdoubleMatrix(XWuuWX_i, *n_cov, *n_cov); */
    

    /* Adding results in OmegaMatrix  */
    for(i=0; i < *n_cov; i++) {
      for(j=0; j < *n_cov; j++){
	OmegaMatrix[i][j] +=  XWuuWX_i[i][j];
      }
    }



    FreeMatrix(X_i, count);
    FreeMatrix(tX_i, *n_cov);
    free(u_i);
    FreeMatrix(uu_i, count);
    FreeMatrix(XWuu_i, *n_cov);
    FreeMatrix(XWuuWX_i, *n_cov);
    free(w_i);
    FreeMatrix(W_i, count);
    FreeMatrix(XW_i, *n_cov);
    FreeMatrix(WX_i, count);



  }
   
  /* storing results */

  itemp = 0;
  for (i = 0; i < *n_cov; i++) {
    for (j = 0; j < *n_cov; j++) {
      Omega_DiD_HAC[itemp++] = OmegaMatrix[i][j];
    }
  }

  FreeMatrix(X, *len_data);
  FreeMatrix(OmegaMatrix, *n_cov);
  free(u);
  free(w);

}






void OmegaDiDHAC2(int* len_data, int* n_cov,
		  int* unit_index, int* len_uniq_u_index,
		  double* Xtilde,
		  double* utilde,
		  double* W,
		  double* Omega_DiD_HAC) {
   
   
  /* temporary storages */
  int i, j, k, itemp;
  int n;
  int count;

  double **X;
  double **OmegaMatrix;
  double *u;
  double *w;

  /* defining vectors and matricies */
  X = doubleMatrix(*len_data, *n_cov);
  OmegaMatrix = doubleMatrix(*n_cov, *n_cov);
  u = doubleArray(*len_data);
  w = doubleArray(*len_data);


  /* packing X */
  itemp = 0;
  for (j = 0; j < *n_cov; j++) {
    for (i = 0; i < *len_data; i++) {
      X[i][j] = Xtilde[itemp++];
    }
  }
 
  /* packing OmegaMatrix */
  itemp = 0;
  for (j = 0; j < *n_cov; j++) {
    for (i = 0; i < *n_cov; i++) {
      OmegaMatrix[i][j] = 0;
    }
  }
 
  /* packing u */
  for (i = 0; i < *len_data; i++) {
    u[i] = utilde[i];
  }

  /* packing w */
  for (i = 0; i < *len_data; i++) {
    w[i] = W[i];
  }
   

  /* PdoubleArray(u, *len_data); */


  /* calculating t(X.tilda.i) %*% W_i %*%  u.tilde.i %*% t(u.tilde.i) %*% W_i %*% X.tilde.i */
  /* Rprintf("unit length is  %d\n", *len_uniq_u_index); */
  Rprintf("done up to here\n");

     
  for (n = 0; n < *len_uniq_u_index; n++) {
    count = 0;
    /* Rprintf("unit is  %d\n", (n+1)); */
    /* calculating t(X.tilda.i) %*% X.tilde.i */

    for (k = 0; k < *len_data; k++) {
      if (unit_index[k] == (n+1)) {
	count ++;
      }
    }

    /* Rprintf("count is  %d\n", count);  */
   
    double **X_i;
    double *w_i;
    double *u_i;
    double *wu_i;
    double **Xwu_i;
    double **XwuuwX_i;


    X_i = doubleMatrix(count, *n_cov);
    w_i = doubleArray(count);
    u_i = doubleArray(count);
    wu_i = doubleArray(count);
    Xwu_i = doubleMatrix(*n_cov, count);
    XwuuwX_i = doubleMatrix(*n_cov, *n_cov);
      
      
    /* packing X_i , w_i, and u_i corresponding each unit i*/
    itemp = 0;
    for (k = 0; k < *len_data; k++) {
      if (unit_index[k] == (n+1)) {
	for (j = 0; j < *n_cov; j++) {
	  X_i[itemp][j] = X[k][j];
	  w_i[itemp] = w[k];
	  u_i[itemp] = u[k];
	}
	itemp++;
      }
    }


   
       
    /* packing wu_i */
    for(i=0; i < count; i++) {
      wu_i[i] = w_i[i] * u_i[i];
    }
     
    /* initialize Xwu_i */
    for(i=0; i < *n_cov; i++) {
      for(j=0; j < count; j++){
	Xwu_i[i][j] = 0;
      }
    }
    /* packing Xwu_i */
    for(i=0; i < *n_cov; i++) {
      for(j=0; j < count; j++){
	/* packing Xwu_i */
	Xwu_i[i][j] = X_i[j][i] * wu_i[j];
      }
    }

    PdoubleMatrix(Xwu_i, *n_cov, count);
    

    /* packing XwuuwX_i */
    for(i=0; i < *n_cov; i++) {
      for(j=0; j < *n_cov; j++){
	XwuuwX_i[i][j] = 0;
      }
    }
    
    for(i=0; i < *n_cov; i++) {
      for(j=0; j < *n_cov; j++){
	for (k = 0; k < count; k++) {
	  XwuuwX_i[i][j] += Xwu_i[i][k] * Xwu_i[j][k];
	}
      }
    }
    
    /* Rprintf("XWuuWX_i for unit %d is \n", (n+1)); */
    /* PdoubleMatrix(XWuuWX_i, *n_cov, *n_cov); */
    

    /* Adding results in OmegaMatrix  */
    for(i=0; i < *n_cov; i++) {
      for(j=0; j < *n_cov; j++){
	OmegaMatrix[i][j] +=  XwuuwX_i[i][j];
      }
    }



    FreeMatrix(X_i, count);
    free(w_i);
    free(u_i);
    free(wu_i);
    FreeMatrix(Xwu_i, *n_cov);
    FreeMatrix(XwuuwX_i, *n_cov);


  }
   
  /* /\* storing results *\/ */

  /* itemp = 0; */
  /* for (i = 0; i < *n_cov; i++) { */
  /*   for (j = 0; j < *n_cov; j++) { */
  /*     Omega_DiD_HAC[itemp++] = OmegaMatrix[i][j]; */
  /*   } */
  /* } */

  FreeMatrix(X, *len_data);
  FreeMatrix(OmegaMatrix, *n_cov);
  free(u);
  free(w);

}









/* Heteroskedasticity constitent Omega */


void OmegaDiDHC(int* len_data, int* n_cov,
		int* unit_index, int* len_uniq_u_index,
		double* Xtilde,
		double* utilde,
		double* W,
		double* Omega_DiD_HC) {


  /* temporary storages */
  int i, j, k, itemp;
  int n;
  int count;

  double **X;
  double **OmegaMatrix;
  double *u;
  double *w;

  /* defining vectors and matricies */
  X = doubleMatrix(*len_data, *n_cov);
  OmegaMatrix = doubleMatrix(*n_cov, *n_cov);
  u = doubleArray(*len_data);
  w = doubleArray(*len_data);

  /* packing X */
  itemp = 0;
  for (j = 0; j < *n_cov; j++) {
    for (i = 0; i < *len_data; i++) {
      X[i][j] = Xtilde[itemp++];
    }
  }
 
  /* packing OmegaMatrix */
  itemp = 0;
  for (j = 0; j < *n_cov; j++) {
    for (i = 0; i < *n_cov; i++) {
      OmegaMatrix[i][j] = 0;
    }
  }
 
  /* packing u */
  for (i = 0; i < *len_data; i++) {
    u[i] = utilde[i];
  }

  /* packing w */
  for (i = 0; i < *len_data; i++) {
    w[i] = W[i];
  }
   

  /* PdoubleArray(u, *len_data); */


  /* calculating t(X.tilda.i) %*% u.tilde.i %*% t(u.tilde.i) %*% X.tilde.i */
   
  /* Rprintf("unit length is  %d\n", *len_uniq_u_index); */

  for (n = 0; n < *len_uniq_u_index; n++) {
    count = 0;
    /* Rprintf("unit is  %d\n", (n+1)); */
    
    /* calculating t(X.tilda.i) %*% X.tilde.i */

    for (k = 0; k < *len_data; k++) {
      if (unit_index[k] == (n+1)) {
	count ++;
      }
    }
   

    /* Rprintf("count is  %d\n", count);  */
      
    double **X_i;
    double **tX_i;
    double *u_i;
    double **uu_i;
    double **XWuu_i;
    double **XWuuWX_i;
    double *u_i2;
    double *w_i;
    double **W_i;
    double **XW_i;
    double **WX_i;


    X_i = doubleMatrix(count, *n_cov);
    tX_i = doubleMatrix(*n_cov, count);
    u_i = doubleArray(count);
    uu_i = doubleMatrix(count, count);
    XWuu_i = doubleMatrix(*n_cov, count);
    XWuuWX_i = doubleMatrix(*n_cov, *n_cov);
    u_i2 = doubleArray(count);
    w_i = doubleArray(count);
    W_i = doubleMatrix(count, count);
    XW_i = doubleMatrix(*n_cov, count);
    WX_i = doubleMatrix(count, *n_cov);
      
      
    /* packing X_i */
    itemp = 0;
    for (k = 0; k < *len_data; k++) {
      if (unit_index[k] == (n+1)) {
	for (j = 0; j < *n_cov; j++) {
	  X_i[itemp][j] = X[k][j];
	}
	itemp++;
      }
    }
   
    /* Rprintf("X_i for unit %d is \n", (n+1)); */
    /* PdoubleMatrix(X_i, count, *n_cov);  */
      
    /* Rprintf("number of observations for unit %d is %d\n", (n+1), count ); */

    /* packing tX_i */

    for (j = 0; j < *n_cov; j++) {
      for (k = 0; k < count; k++) {
	tX_i[j][k] = X_i[k][j];
      }
    }

    /* Rprintf("t(X) for unit %d is \n", (n+1)); */
    /* PdoubleMatrix(tX_i, *n_cov, count); */
 
    /* packing u_i */
    itemp = 0;
    for (k = 0; k < *len_data; k++) {
      if (unit_index[k] == (n+1)) {
	u_i[itemp] = u[k];
	itemp++;
      }
    }
 
    /* packing w_i */
    itemp = 0;
    for (k = 0; k < *len_data; k++) {
      if (unit_index[k] == (n+1)) {
	w_i[itemp] = w[k];
	itemp++;
      }
    }
     
    /* Rprintf("u_i for unit %d is \n", (n+1)); */
    /* PdoubleArray(u_i, count); */
       
    /* packing W_i */
    for(i=0; i < count; i++) {
      for(j=0; j < count; j++){
	if (i == j){
	  W_i[i][j] = w_i[i];
	} else {
	  W_i[i][j] = 0;
	}
      }
    }


    /* packing XW_i */
    for(i=0; i < *n_cov; i++) {
      for(j=0; j < count; j++){
	XW_i[i][j] = 0;
	for (k = 0; k < count; k++) {
	  XW_i[i][j] += tX_i[i][k] * W_i[k][j];
	}
      }
    }
  
    /* packing WX_i */
    for(i=0; i < count; i++) {
      for(j=0; j < *n_cov; j++){
	WX_i[i][j] = 0;
	for (k = 0; k < count; k++) {
	  WX_i[i][j] += W_i[i][k] * X_i[k][j]; 
	}
      }
    }
  


    /* packing u_i2 */

    for (k = 0; k < count; k++) {
      u_i2[k] = u_i[k] * u_i[k];
    }
   



    /* packing uu_i */

    for(i=0; i < count; i++) {
      for(j=0; j < count; j++){
	if (i == j){
	  uu_i[i][j] = u_i2[i];
	} else {
	  uu_i[i][j] = 0;
	}
      }
    }

    /* Rprintf("uu_i for unit %d is \n", (n+1)); */
    /* PdoubleMatrix(uu_i, count, count); */
         

    /* Matrix multiplication */


    for(i=0; i < *n_cov; i++) {
      for(j=0; j < count; j++){
	XWuu_i[i][j] = 0;
	for (k = 0; k < count; k++) {
	  XWuu_i[i][j] += XW_i[i][k] * uu_i[k][j];
	}
      }
    }
  
    /* Rprintf("Xuu_i for unit %d is \n", (n+1)); */
    /* PdoubleMatrix(Xuu_i, *n_cov, count); */
    


    for(i=0; i < *n_cov; i++) {
      for(j=0; j < *n_cov; j++){
	XWuuWX_i[i][j] = 0;
	for (k = 0; k < count; k++) {
	  XWuuWX_i[i][j] += XWuu_i[i][k] * WX_i[k][j];
	}
      }
    }
  
    /* Rprintf("XWuuWX_i for unit %d is \n", (n+1)); */
    /* PdoubleMatrix(XWuuWX_i, *n_cov, *n_cov); */
    

    /* Adding results in OmegaMatrix  */
    for(i=0; i < *n_cov; i++) {
      for(j=0; j < *n_cov; j++){
	OmegaMatrix[i][j] +=  XWuuWX_i[i][j];
      }
    }



    FreeMatrix(X_i, count);
    FreeMatrix(tX_i, *n_cov);
    free(u_i);
    FreeMatrix(uu_i, count);
    FreeMatrix(XWuu_i, *n_cov);
    FreeMatrix(XWuuWX_i, *n_cov);
    free(u_i2);
    free(w_i);
    FreeMatrix(W_i, count);
    FreeMatrix(XW_i, *n_cov);
    FreeMatrix(WX_i, count);



  }
   
  /* storing results */

  itemp = 0;
  for (i = 0; i < *n_cov; i++) {
    for (j = 0; j < *n_cov; j++) {
      Omega_DiD_HC[itemp++] = OmegaMatrix[i][j];
    }
  }

  FreeMatrix(X, *len_data);
  FreeMatrix(OmegaMatrix, *n_cov);
  free(u);
  free(w);

}








void LamdaDID1(int* len_Xtrow, int* len_Xhrow,
	       int* Tunit_index, int* len_uniq_Tu_index,
	       int* Hunit_index, int* len_uniq_Hu_index,
	       double* Xtilde,
	       int* len_Xtcol, /* number of columns of Xtilde */
	       double* utilde,
	       double* Xhat,
	       int* len_Xhcol, /* number of columns of Xhat */
	       double* uhat,
	       double* W,
	       double* LamdaDID1) {


  /* temporary storages */
  int i, j, k, itemp;
  int n;
  int tcount;
  int hcount;
  double **Xt;
  double **Xh;
  double *ut;
  double *uh;
  double *w;
  double **LamdaMatrix;

  /* defining vectors and matricies */
  Xt = doubleMatrix(*len_Xtrow, *len_Xtcol);
  Xh = doubleMatrix(*len_Xhrow, *len_Xhcol);
  ut = doubleArray(*len_Xtrow);
  uh = doubleArray(*len_Xhrow);
  w = doubleArray(*len_Xtrow);
  LamdaMatrix = doubleMatrix(*len_Xhcol, *len_Xtcol);

  /* packing Xt */
  itemp = 0;
  for (j = 0; j < *len_Xtcol; j++) {
    for (i = 0; i < *len_Xtrow; i++) {
      Xt[i][j] = Xtilde[itemp++];
    }
  }



  /* packing Xh */
  itemp = 0;
  for (j = 0; j < *len_Xhcol; j++) {
    for (i = 0; i < *len_Xhrow; i++) {
      Xh[i][j] = Xhat[itemp++];
    }
  }
 
 
  /* packing LamdaMatrix */
  itemp = 0;
  for (j = 0; j < *len_Xtcol; j++) {
    for (i = 0; i < *len_Xhcol; i++) {
      LamdaMatrix[i][j] = 0;
    }
  }
 
  /* packing ut */
  for (i = 0; i < *len_Xtrow; i++) {
    ut[i] = utilde[i];
  }



  /* packing uh */
  for (i = 0; i < *len_Xhrow; i++) {
    uh[i] = uhat[i];
  }

  /* packing w */
  for (i = 0; i < *len_Xtrow; i++) {
    w[i] = W[i];
  }
   

  /* PdoubleArray(u, *len_data); */



   
  /* Rprintf("unit length is  %d\n", *len_uniq_u_index); */

  for (n = 0; n < *len_uniq_Hu_index; n++) {
    tcount = 0;
    hcount = 0;
    /* Rprintf("unit is  %d\n", (n+1)); */
    
    /* calculating t(X.tilda.i) %*% X.tilde.i */

    for (k = 0; k < *len_Xtrow; k++) {
      if (Tunit_index[k] == (n+1)) {
	tcount ++;
      }
    }
   
    for (k = 0; k < *len_Xhrow; k++) {
      if (Hunit_index[k] == (n+1)) {
	hcount ++;
      }
    }
   
    /* Rprintf("count for hat is %d count for tilde is %d\n", hcount, tcount); */


   
    double **Xt_i;
    double **Xh_i;
    double **tXh_i;
   
    double *ut_i;
    double *uh_i;
    double *w_i;
    double **W_i;
    double **WX_i;
    double **uu_i;
    double **uuWX_i;
    double **XuuWX_i;


      
    Xt_i = doubleMatrix(tcount, *len_Xtcol);
    Xh_i = doubleMatrix(hcount, *len_Xhcol);
    tXh_i = doubleMatrix(*len_Xhcol, hcount);
     
    ut_i = doubleArray(tcount);
    uh_i = doubleArray(hcount);
    w_i = doubleArray(tcount);
    W_i = doubleMatrix(tcount, tcount);

    WX_i = doubleMatrix(tcount, *len_Xtcol);
    uu_i = doubleMatrix(hcount, tcount);
    uuWX_i = doubleMatrix(hcount, *len_Xtcol);
    XuuWX_i = doubleMatrix(*len_Xhcol, *len_Xtcol);
      



    /* packing Xt_i */
    itemp = 0;
    for (k = 0; k < *len_Xtrow; k++) {
      if (Tunit_index[k] == (n+1)) {
	for (j = 0; j < *len_Xtcol; j++) {
	  Xt_i[itemp][j] = Xt[k][j];
	}
	itemp++;
      }
    }

    /* PdoubleMatrix(Xt_i, tcount, *len_Xtcol); */
      

    /* packing Xh_i */
    itemp = 0;
    for (k = 0; k < *len_Xhrow; k++) {
      if (Hunit_index[k] == (n+1)) {
	for (j = 0; j < *len_Xhcol; j++) {
	  Xh_i[itemp][j] = Xh[k][j];
	}
	itemp++;
      }
    }

    /* PdoubleMatrix(Xh_i, hcount, *len_Xhcol); */
      

    /* packing tXh_i */

    for (j = 0; j < *len_Xhcol; j++) {
      for (k = 0; k < hcount; k++) {
	tXh_i[j][k] = Xh_i[k][j];
      }
    }


    /* /\* PdoubleMatrix(tXt_i, *len_Xtcol, tcount); *\/ */
   

    /* packing ut_i */
    itemp = 0;
    for (k = 0; k < *len_Xtrow; k++) {
      if (Tunit_index[k] == (n+1)) {
	ut_i[itemp] = ut[k];
	itemp++;
      }
    }
      
    /* packing ut_i */
    itemp = 0;
    for (k = 0; k < *len_Xhrow; k++) {
      if (Hunit_index[k] == (n+1)) {
	uh_i[itemp] = uh[k];
	itemp++;
      }
    }
      
    /* PdoubleArray(uh_i, hcount); */
 

    /* packing w_i */
    itemp = 0;
    for (k = 0; k < *len_Xtrow; k++) {
      if (Tunit_index[k] == (n+1)) {
	w_i[itemp] = w[k];
	itemp++;
      }
    }
     
    /* PdoubleArray(w_i, tcount); */

    /* packing W_i */
    for(i=0; i < tcount; i++) {
      for(j=0; j < tcount; j++){
	if (i == j){
	  W_i[i][j] = w_i[i];
	} else {
	  W_i[i][j] = 0;
	}
      }
    }

    /*  /\* PdoubleMatrix(W_i, tcount, tcount); *\/ */

    /* packing WX_i */
    for(i=0; i < tcount; i++) {
      for(j=0; j < *len_Xtcol; j++){
	WX_i[i][j] = 0;
	for (k = 0; k < tcount; k++) {
	  WX_i[i][j] += W_i[i][k] * Xt_i[k][j];
	}
      }
    }
    /* PdoubleMatrix(WX_i, *len_Xtcol, tcount); */
   



    /*  /\* packing uu_i *\/ */

    for(i=0; i < hcount; i++) {
      for(j=0; j < tcount; j++){
	uu_i[i][j] = uh_i[i] * uh_i[j] ;
      }
    }
      
    /* PdoubleMatrix(uu_i, tcount, hcount); */
   

    /* packing uuWX_i */
    for(i=0; i < hcount; i++) {
      for(j=0; j < *len_Xtcol; j++){
	uuWX_i[i][j] = 0;
	for (k = 0; k < tcount; k++) {
	  uuWX_i[i][j] += uu_i[i][k] * WX_i[k][j];
	}
      }
    }

    /* PdoubleMatrix(uuWX_i, *len_Xtcol, hcount); */
   
 

    /* packing XuuWX_i */
    for(i=0; i < *len_Xhcol; i++) {
      for(j=0; j < *len_Xtcol; j++){
	XuuWX_i[i][j] = 0;
	for (k = 0; k < hcount; k++) {
	  XuuWX_i[i][j] += tXh_i[i][k] * uuWX_i[k][j];
	}
      }
    }

    /* PdoubleMatrix(XWuuX_i, *len_Xtcol, *len_Xhcol); */
   
 

    /* Adding results in LamdaMatrix  */
    for(i=0; i < *len_Xhcol; i++) {
      for(j=0; j < *len_Xtcol; j++){
	LamdaMatrix[i][j] +=  XuuWX_i[i][j];
      }
    }




    FreeMatrix(Xt_i, tcount);
    FreeMatrix(Xh_i, hcount);
    FreeMatrix(tXh_i, *len_Xhcol);

    free(ut_i);
    free(uh_i);
    free(w_i);
    FreeMatrix(W_i, tcount);

    FreeMatrix(WX_i, tcount);
    FreeMatrix(uu_i, tcount);
    FreeMatrix(uuWX_i, hcount);
    FreeMatrix(XuuWX_i, *len_Xhcol);



  }
   


  /* storing results */
   
  itemp = 0;
  for (i = 0; i < *len_Xhcol; i++) {
    for (j = 0; j < *len_Xtcol; j++) {
      LamdaDID1[itemp++] = LamdaMatrix[i][j];
    }
  }

  FreeMatrix(Xt, *len_Xtrow);
  FreeMatrix(Xh, *len_Xhrow);
  free(ut);
  free(uh);
  free(w);
  FreeMatrix(LamdaMatrix, *len_Xhcol);

   
}




























void LamdaDID2(int* len_Xtrow, int* len_Xhrow,
	       int* Tunit_index, int* len_uniq_Tu_index,
	       int* Hunit_index, int* len_uniq_Hu_index,
	       double* Xtilde,
	       int* len_Xtcol, /* number of columns of Xtilde */
	       double* utilde,
	       double* Xhat,
	       int* len_Xhcol, /* number of columns of Xhat */
	       double* uhat,
	       double* W,
	       double* LamdaDID2) {


  /* temporary storages */
  int i, j, k, itemp;
  int n;
  int tcount;
  int hcount;
  double **Xt;
  double **Xh;
  double *ut;
  double *uh;
  double *w;
  double **LamdaMatrix;

  /* defining vectors and matricies */
  Xt = doubleMatrix(*len_Xtrow, *len_Xtcol);
  Xh = doubleMatrix(*len_Xhrow, *len_Xhcol);
  ut = doubleArray(*len_Xtrow);
  uh = doubleArray(*len_Xhrow);
  w = doubleArray(*len_Xtrow);
  LamdaMatrix = doubleMatrix(*len_Xtcol, *len_Xhcol);

  /* packing Xtil */
  itemp = 0;
  for (j = 0; j < *len_Xtcol; j++) {
    for (i = 0; i < *len_Xtrow; i++) {
      Xt[i][j] = Xtilde[itemp++];
    }
  }



  /* packing Xhat */
  itemp = 0;
  for (j = 0; j < *len_Xhcol; j++) {
    for (i = 0; i < *len_Xhrow; i++) {
      Xh[i][j] = Xhat[itemp++];
    }
  }
 
 
  /* packing LamdaMatrix */
  itemp = 0;
  for (j = 0; j < *len_Xhcol; j++) {
    for (i = 0; i < *len_Xtcol; i++) {
      LamdaMatrix[i][j] = 0;
    }
  }
 
  /* packing ut */
  for (i = 0; i < *len_Xtrow; i++) {
    ut[i] = utilde[i];
  }



  /* packing uh */
  for (i = 0; i < *len_Xhrow; i++) {
    uh[i] = uhat[i];
  }

  /* packing w */
  for (i = 0; i < *len_Xtrow; i++) {
    w[i] = W[i];
  }
   

  /* PdoubleArray(u, *len_data); */



   
  /* Rprintf("unit length is  %d\n", *len_uniq_u_index); */

  for (n = 0; n < *len_uniq_Hu_index; n++) {
    tcount = 0;
    hcount = 0;
    /* Rprintf("unit is  %d\n", (n+1)); */
    
    /* calculating t(X.tilda.i) %*% X.tilde.i */

    for (k = 0; k < *len_Xtrow; k++) {
      if (Tunit_index[k] == (n+1)) {
	tcount ++;
      }
    }
   
    for (k = 0; k < *len_Xhrow; k++) {
      if (Hunit_index[k] == (n+1)) {
	hcount ++;
      }
    }
   
    /* Rprintf("count for hat is %d count for tilde is %d\n", hcount, tcount); */


   
    double **Xt_i;
    double **Xh_i;
    double **tXt_i;
    double *ut_i;
    double *uh_i;
    double *w_i;
    double **W_i;
    double **XW_i;
    double **uu_i;
    double **XWuu_i;
    double **XWuuX_i;


      
    Xt_i = doubleMatrix(tcount, *len_Xtcol);
    Xh_i = doubleMatrix(hcount, *len_Xhcol);
    tXt_i = doubleMatrix(*len_Xtcol, tcount);
     
    ut_i = doubleArray(tcount);
    uh_i = doubleArray(hcount);
    w_i = doubleArray(tcount);
    W_i = doubleMatrix(tcount, tcount);

    XW_i = doubleMatrix(*len_Xtcol, tcount);
    uu_i = doubleMatrix(tcount, hcount);
    XWuu_i = doubleMatrix(*len_Xtcol, hcount);
    XWuuX_i = doubleMatrix(*len_Xtcol, *len_Xhcol);
      



    /* packing Xt_i */
    itemp = 0;
    for (k = 0; k < *len_Xtrow; k++) {
      if (Tunit_index[k] == (n+1)) {
	for (j = 0; j < *len_Xtcol; j++) {
	  Xt_i[itemp][j] = Xt[k][j];
	}
	itemp++;
      }
    }

    /* PdoubleMatrix(Xt_i, tcount, *len_Xtcol); */
      

    /* packing Xh_i */
    itemp = 0;
    for (k = 0; k < *len_Xhrow; k++) {
      if (Hunit_index[k] == (n+1)) {
	for (j = 0; j < *len_Xhcol; j++) {
	  Xh_i[itemp][j] = Xh[k][j];
	}
	itemp++;
      }
    }

    /* PdoubleMatrix(Xh_i, hcount, *len_Xhcol); */
      

    /* packing tXt_i */

    for (j = 0; j < *len_Xtcol; j++) {
      for (k = 0; k < tcount; k++) {
	tXt_i[j][k] = Xt_i[k][j];
      }
    }


    /* PdoubleMatrix(tXt_i, *len_Xtcol, tcount); */
   

    /* packing ut_i */
    itemp = 0;
    for (k = 0; k < *len_Xtrow; k++) {
      if (Tunit_index[k] == (n+1)) {
	ut_i[itemp] = ut[k];
	itemp++;
      }
    }
      
    /* packing ut_i */
    itemp = 0;
    for (k = 0; k < *len_Xhrow; k++) {
      if (Hunit_index[k] == (n+1)) {
	uh_i[itemp] = uh[k];
	itemp++;
      }
    }
      
    /* PdoubleArray(uh_i, hcount); */
 

    /* packing w_i */
    itemp = 0;
    for (k = 0; k < *len_Xtrow; k++) {
      if (Tunit_index[k] == (n+1)) {
	w_i[itemp] = w[k];
	itemp++;
      }
    }
     
    /* PdoubleArray(w_i, tcount); */

    /* packing W_i */
    for(i=0; i < tcount; i++) {
      for(j=0; j < tcount; j++){
	if (i == j){
	  W_i[i][j] = w_i[i];
	} else {
	  W_i[i][j] = 0;
	}
      }
    }

    /* PdoubleMatrix(W_i, tcount, tcount); */

    /* packing XW_i */
    for(i=0; i < *len_Xtcol; i++) {
      for(j=0; j < tcount; j++){
	XW_i[i][j] = 0;
	for (k = 0; k < tcount; k++) {
	  XW_i[i][j] += tXt_i[i][k] * W_i[k][j];
	}
      }
    }
    /* PdoubleMatrix(XW_i, *len_Xtcol, tcount); */
   



    /* packing uu_i */

    for(i=0; i < tcount; i++) {
      for(j=0; j < hcount; j++){
	uu_i[i][j] = ut_i[i] * uh_i[j];
      }
    }
      
    /* PdoubleMatrix(uu_i, tcount, hcount); */
   

    /* packing XWuu_i */
    for(i=0; i < *len_Xtcol; i++) {
      for(j=0; j < hcount; j++){
	XWuu_i[i][j] = 0;
	for (k = 0; k < tcount; k++) {
	  XWuu_i[i][j] += XW_i[i][k] * uu_i[k][j];
	}
      }
    }

    /* PdoubleMatrix(XWuu_i, *len_Xtcol, hcount); */
   
 

    /* packing XWuuX_i */
    for(i=0; i < *len_Xtcol; i++) {
      for(j=0; j < *len_Xhcol; j++){
	XWuuX_i[i][j] = 0;
	for (k = 0; k < *len_Xhcol; k++) {
	  XWuuX_i[i][j] += XWuu_i[i][k] * Xh_i[k][j];
	}
      }
    }

    /* PdoubleMatrix(XWuuX_i, *len_Xtcol, *len_Xhcol); */
   
 

    /* Adding results in LamdaMatrix  */
    for(i=0; i < *len_Xtcol; i++) {
      for(j=0; j < *len_Xhcol; j++){
	LamdaMatrix[i][j] +=  XWuuX_i[i][j];
      }
    }




    FreeMatrix(Xt_i, tcount);
    FreeMatrix(Xh_i, hcount);
    FreeMatrix(tXt_i, *len_Xtcol);

    free(ut_i);
    free(uh_i);
    free(w_i);
    FreeMatrix(W_i, tcount);

    FreeMatrix(XW_i, *len_Xtcol);
    FreeMatrix(uu_i, tcount);
    FreeMatrix(XWuu_i, *len_Xtcol);
    FreeMatrix(XWuuX_i, *len_Xtcol);



  }
   


  /* storing results */
   
  itemp = 0;
  for (i = 0; i < *len_Xtcol; i++) {
    for (j = 0; j < *len_Xhcol; j++) {
      LamdaDID2[itemp++] = LamdaMatrix[i][j];
    }
  }

  FreeMatrix(Xt, *len_Xtrow);
  FreeMatrix(Xh, *len_Xhrow);
  free(ut);
  free(uh);
  free(w);
  FreeMatrix(LamdaMatrix, *len_Xtcol);

   
}






