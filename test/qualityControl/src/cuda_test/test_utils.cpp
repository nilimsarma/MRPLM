#include "test_utils.h"

void gpu_XTWY(int_t y_rows, int_t y_cols, double *wts,double *y, double *xtwy){

  int_t i,j;
    
  /* sweep columns (ie chip effects) */

  for (j=0; j < y_cols; j++){
    xtwy[j] = 0.0;
    for (i=0; i < y_rows; i++){
      xtwy[j] += wts[j*y_rows + i]* y[j*y_rows + i];
    }
  }
	
  /* sweep rows  (ie probe effects) */
 
  for (i=0; i < y_rows; i++){
    xtwy[i+y_cols] = 0.0;
    for (j=0; j < y_cols; j++){
      xtwy[i+y_cols] += wts[j*y_rows + i]* y[j*y_rows + i]; 
    }
  }

  for (i=0; i < y_rows-1; i++){
    xtwy[i+y_cols] = xtwy[i+y_cols] - xtwy[y_cols+y_rows-1];
  }
  
}

//Calculate transpose of matrix
void mat_transpose(double* mat, double* mat_trans, int_t y_rows, int_t y_cols){
	int_t i,j;
	for(j=0; j<y_cols; j++){
		for(i=0; i<y_rows; i++){
			mat_trans[i*y_cols+j] = mat[j*y_rows+i];
		}
	}
}

void mat_subtract(double* A, double* B, double* C, int_t y_rows, int_t y_cols){
	int_t i,j;
	for(j=0; j<y_cols; j++){
		for(i=0; i<y_rows; i++){
			C[j*y_rows +i] = A[j*y_rows +i] - B[j*y_rows +i];
		}
	}
}

void mat_negate(double* mat, double* mat_neg, int_t y_rows, int_t y_cols){
	int_t i,j;
	for(j=0; j<y_cols; j++){
		for(i=0; i<y_rows; i++){
			mat_neg[j*y_rows +i]= -mat[j*y_rows +i];
		}
	}
}

void mat_inverse(double* mat, double* mat_inv, int_t y_rows)
{
	double* work = (double*)calloc(y_rows*y_rows,sizeof(double));
	int err = Choleski_inverse(mat, mat_inv, work, (int)y_rows, 0);
	
	if(err) printf("\nCholeski error code = %d\n", err);
	
	free(work);
}

//C(m,n) = A(m,p) * B(p,n)
void mat_mult(const double *A, const double *B, double *C, const int_t m, const int_t p, const int_t n){

	int_t i,j,k;
	for (j = 0; j<n; j++){
		for (i = 0; i<m; i++){
			C[j*m+i] = 0.0;
			for (k = 0; k<p; k++){
				C[j*m+i]+=A[k*m+i]*B[j*p+k];
			}
		}
	}
}
	
//Calcuate on Host -> D from D_temp
void host_D_from_Dtemp(double* h_D, double* h_D_temp, int_t y_rows, int_t y_cols){
	
	int_t i,j;

	double val_temp = 0.0; //h_D_temp[i][row-1]
	memset(h_D, 0, (y_rows-1)*(y_rows-1)*sizeof(double));
	for(i=0; i<CUSTOM_VAL; i++){
		for(j=0; j<(y_rows-1);j++){
			h_D[j*(y_rows-1)+j] += h_D_temp[i*y_rows+j];
		}
		val_temp += h_D_temp[i*y_rows+y_rows-1];
	}

	for(i=0; i<(y_rows-1); i++){
		for(j=0; j<(y_rows-1); j++){
			h_D[i*(y_rows-1)+j] += val_temp;
		}
	}
}

//print matrix
void print_mat(double* A, int_t y_rows, int_t y_cols){
	printf("\n");
	printf("[");
	int_t i,j;
	for (j=0; j<y_cols; j++){
		for(i=0;i<y_rows;i++){
			printf("%e ",A[j*y_rows+i]);
		}
		printf(";");
	}
	printf("]");
	printf("\n");
}

