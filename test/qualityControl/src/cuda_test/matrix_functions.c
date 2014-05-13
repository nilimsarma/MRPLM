#include "matrix_functions.h"

static int use_lapack = 1;

/*********************************************************************
 **
 ** int Choleski_decompose(double *X, double *L, int n)
 **
 ** double *X - a square matrix 
 ** double *L - space already allocated to store Cholesky decomposition
 ** int n - matrix dimension
 ** int lapack - if 0 use LINPACK otherwise LAPACK
 **
 ** RETURNS integer code indicating success or failure 0 means no error, 
 **      non zero indicates failure ie not positive definite
 **
 ** this function choleski decomposes a positive definite symmetric matrix,
 ** on output L will contain the Choleski decomposition in the upper triangle
 **
 **
 *********************************************************************/


static int Choleski_decompose(double *X, double *L, int n, int lapack){
  int i,j,error_code;
  char upper = 'U';

  for (i=0; i < n; i++){
    for (j=0; j < n; j++){
      if (i > j)
	L[j*n+i] = 0.0;
      else {
	L[j*n+i] = X[j*n + i];
      }
    }
  }
#if 0
  if (!lapack){
    dpofa_(L,&n,&n,&error_code);
  } else {
    dpotrf_(&upper,&n,L,&n,&error_code);
  }
#endif    
	error_code = LAPACKE_dpotrf( LAPACK_COL_MAJOR, upper, n, L, n );


  return error_code;
}

/***********************************************************************
 **
 ** int Choleski_2_inverse(double *X, double *Xinv, int n)
 **
 ** double *X - matrix containing choleski decomposition in upper triangle
 ** double *Xinv - on output will contain the inverse
 ** int n - dimension of matrix
 ** int upperonly - if non zero return only the upper triangle of the inverse.
 ** int lapack - use LINPACK if 0 otherwise LAPACK
 **
 ** RETURNS integer code, indicating success 0  or error (non zero) 
 **
 ** this function uses the choleski decomposition of a 
 ** matrix to compute the inverse of a matrix.
 ** typically it would be used in conjunction with the choleski_decompose
 ** function above.
 **
 **
 **********************************************************************/

static int Choleski_2_inverse(double *X, double *Xinv, int n,int upperonly, int lapack){
  
  int i,j ,error_code=0,inverseonly;
  double d =0.0;
  char upper = 'U';
  
  for (i=0; i < n; i++){ 
    /* check for a zero or close to zero diagonal element */ 
    if(fabs(X[i*n+ i]) < 1e-06){
      error_code = 1;
      return error_code;
    }

    for (j=i; j < n; j++){
      Xinv[j*n + i] = X[j*n + i];
    }
  }

  inverseonly = 1;

#if 0
  if (!lapack){
    //dpodi_(Xinv,&n,&n,&d,&inverseonly);
  } else {
    //dpotri_(&upper,&n,Xinv,&n,&error_code);
    
  }
#endif
	error_code = LAPACKE_dpotri( LAPACK_COL_MAJOR, upper, n, Xinv, n );	
	
  if (!upperonly){
    for (i=0; i < n; i++){
      for (j=0; j <= i-1; j++){
		Xinv[j*n+i] = Xinv[i*n+j];
      }
    }
  }
  return error_code;

}

/***********************************************************************
 **
 ** int Choleski_inverse(double *X, double *Xinv, double *work, int n)
 **
 ** double *X - matrix containing choleski decomposition in upper triangle
 ** double *Xinv - on output will contain the inverse
 ** double *work - working space n*n dimension
 ** int n - dimension of matrix
 ** int upperonly - if non zero return only upper triangle of inverse.
 **
 ** RETURNS integer code, indicating success 0  or error (non zero) 
 **
 ** This function will compute the inverse of a positive definite symmetric
 ** matrix using choleski decomposition.
 **
 **********************************************************************/

int Choleski_inverse(double *X, double *Xinv, double *work, int n, int upperonly){

  int error_code;
  
  error_code = Choleski_decompose(X, work, n,use_lapack);
  if (!error_code){
    error_code = Choleski_2_inverse(work, Xinv, n,upperonly,use_lapack);
  }
  return error_code;

}

