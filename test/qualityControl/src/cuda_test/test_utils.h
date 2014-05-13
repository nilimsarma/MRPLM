#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "string.h"
#include "define_t.h"

void gpu_XTWY(int_t y_rows, int_t y_cols, double *wts,double *y, double *xtwy);

void mat_transpose(double* mat, double* mat_trans, int_t y_cols, int_t y_rows);

void mat_subtract(double* A, double* B, double* C, int_t y_cols, int_t y_rows);

void mat_negate(double* mat, double* mat_neg, int_t y_cols, int_t y_rows);

void mat_inverse(double* mat, double* mat_inv, int_t y_rows);

void mat_mult(const double *A, const double *B, double *C, const int_t m, const int_t p, const int_t n);

void host_D_from_Dtemp(double* h_D, double* h_D_temp, int_t y_cols, int_t y_rows);

void print_mat(double* A, int_t y_cols, int_t y_rows);

extern "C" int Choleski_inverse(double *X, double *Xinv, double *work, int n, int upperonly);

