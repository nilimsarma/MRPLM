#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "define_t.h"

void XTWX(long int y_rows, long int y_cols, double *wts, double *xtwx);
void XTWXinv(long int y_rows, long int y_cols,double *xtwx);
void XTWY(long int y_rows, long int y_cols, double *wts,double *y, double *xtwy);
int Choleski_inverse(double *X, double *Xinv, double *work, int n, int upperonly);

