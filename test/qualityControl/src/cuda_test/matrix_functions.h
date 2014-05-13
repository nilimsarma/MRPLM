#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lapacke.h"

static int Choleski_decompose(double *X, double *L, int n, int lapack);
static int Choleski_2_inverse(double *X, double *Xinv, int n,int upperonly, int lapack);
int Choleski_inverse(double *X, double *Xinv, double *work, int n, int upperonly);

