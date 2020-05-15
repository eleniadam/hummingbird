#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#define idx2c(i ,j , ld ) ((( j )*( ld ))+( i ))

/* Functions accessible by R */

SEXP hummingbirdPostAdjustment(SEXP em, SEXP pos, SEXP minCpGs, SEXP minLength, SEXP maxGap);

SEXP hummingbirdEM(SEXP normM, SEXP normUM, SEXP abnormM, SEXP abnormUM, SEXP pos, SEXP binSize);


