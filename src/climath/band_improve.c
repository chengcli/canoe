#include <stdlib.h>
#include "linalg.h"

#undef IMIN
#define IMIN(i,j) ({ \
         const int _i = (int)(i); \
         const int _j = (int)(j); \
         _i < _j ? _i : _j; })

#undef IMAX
#define IMAX(i,j) ({ \
         const int _i = (int)(i); \
         const int _j = (int)(j); \
         _i > _j ? _i : _j; })


/*======================= band_improve() =====================================*/
/*
 * Based on Numerical Recipes in C, Secion 2.5, Iterative Improvement of
 * a Solution to Linear Equations.
 * This is for the band-diagnonal matrix case, and is analogous to lu_improve().
 * AORIG is the original matrix in compact form, whereas A and AL are the
 * matrices returned from band_decomp().
 * NOTE: The functionality of band_multiply() is echoed here because of the
 *       requirement of double precision.
 * Assumes zero-based indexing.
 */

#undef  AORIG
#define AORIG(i,j) aorig[(m1+m2+1)*(i)+(j)]

void band_improve(int     n,
                  int     m1,
                  int     m2,
                  double  *aorig,
                  double  *a,
                  double  *al,
                  int    *index,
                  double  *b,
                  double  *x)
{
  int
    k,j,i,tmploop;
  static int
    n_max = 0;
  double
   *r;
  double
    sdp;  /* NOTE: sdp must be double precision. */
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="band_improve";

  /*
   * Allocate memory.
   */
  if (n_max == 0) {
    n_max = n;
    r     = (double*)malloc(n_max*sizeof(double));
  }
  else if (n > n_max) {
    free(r);
    n_max = n;
    r     = (double*)malloc(n_max*sizeof(double));
  }

  /*
   * The band-diagonal indexing is as in band_multiply().
   */
  for (i = 0; i < n; i++) {
    k       = i-m1;
    tmploop = IMIN(m1+m2+1,n-k);
    sdp     = -b[i];
    for (j = IMAX(0,-k); j < tmploop; j++) {
      sdp += AORIG(i,j)*x[j+k];
    }
    r[i] = sdp;
  }

  band_back_sub(n,m1,m2,a,al,index,r);

  for (i = 0; i < n; i++) {
    x[i] -= r[i];
  }

  free(r);

  return;
}

/*======================= end of band_improve() ==============================*/
