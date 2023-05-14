#include <stdio.h>
#include <math.h>

/*======================= band_decomp() ======================================*/

/*
 * Decompose a banded matrix.
 * Adapted from Numerical Recipes in C, 2nd ed., p. 53.
 * The input matrix must be stored in compact form, as described on p. 52.
 * Assumes zero-based indexing.
 */

#undef  SWAP
#define SWAP(a,b) {tmp=(a);(a)=(b);(b)=tmp;}
#undef  TINY
#define TINY 1.e-20

#undef  A
#define A(i,j) a[(m1+m2+1)*(i)+(j)]
#undef  AL
#define AL(i,j) al[m1*(i)+(j)]

void band_decomp(int     n,
                 int     m1,
                 int     m2,
                 double  *a,
                 double  *al,
                 int    *index,
                 double  *d)
{
  int
    i,j,k,l,
    mm;
  double
    tmp;
  static int
    warned = 0;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="band_decomp";

  mm = m1+m2+1;
  l  = m1;
  for (i = 0; i < m1; i++) {
    for (j = m1-i; j < mm; j++) {
      A(i,j-l) = A(i,j);
    }
    l--;
    for (j = mm-l-1; j < mm; j++) {
      A(i,j) = 0.;
    }
  }
  *d = 1.;
  l  = m1;
  for (k = 0; k < n; k++) {
    tmp = A(k,0);
    i   = k;
    if (l < n) {
      l++;
    }
    for (j = k+1; j < l; j++) {
      if (fabs(A(j,0)) > fabs(tmp)) {
        tmp = A(j,0);
        i   = j;
      }
    }
    index[k] = i;
    if (tmp == 0.) {
      /*
       * Matrix is algorithmically singular.
       */
      A(k,0) = TINY;
      if (!warned) {
        fprintf(stderr,"**warning: %s, matrix is algorithmically singular (future warnings suppressed)\n",dbmsname);
        warned = 1;
      }
    }
    if (i != k) {
      *d = -(*d);
      for (j = 0; j < mm; j++) {
        SWAP(A(k,j),A(i,j))
      }
    }
    for (i = k+1; i < l; i++) {
      tmp          = A(i,0)/A(k,0);
      AL(k,i-k-1) = tmp;
      for (j = 1; j < mm; j++) {
        A(i,j-1) = A(i,j)-tmp*A(k,j);
      }
      A(i,mm-1) = 0.;
    }
  }

  return;
}

/*======================= end of band_decomp() ===============================*/
