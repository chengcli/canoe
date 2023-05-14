#include <stdlib.h>
#include <stdio.h>
#include "linalg.h"

#define WITHOUT_PIVOTING 0
#define WITH_PIVOTING    1

/*======================= tridiag() =========================================*/

/*
 * Solve a tridiagnonal matrix system.
 * Adapted from Numerical Recipes in C, 2nd ed., pp. 51-54.
 * If pivot_type = WITH_PIVOTING, use band_decomp() and band_back_sub().
 * Assumes zero-based indexing.
 */

#undef  AA
#define AA(i,j) aa[(m1+m2+1)*(i)+(j)]
#undef  AAORIG
#define AAORIG(i,j) aaorig[(m1+m2+1)*(i)+(j)]
#undef  AAL
#define AAL(i,j) aal[m1*(i)+(j)]

void tridiag(int    n,
             double *a,
             double *b,
             double *c,
             double *r,
             double *u,
             int    pivot_type)
{
  int
    j,
    m1,m2,mm,
    *index;
  double
    bet,
    *gam,
    *aa,
    *aaorig,
    *aal,
     d;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="tridiag";

  /*
   * Check validity of n.
   */
  if (n <= 0) {
    fprintf(stderr,"**error:%s, n = %d\n",dbmsname,n);
    exit(1);
  }

  if (pivot_type == WITHOUT_PIVOTING) {
    /* Allocate memory. */
    gam = (double*)malloc(n*sizeof(double));

    if (b[0] == 0.0) {
      fprintf(stderr,"**error:%s,b[0] = 0\n"
                     "Rewrite equations as a set of order n-1, with u[1]\n"
                     "trivially eliminated\n",dbmsname);
      exit(1);
    }
    bet  = b[0];
    u[0] = r[0]/bet;
    for (j = 1; j < n; j++) {
      gam[j] = c[j-1]/bet;
      bet    = b[j]-a[j]*gam[j];
      if (bet == 0.) {
        /* 
         * Encountered a zero pivot. 
         * Try again using pivot_type = WITH_PIVOTING.
         */
        fprintf(stderr,"Warning: tridiag(): retrying with pivoting.\n");
        /* Free allocated memory. */
        free(gam);
        tridiag(n,a,b,c,r,u,WITH_PIVOTING);
        return;
      }
      u[j] =(r[j]-a[j]*u[j-1])/bet;
    }

    /* Backsubstitution: */
    for (j = n-2; j >= 0; j--) {
      u[j] -= gam[j+1]*u[j+1];
    }
    /* Free allocated memory. */
    free(gam);
    return;
  }
  else if (pivot_type == WITH_PIVOTING) {
    /*
     * Use band_decomp() and band_back_sub().
     */
    m1 = 1;
    m2 = 1;
    mm = m1+m2+1;
    /*
     * Allocate memory.
     */
    aa     = (double*)malloc(n*(m1+m2+1)*sizeof(double));
    aaorig = (double*)malloc(n*(m1+m2+1)*sizeof(double));
    aal    = (double*)malloc(n*m1*sizeof(double));
    index  = (int*)malloc(n*sizeof(int));

    /*
     * Load matrix AA and keep copy AAORIG.
     */
    for (j = 0; j < n; j++) {
      AA(j,m1+1) = AAORIG(j,m1+1) = c[j];
      AA(j,m1  ) = AAORIG(j,m1  ) = b[j];
      AA(j,m1-1) = AAORIG(j,m1-1) = a[j];
    }
    
    band_decomp(n,m1,m2,aa,aal,index,&d);

    /* 
     * Since tridiag() does not overwrite the input rhs vector, r,
     * with the answer, u, but band_back_sub() does, copy r into u
     * before calling band_back_sub().
     */
    for (j = 0; j < n; j++) {
      u[j] = r[j];
    }

    band_back_sub(n,m1,m2,aa,aal,index,u);

    /*
     *  Reduce roundoff errors with call to band_improve().
     */
    band_improve(n,m1,m2,aaorig,aa,aal,index,r,u);

    /*
     * Free allocated memory.
     */
    free(aa);
    free(aaorig);
    free(aal);
    free(index);

    return;
  }
  else {
    fprintf(stderr,"**error:%s, unrecognized pivot_type=%d\n",dbmsname,pivot_type);
    exit(1);
  }
}

/*======================= end of tridiag() ===================================*/
