/*======================= band_back_sub() ====================================*/

/*
 * Back substitute for a banded matrix using the output from band_decomp().
 * Adapted from Numerical Recipes in C, 2nd ed., p. 54.
 * Assumes zero-based indexing.
 */

#undef  SWAP
#define SWAP(a,b) {tmp=(a);(a)=(b);(b)=tmp;}

#undef  A
#define A(i,j) a[(m1+m2+1)*(i)+(j)]
#undef  AL
#define AL(i,j) al[m1*(i)+(j)]

void band_back_sub(int     n,
                   int     m1,
                   int     m2,
                   double  *a,
                   double  *al,
                   int    *index,
                   double  *b)
{
  int
    i,k,l,
    mm;
  double
    tmp;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="band_back_sub";

  mm = m1+m2+1;
  l  = m1;
  for (k = 0; k < n; k++) {
    i = index[k];
    if (i != k) {
      SWAP(b[k],b[i]);
    }
    if (l < n) {
      l++;
    }
    for (i = k+1; i < l; i++) {
      b[i] -= AL(k,i-k-1)*b[k];
    }
  }
  l = 1;
  for (i = n-1; i >= 0; i--) {
    tmp = b[i];
    for (k = 1; k < l; k++) {
      tmp -= A(i,k)*b[k+i];
    }
    b[i] = tmp/A(i,0);
    if (l < mm) {
      l++;
    }
  }

  return;
}

/*======================= end of band_back_sub() =============================*/
