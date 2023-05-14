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

/*======================= band_multiply() ====================================*/
/*
 *  Compute matrix muliplication b = A.x, where A is in the compact-storage
 *  form of a band-diagonal matrix. 
 *  Based on Numerical Recipes in C, banmul(), p. 52.
 *  Assumes zero-based indexing.
 */

#undef  A
#define A(i,j) a[(m1+m2+1)*(i)+(j)]

void band_multiply(int     n,
                   int     m1,
                   int     m2,
                   double  *a,
                   double  *x,
                   double  *b)
{
  int
    i,j,k,tmploop;

  for (i = 0; i < n; i++) {
    k       = i-m1;
    tmploop = IMIN(m1+m2+1,n-k);
    b[i]    = 0.;
    for (j = IMAX(0,-k); j < tmploop; j++) {
      b[i] += A(i,j)*x[j+k];
    }
  }

  return;
}
                   
/*======================= end of band_multiply() =============================*/
