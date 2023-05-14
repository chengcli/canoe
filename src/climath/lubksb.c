/*! Solves the set of n linear equations A X = B. Here
 * a[0..n-1][0..n-1] is input, not as the matrix A but rather as its LU decomposition, 
                 determined by the routine ludcmp.
 * indx[0..n-1] is input as the permutation vector returned by ludcmp.
 * b[0..n-1] is input as the right-hand side vector B, and returns with the solution vector X.
 * a, n, and indx are not modified by this routine and can be left in place for 
 * successive calls with different right-hand sides b. 
 * This routine takes into account the possibility that b will begin with many zero elements, 
 * so it is efficient for use in matrix inversion.
 * adapted from Numerical Recipes in C, 2nd Ed., p. 47.
 */
void lubksb(double **a, int n, int *indx, double *b)
{
	int i, ii = 0, ip, j;
	double sum;

	for (i = 0; i < n; i++) {
		ip = indx[i];
		sum = b[ip];
		b[ip] = b[i];
		if (ii)
      for (j = ii-1; j < i; j++) sum -= a[i][j] * b[j];
		else if (sum) ii = i+1;
		b[i] = sum;
	}
	for (i = n-1; i >= 0; i--) {
		sum = b[i];
		for (j = i + 1; j < n; j++) sum -= a[i][j] * b[j];
		b[i] = sum / a[i][i];
	}
}
