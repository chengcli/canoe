template<int N>
int ludcmp(double a[N][N], int indx[N], double vv[N])
{
	int i, imax, j, k, d;
	double big, dum, sum, temp;

	d = 1;
	for (i = 0; i < N; i++) {
		big = 0.0;
		for (j = 0; j < N; j++)
      if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) {
			fprintf(stderr, "Singular matrix in routine ludcmp");
			exit(1);
		}
		vv[i] = 1.0 / big;
	}
	for (j = 0; j < N; j++) {
		for (i = 0; i < j; i++) {
			sum = a[i][j];
			for (k = 0; k < i; k++) sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
		}
		big = 0.0;
		for (i = j; i < N; i++) {
			sum = a[i][j];
			for (k = 0; k < j; k++)
				sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
			if ((dum = vv[i] * fabs(sum)) >= big) {
				big = dum;
				imax = i;
			}
		}
		if (j != imax) {
			for (k = 0; k < N; k++){
				dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = dum;
			}
			d = -d;
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		if (j != N - 1) {
			dum = (1.0 / a[j][j]);
			for (i = j + 1; i < N; i++) a[i][j] *= dum;
		}
	}

  return d;
}

template<int N>
void lubksb(double a[N][N], int indx[N], double b[N])
{
	int i, ii = 0, ip, j;
	double sum;

	for (i = 0; i < N; i++) {
		ip = indx[i];
		sum = b[ip];
		b[ip] = b[i];
		if (ii)
      for (j = ii-1; j < i; j++) sum -= a[i][j] * b[j];
		else if (sum) ii = i+1;
		b[i] = sum;
	}
	for (i = N-1; i >= 0; i--) {
		sum = b[i];
		for (j = i + 1; j < N; j++) sum -= a[i][j] * b[j];
		b[i] = sum / a[i][i];
	}
}
