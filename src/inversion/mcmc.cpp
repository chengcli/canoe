// C/C++
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

#if FITSOUTPUT
extern "C" {
#include <fitsio.h>
}
#endif

// climath
#include <climath/core.h>

// utils
#include <utils/ndarrays.hpp>

// canoe
#include <configure.hpp>

// inversion
#include "mcmc.hpp"

void _choldc(double **a, int n, double *p) {
  int i, j, k;
  double sum;

  for (i = 0; i < n; ++i)
    for (j = i; j < n; ++j) {
      for (sum = a[i][j], k = i - 1; k >= 0; --k) {
        sum -= a[i][k] * a[j][k];
      }
      if (i == j) {
        assert(sum > 0.);
        p[i] = sqrt(sum);
      } else {
        a[j][i] = sum / p[i];
      }
    }
}

// Compute tau directly only if tau < TAUMAX.
// Otherwise compute tau using the pairwise sum series
#define TAUMAX 2

// Compute autocovariances up to lag s = WINMULT*TAU
#define WINMULT 5

// The autocovariance array is double C[MAXLAG+1]
// so that C[s] makes sense for s = MAXLAG.
#define MAXLAG TAUMAX *WINMULT

// Stop and print an error message if the array is shorter
// than MINFAC * MAXLAG.
#define MINFAC 5

int _acor(double *mean, double *sigma, double *tau, double *X, int L) {
  *mean = 0.;  // Compute the mean of X ...
  for (int i = 0; i < L; i++) *mean += X[i];
  *mean = *mean / L;
  for (int i = 0; i < L; i++) X[i] -= *mean;  //  ... and subtract it away.

  if (L < MINFAC * MAXLAG) {
    // printf("Acor error 1: The autocorrelation time is too long relative to
    // the variance.\n");
    return 1;
  }

  double C[MAXLAG + 1];
  for (int s = 0; s <= MAXLAG; s++)
    C[s] =
        0.;  // Here, s=0 is the variance, s = MAXLAG is the last one computed.

  int iMax = L - MAXLAG;  // Compute the autocovariance function . . .
  for (int i = 0; i < iMax; i++)
    for (int s = 0; s <= MAXLAG; s++)
      C[s] += X[i] * X[i + s];  // ...  first the inner products ...
  for (int s = 0; s <= MAXLAG; s++)
    C[s] = C[s] / iMax;  // ...  then the normalization.

  double D =
      C[0];  // The "diffusion coefficient" is the sum of the autocovariances
  for (int s = 1; s <= MAXLAG; s++)
    D += 2 *
         C[s];  // The rest of the C[s] are double counted since C[-s] = C[s].
  *sigma = sqrt(
      D / L);  // The standard error bar formula, if D were the complete sum.
  *tau = D / C[0];  // A provisional estimate, since D is only part of the
                    // complete sum.

  if (*tau * WINMULT < MAXLAG) {
    return 0;  // Stop if the D sum includes the given multiple of tau.
               // This is the self consistent window approach.

  } else {  // If the provisional tau is so large that we don't think tau
            // is accurate, apply the acor procedure to the pairwase sums
            // of X.
    int Lh = L / 2;  // The pairwise sequence is half the length (if L is even)
    double newMean;  // The mean of the new sequence, to throw away.
    int j1 = 0;
    int j2 = 1;
    for (int i = 0; i < Lh; i++) {
      X[i] = X[j1] + X[j2];
      j1 += 2;
      j2 += 2;
    }
    _acor(&newMean, sigma, tau, X, Lh);
    D = .25 * (*sigma) * (*sigma) *
        L;  // Reconstruct the fine time series numbers from the coarse series
            // numbers.
    *tau = D / C[0];       // As before, but with a corrected D.
    *sigma = sqrt(D / L);  // As before, again.
  }

  return 0;
}

#undef TAUMAX
#undef WINMULT
#undef MAXLAG
#undef MINFAC

void mcmc_alloc(mcmc_recs *recs, int nstep, int nwalker, int ndim, int nvalue) {
  recs->nstep = nstep;
  recs->nwalker = nwalker;
  recs->ndim = ndim;
  recs->nvalue = nvalue;

  NewCArray(recs->par, nstep, nwalker, ndim);
  NewCArray(recs->val, nstep, nwalker, nvalue);
  NewCArray(recs->lnp, nstep, nwalker);
  NewCArray(recs->newstate, nstep, nwalker);

  memset(recs->par[0][0], 0, nstep * nwalker * ndim * sizeof(double));
  memset(recs->val[0][0], 0, nstep * nwalker * nvalue * sizeof(double));
  memset(recs->lnp[0], 0, nstep * nwalker * sizeof(double));
  memset(recs->newstate[0], 0, nstep * nwalker * sizeof(int));

  recs->opt_par = new double[ndim];
  recs->opt_val = new double[nvalue];

  recs->cur = 0;
  recs->accept = 0;
  recs->reset = 0;
  recs->opt_lnp = -1.E10;
}

void mcmc_free(mcmc_recs *recs) {
  FreeCArray(recs->par);
  FreeCArray(recs->val);
  FreeCArray(recs->lnp);
  FreeCArray(recs->newstate);
  delete[] recs->opt_par;
  delete[] recs->opt_val;
}

void mcmc_statistics(double *mean, double *sigma, double *tau,
                     mcmc_recs *recs) {
  double *x = new double[recs->cur];

  for (size_t p = 0; p < recs->ndim; ++p) {
    for (size_t t = 0; t < recs->cur; ++t) {
      x[t] = 0.;
      for (size_t k = 0; k < recs->nwalker; ++k) x[t] += recs->par[t][k][p];
      x[t] /= recs->nwalker;
    }

    int err = _acor(mean + p, sigma + p, tau + p, x, recs->cur);

    if (err != 0) tau[p] = -1.;
  }

  delete[] x;
}

void mcmc_save_fits(char const *fname, mcmc_opts *opts, mcmc_recs *recs,
                    int include_last) {
  double ***buf;
  int ***ibuf;

  double ***par, ***val, **lnp;
  int **newstate;

  // indices
  int rank = 0, size = 1;
  int nwalker = recs->nwalker;
  int ndim = recs->ndim;
  int nvalue = recs->nvalue;
  int cur = include_last ? recs->cur : recs->cur - 1;

  int accept = recs->accept;
  int reset = recs->reset;
  double opt_lnp = recs->opt_lnp;

  // do not make outputs if cur <= 0
  if (cur <= 0) return;

#ifdef MPI_PARALLEL
  MPI_Comm_rank(opts->mpi_comm, &rank);
  MPI_Comm_size(opts->mpi_comm, &size);

  NewCArray(par, cur, size * nwalker, ndim);
  NewCArray(val, cur, size * nwalker, nvalue);
  NewCArray(lnp, cur, size * nwalker);
  NewCArray(newstate, cur, size * nwalker);

  // 1. gather walker positions
  NewCArray(buf, size, cur * nwalker, ndim);
  MPI_Gather(**recs->par, cur * nwalker * ndim, MPI_DOUBLE, **buf,
             cur * nwalker * ndim, MPI_DOUBLE, 0, opts->mpi_comm);
  if (rank == 0) {  // reorder dimension
    for (int t = 0; t < cur; ++t)
      for (int r = 0; r < size; ++r)
        for (int k = 0; k < nwalker; ++k)
          for (int d = 0; d < ndim; ++d)
            par[t][r * nwalker + k][d] = buf[r][t * nwalker + k][d];
  }
  FreeCArray(buf);

  // 2. gather function values
  NewCArray(buf, size, cur * nwalker, nvalue);
  MPI_Gather(**recs->val, cur * nwalker * nvalue, MPI_DOUBLE, **buf,
             cur * nwalker * nvalue, MPI_DOUBLE, 0, opts->mpi_comm);
  if (rank == 0) {  // reorder dimension
    for (int t = 0; t < cur; ++t)
      for (int r = 0; r < size; ++r)
        for (int k = 0; k < nwalker; ++k)
          for (int d = 0; d < nvalue; ++d)
            val[t][r * nwalker + k][d] = buf[r][t * nwalker + k][d];
  }
  FreeCArray(buf);

  // 3. gather log probability
  NewCArray(buf, size, cur, nwalker);
  MPI_Gather(*recs->lnp, cur * nwalker, MPI_DOUBLE, **buf, cur * nwalker,
             MPI_DOUBLE, 0, opts->mpi_comm);
  if (rank == 0) {  // reorder dimension
    for (int t = 0; t < cur; ++t)
      for (int r = 0; r < size; ++r)
        for (int k = 0; k < nwalker; ++k)
          lnp[t][r * nwalker + k] = buf[r][t][k];
  }
  FreeCArray(buf);

  // 4. reduce optimal log probability
  MPI_Reduce(&recs->opt_lnp, &opt_lnp, 1, MPI_DOUBLE, MPI_MAX, 0,
             opts->mpi_comm);

  // 5. reduce accepted states
  MPI_Reduce(&recs->accept, &accept, 1, MPI_INT, MPI_SUM, 0, opts->mpi_comm);

  // 6. reduce resetted states
  MPI_Reduce(&recs->reset, &reset, 1, MPI_INT, MPI_SUM, 0, opts->mpi_comm);

  // 7. gather new state indicator
  NewCArray(ibuf, size, cur, nwalker);
  MPI_Gather(*recs->newstate, cur * nwalker, MPI_INT, **ibuf, cur * nwalker,
             MPI_INT, 0, opts->mpi_comm);
  if (rank == 0) {  // reorder dimension
    for (int t = 0; t < cur; ++t)
      for (int r = 0; r < size; ++r)
        for (int k = 0; k < nwalker; ++k)
          newstate[t][r * nwalker + k] = ibuf[r][t][k];
  }
  FreeCArray(ibuf);
#endif

  // 6. write FITS file
#if FITSOUTPUT
  if (rank == 0) {
    fitsfile *fp;
    int status = 0;
    int64_t naxis = 3;
    int64_t naxes[3] = {ndim, size * nwalker, cur};
    int64_t fpixel = 1;
    int64_t nelements = naxes[0] * naxes[1] * naxes[2];

    fits_create_file(&fp, fname, &status);
    if (status) fits_report_error(stderr, status);

    // 6.1 create parameter fields
    fits_create_img(fp, DOUBLE_IMG, naxis, naxes, &status);
    if (status) fits_report_error(stderr, status);

      // 6.2 write parameter fields
#ifdef MPI_PARALLEL
    fits_write_img(fp, TDOUBLE, fpixel, nelements, **par, &status);
#else
    fits_write_img(fp, TDOUBLE, fpixel, nelements, **recs->par, &status);
#endif
    if (status) fits_report_error(stderr, status);

    // 6.3 write headers
    char key[80] = "C.Li";
    fits_write_key(fp, TSTRING, "CREATOR", key, "file created by Cheng Li",
                   &status);

    strcpy(key, "par");
    fits_write_key(fp, TSTRING, "VAR", key, "retrieved parameters", &status);

    strcpy(key, "MCMC");
    fits_write_key(fp, TSTRING, "METHOD", key, "retrieval method", &status);

    fits_write_key(fp, TDOUBLE, "A", &opts->a, "stretch move parameter",
                   &status);
    fits_write_key(fp, TINT, "P", &opts->p, "walk move parameter", &status);

    fits_write_key(fp, TDOUBLE, "LOGP", &opt_lnp, "optimal log probability",
                   &status);
    fits_write_key(fp, TINT, "RESET", &reset, "number of resetted states",
                   &status);
    fits_write_key(fp, TINT, "ACCEPT", &accept, "number of accepted states",
                   &status);

    // 6.4 create function value fields, 2nd HDU
    naxes[0] = nvalue;
    nelements = naxes[0] * naxes[1] * naxes[2];
    fits_create_img(fp, DOUBLE_IMG, naxis, naxes, &status);
    if (status) fits_report_error(stderr, status);

      // write function values
#ifdef MPI_PARALLEL
    fits_write_img(fp, TDOUBLE, fpixel, nelements, **val, &status);
#else
    fits_write_img(fp, TDOUBLE, fpixel, nelements, **recs->val, &status);
#endif
    if (status) fits_report_error(stderr, status);

    strcpy(key, "val");
    fits_write_key(fp, TSTRING, "VAR", key, "values of forward model", &status);

    // 6.5 create log probability fields, 3rd HDU
    naxis = 2;
    nelements = naxes[1] * naxes[2];
    fits_create_img(fp, DOUBLE_IMG, naxis, naxes + 1, &status);
    if (status) fits_report_error(stderr, status);

    strcpy(key, "lnp");
    fits_write_key(fp, TSTRING, "VAR", key, "log probability", &status);

    // write log probability
#ifdef MPI_PARALLEL
    fits_write_img(fp, TDOUBLE, fpixel, nelements, *lnp, &status);
#else
    fits_write_img(fp, TDOUBLE, fpixel, nelements, *recs->lnp, &status);
#endif
    if (status) fits_report_error(stderr, status);

    // 6.6 create newstate masks, 4th HDU
    fits_create_img(fp, SHORT_IMG, naxis, naxes + 1, &status);
    if (status) fits_report_error(stderr, status);

    strcpy(key, "inew");
    fits_write_key(fp, TSTRING, "VAR", key, "new state indicator", &status);

    // write new state indicator
#ifdef MPI_PARALLEL
    fits_write_img(fp, TINT, fpixel, nelements, *newstate, &status);
#else
    fits_write_img(fp, TINT, fpixel, nelements, *recs->newstate, &status);
#endif
    if (status) fits_report_error(stderr, status);

    fits_close_file(fp, &status);
  }
#endif  // FITSOUTPUT

#ifdef MPI_PARALLEL
  FreeCArray(par);
  FreeCArray(val);
  FreeCArray(lnp);
  FreeCArray(newstate);
#endif
}

void mcmc_load_fits(char const *fname, mcmc_opts *opts, mcmc_recs *recs,
                    int alloc) {
  int status = 0, hdutype;
  int nstep, nwalker, ndim, nvalue;
  int64_t dims[3] = {1, 1, 1};
  int64_t fpixel[3] = {1, 1, 1};

#if FITSOUTPUT
  fitsfile *fp;

  // open FITS file
  fits_open_file(&fp, fname, READONLY, &status);
  if (status) {
    fits_report_error(stderr, status);
    exit(1);
  }

  // move to the first hdu, containing the parameters
  fits_movabs_hdu(fp, 1, &hdutype, &status);

  // read img dimensions
  fits_get_img_size(fp, 3, dims, &status);
  nstep = dims[2];
  nwalker = dims[1];
  ndim = dims[0];

  // move to the second hdu, containing the values of the forward model
  fits_movabs_hdu(fp, 2, &hdutype, &status);
  fits_get_img_size(fp, 3, dims, &status);
  nvalue = dims[0];

  // allocate memory
  if (alloc) {
    mcmc_alloc(recs, nstep, nwalker, ndim, nvalue);
  } else {
    assert(recs->nstep >= nstep);
    assert(recs->nwalker == nwalker);
    assert(recs->ndim == ndim);
    assert(recs->nvalue == nvalue);
  }

  // move to first hdu, read positions
  fits_movabs_hdu(fp, 1, &hdutype, &status);
  fits_read_pix(fp, TDOUBLE, fpixel, nstep * nwalker * ndim, NULL, **recs->par,
                NULL, &status);

  // read keywords
  fits_read_key(fp, TDOUBLE, "LOGP", &recs->opt_lnp, NULL, &status);
  fits_read_key(fp, TDOUBLE, "A", &opts->a, NULL, &status);
  fits_read_key(fp, TINT, "P", &opts->p, NULL, &status);
  fits_read_key(fp, TINT, "RESET", &recs->reset, NULL, &status);
  fits_read_key(fp, TINT, "ACCEPT", &recs->accept, NULL, &status);

  // move to second hdu, read values
  fits_movabs_hdu(fp, 2, &hdutype, &status);
  fits_read_pix(fp, TDOUBLE, fpixel, nstep * nwalker * nvalue, NULL,
                **recs->val, NULL, &status);

  // move to third hdu, read lnp
  fits_movabs_hdu(fp, 3, &hdutype, &status);
  fits_read_pix(fp, TDOUBLE, fpixel, nstep * nwalker, NULL, *recs->lnp, NULL,
                &status);

  // move to forth hdu, read new state indicator
  fits_movabs_hdu(fp, 4, &hdutype, &status);
  // fits_read_pix(fp, TSHORT, fpixel, nstep*nwalker, NULL,
  fits_read_pix(fp, TINT, fpixel, nstep * nwalker, NULL, *recs->newstate, NULL,
                &status);

  if (status) {
    fits_report_error(stderr, status);
    exit(1);
  }

  recs->cur = nstep;

  fits_close_file(fp, &status);
#endif
}

void mcmc_append_recs(mcmc_recs *dst, mcmc_recs *src) {
  assert(dst->ndim == src->ndim);
  assert(dst->nvalue == src->nvalue);
  assert(dst->nwalker == src->nwalker);
  assert(dst->nstep >= dst->cur + src->cur);

  int size_par = src->cur * src->nwalker * src->ndim;
  int size_val = src->cur * src->nwalker * src->nvalue;
  int size_lnp = src->cur * src->nwalker;
  int cur = dst->cur;

  memcpy(*dst->par[cur], **src->par, size_par * sizeof(double));
  memcpy(*dst->val[cur], **src->val, size_val * sizeof(double));
  memcpy(dst->lnp[cur], *src->lnp, size_lnp * sizeof(double));
  memcpy(dst->newstate[cur], *src->newstate, size_lnp * sizeof(int));

  if (dst->opt_lnp < src->opt_lnp) {
    memcpy(dst->opt_par, src->opt_par, src->ndim * sizeof(double));
    memcpy(dst->opt_val, src->opt_val, src->nvalue * sizeof(double));
    dst->opt_lnp = src->opt_lnp;
  }

  dst->cur += src->cur;
  dst->accept += src->accept;
  dst->reset += src->reset;
}

void mcmc_report(mcmc_opts *opts, mcmc_recs *recs, char const *mode) {
  int rank = 0, size = 1;
  double **par, *lnp, **opt_par, **opt_val;
  int accept;
  double *mean, *sigma, *tau, **cov;

  struct {
    double value;
    int rank;
  } opt_lnp;

  // indices
  int nwalker = recs->nwalker;
  int ndim = recs->ndim;
  int nvalue = recs->nvalue;
  int cur = recs->cur;

  // local statistics
  mean = new double[ndim];
  sigma = new double[ndim];
  tau = new double[ndim];
  NewCArray(cov, ndim, ndim);
  mcmc_statistics(mean, sigma, tau, recs);

  // amend the calculation of standard deviation
  for (int d = 0; d < ndim; ++d) {
    sigma[d] = 0.;
    for (int t = 0; t < cur; ++t)
      for (int k = 0; k < nwalker; ++k)
        sigma[d] += sqr(recs->par[t][k][d] - mean[d]);
  }

  // calculate correlation coefficients
  for (int i = 0; i < ndim; ++i)
    for (int j = 0; j < ndim; ++j) {
      cov[i][j] = 0;
      for (int t = 0; t < cur; ++t)
        for (int k = 0; k < nwalker; ++k)
          cov[i][j] +=
              (recs->par[t][k][i] - mean[i]) * (recs->par[t][k][j] - mean[j]);
    }

#ifdef MPI_PARALLEL
  MPI_Comm_rank(opts->mpi_comm, &rank);
  MPI_Comm_size(opts->mpi_comm, &size);

  lnp = new double[nwalker * size];
  NewCArray(par, nwalker * size, ndim);
  NewCArray(opt_par, size, ndim);
  NewCArray(opt_val, size, nvalue);

  // 1. gather walker positions
  MPI_Gather(*recs->par[cur - 1], nwalker * ndim, MPI_DOUBLE, *par,
             nwalker * ndim, MPI_DOUBLE, 0, opts->mpi_comm);

  // 2. gather log probability
  MPI_Gather(recs->lnp[cur - 1], nwalker, MPI_DOUBLE, lnp, nwalker, MPI_DOUBLE,
             0, opts->mpi_comm);

  // 3. gather optimal positions
  MPI_Gather(recs->opt_par, ndim, MPI_DOUBLE, *opt_par, ndim, MPI_DOUBLE, 0,
             opts->mpi_comm);

  // 4. gather optimal forward result
  MPI_Gather(recs->opt_val, nvalue, MPI_DOUBLE, *opt_val, nvalue, MPI_DOUBLE, 0,
             opts->mpi_comm);

  // 5. reduce optimal log probability with rank
  opt_lnp.value = recs->opt_lnp;
  opt_lnp.rank = rank;
  if (rank == 0)
    MPI_Reduce(MPI_IN_PLACE, &opt_lnp, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0,
               opts->mpi_comm);
  else
    MPI_Reduce(&opt_lnp, &opt_lnp, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0,
               opts->mpi_comm);

  // 6. reduce total number of accepted states
  MPI_Reduce(&recs->accept, &accept, 1, MPI_INT, MPI_SUM, 0, opts->mpi_comm);

  // 7. reduce mean value
  MPI_Allreduce(MPI_IN_PLACE, mean, ndim, MPI_DOUBLE, MPI_SUM, opts->mpi_comm);
  for (int d = 0; d < ndim; ++d) mean[d] /= size;

  // 8. reduce standard deviation
  if (rank == 0)
    MPI_Reduce(MPI_IN_PLACE, sigma, ndim, MPI_DOUBLE, MPI_SUM, 0,
               opts->mpi_comm);
  else
    MPI_Reduce(sigma, sigma, ndim, MPI_DOUBLE, MPI_SUM, 0, opts->mpi_comm);

  // 9. reduce correlation coefficients
  if (rank == 0)
    MPI_Reduce(MPI_IN_PLACE, *cov, ndim * ndim, MPI_DOUBLE, MPI_SUM, 0,
               opts->mpi_comm);
  else
    MPI_Reduce(*cov, *cov, ndim * ndim, MPI_DOUBLE, MPI_SUM, 0, opts->mpi_comm);

  // 10. reduce autocorrelation time
  if (rank == 0)
    MPI_Reduce(MPI_IN_PLACE, tau, ndim, MPI_DOUBLE, MPI_SUM, 0, opts->mpi_comm);
  else
    MPI_Reduce(tau, tau, ndim, MPI_DOUBLE, MPI_SUM, 0, opts->mpi_comm);
  for (int d = 0; d < ndim; ++d) tau[d] /= size;
#else
  par = recs->par[cur - 1];
  lnp = recs->lnp[cur - 1];
  opt_par = new double *[1];
  opt_par[0] = recs->opt_par;
  opt_val = new double *[1];
  opt_val[0] = recs->opt_val;
  opt_lnp.value = recs->opt_lnp;
  opt_lnp.rank = rank;
  accept = recs->accept;
#endif

  // final step for standard deviation
  for (int d = 0; d < ndim; ++d)
    sigma[d] = sqrt(sigma[d] / (nwalker * size * cur));

  // final step for correlation coefficients
  for (int i = 0; i < ndim; ++i)
    for (int j = 0; j < ndim; ++j) cov[i][j] /= size * nwalker * cur;

  // report
  if (rank == 0) {
    FILE *fout = fopen(opts->logfile, mode);

    fprintf(fout, "#### MCMC Iteration: No. %d ####\n", recs->reset + cur);
    // 1. walkers positions
    fprintf(fout, "1. Current positions:\n");
    fprintf(fout, "%3s", " ");
    for (int d = 0; d < ndim; ++d) fprintf(fout, "%9s%-2d", "D", d + 1);
    fprintf(fout, "\n");
    for (int k = 0; k < nwalker * size; k++) {
      fprintf(fout, "%3d", static_cast<int>(k) + 1);
      for (int d = 0; d < ndim; ++d) fprintf(fout, "%11.4g", par[k][d]);
      fprintf(fout, "  |%9.4f\n", lnp[k]);
      fflush(stdout);
    }

    // 2. optimal position
    fprintf(fout, "2. Optimal parameters:\n");
    fprintf(fout, "%3s", " ");
    for (int d = 0; d < ndim; ++d)
      fprintf(fout, "%11.4g", opt_par[opt_lnp.rank][d]);
    fprintf(fout, "\n");

    // 3. mean value
    fprintf(fout, "4. Mean value:\n");
    fprintf(fout, "%3s", " ");
    for (int d = 0; d < ndim; ++d) fprintf(fout, "%11.4g", mean[d]);
    fprintf(fout, "\n");

    // 4. standard deviation
    fprintf(fout, "4. Standard deviation:\n");
    fprintf(fout, "%3s", " ");
    for (int d = 0; d < ndim; ++d) fprintf(fout, "%11.4g", sigma[d]);
    fprintf(fout, "\n");

    // 5. correlation coefficient
    fprintf(fout, "5. Correlation coefficient:\n");
    for (int i = 0; i < ndim; ++i) {
      fprintf(fout, "%1s%-2d", "D", i + 1);
      for (int j = 0; j < ndim; ++j)
        fprintf(fout, "%11.4g", cov[i][j] / (sigma[i] * sigma[j]));
      fprintf(fout, "\n");
    }

    // 6. autocorrelation time
    fprintf(fout, "6. Autocorrelation time:\n");
    fprintf(fout, "%3s", " ");
    for (int d = 0; d < ndim; ++d) fprintf(fout, "%11.4g", tau[d]);
    fprintf(fout, "\n");

    // 7. optimal model results
    fprintf(fout, "7. Optimal model results:\n");
    for (int d = 0; d < nvalue; ++d) {
      fprintf(fout, "%11.4f", opt_val[opt_lnp.rank][d]);
      if ((d + 1) % 6 == 0) fprintf(fout, "\n");
    }
    if (nvalue % 6 != 0) fprintf(fout, "\n");

    // 8. optimal probability
    fprintf(fout, "8. Optimal log probability:\n");
    fprintf(fout, "%11.4g\n", opt_lnp.value);

    // 9. acceptance fraction
    fprintf(fout, "9. Acceptance fraction:\n");
    fprintf(fout, "%11.4g",
            1. * accept / ((recs->reset + cur) * nwalker * size));
    fprintf(fout, "\n##############################\n");

    fclose(fout);
  }

#ifdef MPI_PARALLEL
  delete[] lnp;
  FreeCArray(par);
  FreeCArray(opt_par);
  FreeCArray(opt_val);
#else
  delete[] opt_par;
  delete[] opt_val;
#endif

  delete[] mean;
  delete[] sigma;
  delete[] tau;
  FreeCArray(cov);
}

void mcmc_init(ObjectiveFunction_t lnprob, double **par, mcmc_opts *opts,
               mcmc_recs *recs, void *obj) {
  int nwalker = recs->nwalker;
  int ndim = recs->ndim;
  int nval = recs->nvalue;

  // make sure that the start points are valid
  for (int k = 0; k < nwalker; ++k) {
    recs->lnp[0][k] = lnprob(par[k], recs->val[0][k], ndim, nval, obj);

    int niter = 0;
    while (std::isnan(recs->lnp[0][k]) &&
           (niter++ < 10)) {  // if point (k) is invalid
      mcmc_stretch_move(par[k], par, k, nwalker, ndim, opts);
      recs->lnp[0][k] = lnprob(par[k], recs->val[0][k], ndim, nval, obj);
    }

    if (niter >= 10) {
      std::cerr << "Starting point iteration > 10 times" << std::endl;
      std::exit(1);
    }

    // transfer input parameters to records
    for (int d = 0; d < ndim; ++d) recs->par[0][k][d] = par[k][d];

    if (recs->lnp[0][k] > recs->opt_lnp) {
      recs->opt_lnp = recs->lnp[0][k];
      for (int d = 0; d < ndim; ++d) recs->opt_par[d] = recs->par[0][k][d];
      for (int d = 0; d < nval; ++d) recs->opt_val[d] = recs->val[0][k][d];
    }
    recs->newstate[0][k] = 1;
  }

  recs->cur++;
  recs->accept += nwalker;

  // make initial output
  mcmc_report(opts, recs, "w");
}

void mcmc_advance(ObjectiveFunction_t lnprob, mcmc_opts *opts, mcmc_recs *recs,
                  void *obj) {
  int cur = recs->cur;
  int ndim = recs->ndim;
  int nwalker = recs->nwalker;
  int nval = recs->nvalue;

  double *par = new double[ndim];

  unsigned int seed = time(NULL);
  for (int k = 0; k < nwalker; ++k) {
    double zz =
        mcmc_stretch_move(par, recs->par[cur - 1], k, nwalker, ndim, opts);
    // mcmc_walk_move(par_all + k*np, par, np, k, nwalker, zz + k, opts);

    recs->lnp[cur][k] = lnprob(par, recs->val[cur][k], ndim, nval, obj);

    double lnp0 = recs->lnp[cur - 1][k], lnp1 = recs->lnp[cur][k],
           pdiff = (ndim - 1.) * log(zz) + lnp1 - lnp0;

    if (pdiff > log(1. * rand_r(&seed) / RAND_MAX)) {  // accept this position
      for (int d = 0; d < ndim; ++d) recs->par[cur][k][d] = par[d];

      if (lnp1 > recs->opt_lnp) {  // new best position
        recs->opt_lnp = lnp1;
        for (int d = 0; d < ndim; ++d) recs->opt_par[d] = recs->par[cur][k][d];
        for (int d = 0; d < nval; ++d) recs->opt_val[d] = recs->val[cur][k][d];
      }

      recs->accept++;
      recs->newstate[cur][k] = 1;
    } else {  // do not accept this position
      for (int d = 0; d < ndim; ++d)
        recs->par[cur][k][d] = recs->par[cur - 1][k][d];
      for (int d = 0; d < nval; ++d)
        recs->val[cur][k][d] = recs->val[cur - 1][k][d];
      recs->lnp[cur][k] = recs->lnp[cur - 1][k];
      recs->newstate[cur][k] = 0;
    }
  }

  recs->cur++;
  if (recs->cur % opts->print == 0) mcmc_report(opts, recs, "a");
  delete[] par;
}

double mcmc_stretch_move(double *newp, double **oldp, int iwalker, int nwalker,
                         int ndim, mcmc_opts *opts) {
  int irank = 0, size = 1;
  double ***par;

#ifdef MPI_PARALLEL
  MPI_Comm_rank(opts->mpi_comm, &irank);
  MPI_Comm_size(opts->mpi_comm, &size);
#endif

  NewCArray(par, size, nwalker, ndim);

#ifdef MPI_PARALLEL
  MPI_Allgather(*oldp, nwalker * ndim, MPI_DOUBLE, **par, nwalker * ndim,
                MPI_DOUBLE, opts->mpi_comm);
#else
  for (int k = 0; k < nwalker; ++k)
    for (int d = 0; d < ndim; ++d) par[0][k][d] = oldp[k][d];
#endif

  // strech step
  unsigned int seed = time(NULL);
  double r = 1. * rand_r(&seed) / RAND_MAX;
  // *zz = sqr((opts->a - 1.) * r + 1.) / opts->a;
  double zz = ((opts->a - 1.) * r + 1.) * ((opts->a - 1.) * r + 1) / opts->a;

  int jwalker;                       // pick a waker
  int jrank = rand_r(&seed) % size;  // pick a rank
  if (jrank == irank) {              // same rank
    jwalker = rand_r(&seed) % (nwalker - 1);
    if (jwalker >= iwalker) jwalker++;  // avoid myself
  } else {                              // another rank
    jwalker = rand_r(&seed) % nwalker;
  }

  for (int d = 0; d < ndim; ++d)
    newp[d] = (1. - zz) * par[jrank][jwalker][d] + zz * par[irank][iwalker][d];

  FreeCArray(par);

  return zz;
}

void mcmc_walk_move(double *newp, double **oldp, int k, int nwalker, int np,
                    mcmc_opts *opts) {
  int *sub = new int[opts->p];
  double **cov;
  double *mean = new double[np], *diag = new double[np];

  NewCArray(cov, np, np);

  // draw unique sub-sample of size opts->p
  bool redraw;
  assert(opts->p < static_cast<int>(nwalker));
  assert(opts->p > static_cast<int>(np));

  unsigned int seed = time(NULL);

  for (int i = 0; i < opts->p; ++i) {
    do {
      int s = rand_r(&seed) % (nwalker - 1);
      if (s >= static_cast<int>(k)) s++;
      redraw = false;
      for (int j = 0; j < i; ++j)
        if (sub[j] == s) redraw = true;
      if (!redraw) sub[i] = s;
    } while (redraw);
  }

  // calculate mean
  for (size_t p = 0; p < np; ++p) {
    mean[p] = 0.;
    for (int q = 0; q < opts->p; ++q) mean[p] += oldp[sub[q]][p];
    mean[p] /= opts->p;
  }

  // calculate covariance matrix of this sample
  for (size_t i = 0; i < np; ++i) {
    for (size_t j = 0; j < np; ++j) {
      cov[i][j] = 0.;
      for (int q = 0; q < opts->p; ++q)
        cov[i][j] += (oldp[sub[q]][i] - mean[i]) * (oldp[sub[q]][j] - mean[j]);
      cov[i][j] /= opts->p;
    }
  }

  // Cholesky decomposition: cov = L * L'
  // L is stored in the lower triangle of A, with diag being the diagonal
  // elements
  _choldc(cov, np, diag);

  // Box-Muller method to generate standard multi-variate normal
  // https://en.wikipedia.org/wiki/Normal_distribution#Generating_values_from_normal_distribution
  // mean is reused to store the standard normal variable

  for (size_t i = 0; i < np; ++i) {
    double u = 1. * rand_r(&seed) / RAND_MAX;
    double v = 1. * rand_r(&seed) / RAND_MAX;

    mean[i] = sqrt(-2. * log(u)) * cos(2. * M_PI * v);
  }

  // multi-variate gaussian with cov: Y = P + L * X, x is a standard normal
  for (size_t i = 0; i < np; ++i) {
    newp[i] = oldp[k][i];
    for (size_t j = 0; j < i; ++j) newp[i] += cov[i][j] * mean[j];
    newp[i] += diag[i] * mean[i];
  }

  delete[] diag;
  delete[] mean;
  delete[] sub;

  FreeCArray(cov);
}
