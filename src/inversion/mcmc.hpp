#ifndef MCMC_HPP_
#define MCMC_HPP_
#include <cstddef>

// Athena++ headers
//#include "../defs.hpp"

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif


/**@file
 * @brief This file contains data structure and subroutines implementing Markov
 * Chain Monte Carlo sampler
 *
 * **Author** : Cheng Li, California Institute of Technology <br>
 * **Contact** : cli@gps.caltech.edu <br>
 * **Revision history** : <br>
 * - June 21 2015, initial document
 * - July 18 2016, remove MathHelper dependence
 */

/** Given a positive definite symmetric matrix a[0..n-1][0..n-1], this
 * subroutine constructs its Cholesky decomposition: A = L*L'. On input, only
 * the upper triangle of A need be given; it is not modified. The Cholesky
 * factor L is returned in the lower triangle of A, except for its diagonal
 * elements which are returned in p[0..n-1].
 * see Numerical Recipes in C, pg 97
 */
void _choldc(
        double **a, 
        int n, 
        double *p
        );

/** calculate auto-correlation time given a time series 
 * Written by Jonathan Goodman, March 2009, goodman@cims.nyu.edu
 * Modified by Cheng Li, June 2016, cli@gps.caltech.edu
 */
int _acor(
        double *mean,
        double *sigma,
        double *tau,
        double *X,
        int L
        );

/** mcmc options */
struct mcmc_opts {
  double  a, b, c, d;
  int     p, q, r;

  int     print;
  char    logfile[80];

#ifdef MPI_PARALLEL
  MPI_Comm    mpi_comm;
#endif
};

/** mcmc chains and records */
struct mcmc_recs {
  int    ndim,       /**< dimension of state vector */
         nvalue,     /**< dimension of measurement  */
         nwalker,    /**< number of walkers */
         nstep;      /**< number of steps */

  double ***par,      /**< state vector in each chain */
         ***val,      /**< forward result in each chain */
         **lnp;      /**< log probability of each state */

  double *opt_par,   /**< optimal parameter */
         *opt_val;   /**< optimal forward result */

  int    **newstate; /**< 0/1 whether the current state is a new state */

  int    cur,        /**< current chain length */
         accept,     /**< number of accepted states */
         reset;      /**< reset chain length */
  double opt_lnp;    /**< optimal log probability */
};

/** allocate memory for mcmc chains */
void mcmc_alloc(
  mcmc_recs *recs,
  int nstep,
  int nwalker,
  int ndim,
  int nvalue
);

/** free memory */
void mcmc_free(
  mcmc_recs *recs
);

/** calculate autocorrelation time */
void mcmc_statistics(
  double *mean,
  double *sigma,
  double *tau,
  mcmc_recs *recs
);

/** save mcmc chains to a FITS file */
void mcmc_save_fits(
  char const *fname,
  mcmc_opts *opts, 
  mcmc_recs *recs,
  int include_last = false
);

/** load mcmc chains from a FITS file */
void mcmc_load_fits(
  char const *fname,
  mcmc_opts *opts,
  mcmc_recs *recs,
  int alloc = true
);

/** copy mcmc chains */
void mcmc_append_recs(
  mcmc_recs *dst,
  mcmc_recs *src
);

// report MCMC sampler result
void mcmc_report(
  mcmc_opts *opts,
  mcmc_recs *recs,
  char const *mode
);

typedef double (*ObjectiveFunction_t)(double *, double*, int, int, void*);

void mcmc_init(
  ObjectiveFunction_t lnprob,
  double **par,
  mcmc_opts *opts,
  mcmc_recs *recs,
  void *obj
);

/** implement Ensemble Samplers with Affine Invariance, Goodman and Weare (2010)
 * @param lnprob I, function calculating log probabiligy
 * @param par I/O, parameter vector
 * @param val I, measurement vector
 * @param np I, size of input parameter vector
 * @param nv I, size of measurement vector
 * @param nwalker I, number of walkers
 * @param nstep I, number of steps
 * @param opts I, mcmc options
 * @param recs O, mcmc chains
 */
void mcmc_advance(
  ObjectiveFunction_t lnprob,
  mcmc_opts *opts,
  mcmc_recs *recs,
  void *obj
);

double mcmc_stretch_move(
  double *newp,
  double **oldp,
  int iwalker,
  int nwalker,
  int ndim,
  mcmc_opts *opts
);

void mcmc_walk_move(
  double *newp,
  double **oldp,
  int k,
  int nwalker,
  int np,
  mcmc_opts *opts
);

#endif
