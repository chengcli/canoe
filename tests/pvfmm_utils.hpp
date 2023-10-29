#ifndef _PVFMM_UTILS_HPP_
#define _PVFMM_UTILS_HPP_

#include <mpi.h>

enum DistribType { UnifGrid, RandUnif, RandGaus, RandElps, RandSphr };

template <class Real_t>
std::vector<Real_t> point_distrib(DistribType dist_type, size_t N,
                                  MPI_Comm comm) {
  int np, myrank;
  MPI_Comm_size(comm, &np);
  MPI_Comm_rank(comm, &myrank);
  static size_t seed = myrank + 1;
  seed += np;
  srand48(seed);

  std::vector<Real_t> coord;
  switch (dist_type) {
    case UnifGrid: {
      size_t NN = (size_t)round(pow((double)N, 1.0 / 3.0));
      size_t N_total = NN * NN * NN;
      size_t start = myrank * N_total / np;
      size_t end = (myrank + 1) * N_total / np;
      for (size_t i = start; i < end; i++) {
        coord.push_back((Real_t)(((i / 1) % NN) + 0.5) / NN);
        coord.push_back((Real_t)(((i / NN) % NN) + 0.5) / NN);
        coord.push_back((Real_t)(((i / (NN * NN)) % NN) + 0.5) / NN);
      }
    } break;
    case RandUnif: {
      size_t N_total = N;
      size_t start = myrank * N_total / np;
      size_t end = (myrank + 1) * N_total / np;
      size_t N_local = end - start;
      coord.resize(N_local * 3);

      for (size_t i = 0; i < N_local * 3; i++) coord[i] = ((Real_t)drand48());
    } break;
    case RandGaus: {
      Real_t e = sctl::const_e<Real_t>();
      Real_t log_e = sctl::log(e);
      size_t N_total = N;
      size_t start = myrank * N_total / np;
      size_t end = (myrank + 1) * N_total / np;

      for (size_t i = start * 3; i < end * 3; i++) {
        Real_t y = -1;
        while (y <= 0 || y >= 1) {
          Real_t r1 = (Real_t)(sqrt(-2 * log(drand48()) / log_e) *
                               cos(2 * M_PI * drand48()));
          Real_t r2 = (Real_t)pow(0.5, i * 10 / N_total);
          y = (Real_t)0.5 + r1 * r2;
        }
        coord.push_back(y);
      }
    } break;
    case RandElps: {
      size_t N_total = N;
      size_t start = myrank * N_total / np;
      size_t end = (myrank + 1) * N_total / np;
      size_t N_local = end - start;
      coord.resize(N_local * 3);

      const Real_t r = (Real_t)0.45;
      const Real_t center[3] = {0.5, 0.5, 0.5};
      for (size_t i = 0; i < N_local; i++) {
        Real_t* y = &coord[i * 3];
        Real_t phi = (Real_t)(2 * M_PI * drand48());
        Real_t theta = (Real_t)(M_PI * drand48());
        y[0] = center[0] + (Real_t)(0.25 * r * sin(theta) * cos(phi));
        y[1] = center[1] + (Real_t)(0.25 * r * sin(theta) * sin(phi));
        y[2] = center[2] + r * (Real_t)cos(theta);
      }
    } break;
    case RandSphr: {
      size_t N_total = N;
      size_t start = myrank * N_total / np;
      size_t end = (myrank + 1) * N_total / np;
      size_t N_local = end - start;
      coord.resize(N_local * 3);

      const Real_t center[3] = {0.5, 0.5, 0.5};
      for (size_t i = 0; i < N_local; i++) {
        Real_t* y = &coord[i * 3];
        Real_t r = 1;
        while (r > 0.5 || r == 0) {
          y[0] = (Real_t)drand48();
          y[1] = (Real_t)drand48();
          y[2] = (Real_t)drand48();
          r = (Real_t)sqrt((y[0] - center[0]) * (y[0] - center[0]) +
                           (y[1] - center[1]) * (y[1] - center[1]) +
                           (y[2] - center[2]) * (y[2] - center[2]));
          y[0] = center[0] + (Real_t)0.45 * (y[0] - center[0]) / r;
          y[1] = center[1] + (Real_t)0.45 * (y[1] - center[1]) / r;
          y[2] = center[2] + (Real_t)0.45 * (y[2] - center[2]) / r;
        }
      }
    } break;
    default:
      break;
  }
  return coord;
}

#endif  // _PVFMM_UTILS_HPP_
