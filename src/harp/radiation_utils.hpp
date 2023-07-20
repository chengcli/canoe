#ifndef SRC_HARP_RADIATION_UTILS_HPP_
#define SRC_HARP_RADIATION_UTILS_HPP_

// C/C++ header
#include <athena/athena.hpp>
#include <string>
#include <vector>

#include "radiation.hpp"

/*void WriteTopFlux(std::string fname) const;
void WriteTopRadiance(std::string fname) const;
void WriteOpticalDepth(std::string fname) const;
void WriteHeatingRate(std::string fname, AthenaArray<Real> const& flux,
      AthenaArray<Real> const& hr, Real const* level); */

void read_radiation_directions(std::vector<Direction> *ray, std::string str);
void set_radiation_flags(uint64_t *flags, std::string str);

void get_phase_momentum(Real *pmom, int iphas, Real gg, int npmom);
void packSpectralProperties(Real *buf, Real const *tau, Real const *ssa,
                            Real const *pmom, int nlayer, int npmom);
void unpackSpectralProperties(Real *tau, Real *ssa, Real *pmom, Real const *buf,
                              int slyr, int npmom, int nblocks,
                              int npmom_max = -1);

#endif  // SRC_HARP_RADIATION_UTILS_HPP_
