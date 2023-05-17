#ifndef RADIATION_UTILS_HPP
#define RADIATION_UTILS_HPP

// C/C++ header
#include <string>

#include <athena/athena.hpp>

#include "radiation.hpp"


/*void WriteTopFlux(std::string fname) const;
void WriteTopRadiance(std::string fname) const;
void WriteOpticalDepth(std::string fname) const;
void WriteHeatingRate(std::string fname, AthenaArray<Real> const& flux,
      AthenaArray<Real> const& hr, Real const* level); */

void read_radiation_directions(std::vector<Direction>& ray, std::string str);
void set_radiation_flags(uint64_t *flags, std::string str);

void getPhaseHenyeyGreenstein(Real *pmom, int iphas, Real gg, int npmom);
void packSpectralProperties(Real *buf, Real const *tau, Real const *ssa, Real const* pmom, int nlayer, int npmom);
void unpackSpectralProperties(Real *tau, Real *ssa, Real *pmom, Real const *buf, int slyr, int npmom, int nblocks, int npmom_max = -1);

#endif
