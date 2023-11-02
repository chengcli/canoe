#ifndef SRC_ASTRO_CELESTRIAL_BODY_HPP_
#define SRC_ASTRO_CELESTRIAL_BODY_HPP_

// C/C++
#include <memory>
#include <string>

// athena
#include <athena/athena.hpp>

class ParameterInput;
struct float_triplet;
struct Direction;

class CelestrialBody {
 public:
  // data
  CelestrialBody *parent;
  std::string name;
  Real re;       // equatorial radius [km -> m]
  Real rp;       // polar radius [km -> m]
  Real obliq;    // obliquity [deg -> rad]
  Real spinp;    // spin period [day -> s]
  Real orbit_a;  // orbital semi-major axis [au -> m]
  Real orbit_e;  // orbital eccentricity [1]
  Real orbit_i;  // orbital inclination to ecliptic [deg -> rad]
  Real orbit_p;  // orbital period [day -> s]
  Real equinox;  // vernal equinox
  Real grav_eq;  // equatorial gravity at surface [m/s^2]

  // functions
  explicit CelestrialBody(ParameterInput *pin);
  CelestrialBody(ParameterInput *pin, std::string myname);
  ~CelestrialBody();

  void ReadSpectraFile(std::string sfile);
  Direction ParentZenithAngle(Real time, Real colat, Real lon);
  Real ParentInsolationFlux(Real wav, Real dist_au);
  Real ParentInsolationFlux(Real wlo, Real whi, Real dist_au);
  Real ParentDistanceInAu(Real time);

 protected:
  void readCelestrialData(ParameterInput *pin, std::string myname);

  // emission spectra data
  float_triplet *spec_;
  int nspec_;
  int il_;  // search pointer
};

using CelestrialBodyPtr = std::shared_ptr<CelestrialBody>;

Real GetGravity(char const *name, Real pclat);

#endif  // SRC_ASTRO_CELESTRIAL_BODY_HPP_
