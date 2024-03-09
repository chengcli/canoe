#ifndef SRC_ASTRO_CELESTRIAL_BODY_HPP_
#define SRC_ASTRO_CELESTRIAL_BODY_HPP_

// C/C++
#include <memory>
#include <string>

// athena
#include <athena/athena.hpp>

// canoe
#include <virtual_groups.hpp>

class ParameterInput;
struct float_triplet;
struct Direction;
class Meshblock;

class CelestrialBody : public NamedGroup {
 public:
  // data
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
  CelestrialBody(ParameterInput *pin, std::string name);
  ~CelestrialBody();

  Direction ParentZenithAngle(Real time, Real colat, Real lon) const;
  Real ParentInsolationFlux(Real wav, Real dist_au) const;
  Real ParentInsolationFlux(Real wlo, Real whi, Real dist_au) const;
  Real ParentDistanceInAu(Real time) const;

  bool HasParentFlux() const {
    if (parent_) {
      return !parent_->spec_.empty();
    }
    return false;
  }

 protected:
  void loadOrbitalData(ParameterInput *pin);
  void loadSpectralData(std::string sfile);

 protected:
  //! pointer to parent body
  std::shared_ptr<CelestrialBody> parent_;

  //! emission spectra data
  std::vector<float_triplet> spec_;

 private:
  //! temporary search point
  mutable int il_ = -1;
};

using CelestrialBodyPtr = std::shared_ptr<CelestrialBody>;

Real GetGravity(char const *name, Real pclat);

class PlanetFactory {
 public:
  static CelestrialBodyPtr CreateFrom(MeshBlock *pmb, ParameterInput *pin);
};

#endif  // SRC_ASTRO_CELESTRIAL_BODY_HPP_
