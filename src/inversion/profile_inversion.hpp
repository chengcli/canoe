#ifndef PROFILE_INVERSION_HPP
#define PROFILE_INVERSION_HPP

// C/C++ headers
#include <vector>

// harp2 headers
#include <configure.hpp>
#include "inversion.hpp"

class ProfileInversion : public Inversion {
public:
  ProfileInversion(MeshBlock *pmb, ParameterInput *pin, std::string name);
  ~ProfileInversion();

  size_t samples() const
  {
    return plevel_.size() - 2;
  }

  void InitializePositions() override;

  void UpdateHydro(Hydro *phydro, ParameterInput *pin) const override;

	Real LogPosteriorProbability(Radiation *prad, Hydro *phydro,
    Real const *par, Real *val, int k) const override;

	void UpdateProfiles(Hydro *phydro, Real **XpSample, int k, int jl, int ju) const;

  void ConvectiveAdjustment(Hydro *phydro, int k, int ju) const;

	virtual Real LogPriorProbability(Real **XpSample) const;

  int getX2Span() const override
  {
    return idx_.size() + 1;
  }

protected:
  // pressure levels
	std::vector<Real>   plevel_;

  // inversion variable id
	std::vector<int>    idx_;

	// hyper-parameters
	Real                chi_;
	Real                Xstd_[1+NVAPOR];
	Real                Xlen_[1+NVAPOR];
};

class VLAProfileInversion : public ProfileInversion {
public:
  // functions
  VLAProfileInversion(MeshBlock *pmb, ParameterInput *pin):
    ProfileInversion(pmb, pin, "profile")
  {}

  ~VLAProfileInversion() {}

  void CalculateFitTarget(Radiation const *prad,
    Real *val, int nvalue, int k, int j) const override;
};

class JunoProfileInversion : public ProfileInversion {
public:
  // functions
  JunoProfileInversion(MeshBlock *pmb, ParameterInput *pin):
    ProfileInversion(pmb, pin, "profile")
  {}

  ~JunoProfileInversion() {}

  void CalculateFitTarget(Radiation const *prad,
    Real *val, int nvalue, int k, int j) const override;
};

#endif
