#ifndef SRC_OPACITY_NITROGEN_CIA_HPP_
#define SRC_OPACITY_NITROGEN_CIA_HPP_

// C/C++
#include <string>
#include <vector>

// harp
#include <harp/absorber.hpp>

class N2N2CIA : public Absorber {
 public:
  N2N2CIA(SpeciesNames const& species, ParameterMap params)
      : Absorber("N2-N2", species, params) {}
  virtual ~N2N2CIA() {}

  void LoadCoefficient(std::string fname, size_t bid) override;

  Real GetAttenuation(Real wave1, Real wave2,
                      Variable const& var) const override {
    return getAttenuation1((wave1 + wave2) / 2., var);
  }

 protected:
  Real getAttenuation1(Real wave, Variable const& var) const;

  size_t rt_len_[2];
  size_t fd1_len_[2];
  size_t fd2_len_[2];

  std::vector<Real> rt_axis_;
  std::vector<Real> fd1_axis_;
  std::vector<Real> fd2_axis_;
  std::vector<Real> rt_;
  std::vector<Real> fd1_;
  std::vector<Real> fd2_;
};

#endif  // SRC_OPACITY_NITROGEN_CIA_HPP_
