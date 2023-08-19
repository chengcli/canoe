#ifndef SRC_OPACITY_OXYGEN_CIA_HPP_
#define SRC_OPACITY_OXYGEN_CIA_HPP_

// C/C++
#include <vector>

// harp
#include <harp/absorber.hpp>

class O2O2CIA : public Absorber {
 public:
  O2O2CIA(std::vector<std::string> const& species, ParameterMap params)
      : Absorber("O2-O2", species, params) {
    category_ = "cia";
  }

  virtual ~O2O2CIA() {}

  void LoadCoefficient(std::string fname, size_t bid) override;

  Real GetAttenuation(Real wave1, Real wave2,
                      AirParcel const& var) const override {
    return getAttenuation1((wave1 + wave2) / 2., var);
  }

 protected:
  Real getAttenuation1(Real wave, AirParcel const& var) const;

  size_t fd_len_[2];
  size_t a1dg_x3sg00_len_[2];
  size_t a1dg_x3sg10_len_[2];
  size_t ab_len_[2];
  size_t other_len_[2];

  std::vector<Real> fd_axis_;
  std::vector<Real> a1dg_x3sg00_axis_;
  std::vector<Real> a1dg_x3sg10_axis_;
  std::vector<Real> ab_axis_;
  std::vector<Real> other_axis_;
  std::vector<Real> fd_;
  std::vector<Real> a1dg_x3sg00_;
  std::vector<Real> a1dg_x3sg10_;
  std::vector<Real> ab_;
  std::vector<Real> other_;
};

#endif  //  SRC_OPACITY_OXYGEN_CIA_HPP_
