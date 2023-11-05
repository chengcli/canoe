/**@file
 * @brief This file contains declaration of Absorber
 *
 * **Author** : Cheng Li, California Institute of Technology <br>
 * **Contact** : cli@gps.caltech.edu <br>
 * **Revision history** :
 * - June 21 2016, start documenting this file
 * - July 28 2016, merge scatter into absorber
 * - June 24 2017, adapt to Athena++ framework
 * - Aug 1 2019, add annotation and protection
 */

#ifndef SRC_OPACITY_MICROWAVE_MWR_ABSORBERS_HPP_
#define SRC_OPACITY_MICROWAVE_MWR_ABSORBERS_HPP_

// C++
#include <string>
#include <vector>

// opacity
#include <opacity/absorber.hpp>

class MeshBlock;

namespace GiantPlanets {

class MwrAbsorberCIA : public Absorber {
 public:
  MwrAbsorberCIA() : Absorber("CIA") {}

  void CheckFail() const override;
  Real GetAttenuation(Real wave1, Real wave2,
                      AirParcel const& var) const override;
};

class MwrAbsorberNH3 : public Absorber {
 public:
  MwrAbsorberNH3() : Absorber("NH3") {}

  MwrAbsorberNH3& SetModelHanley() {
    model_name_ = "Hanley09";
    return *this;
  }
  MwrAbsorberNH3& SetModelBellotti() {
    model_name_ = "Bellotti16";
    return *this;
  }
  MwrAbsorberNH3& SetModelBellottiSwitch() {
    model_name_ = "BellottiSwitch16";
    return *this;
  }
  MwrAbsorberNH3& SetModelDevaraj() {
    model_name_ = "Devaraj";
    return *this;
  }
  MwrAbsorberNH3& SetModelRadtran() {
    model_name_ = "Radtran";
    return *this;
  }

  void CheckFail() const override;
  Real GetAttenuation(Real wave1, Real wave2,
                      AirParcel const& var) const override;
};

class MwrAbsorberPH3 : public Absorber {
 public:
  MwrAbsorberPH3() : Absorber("PH3") {}

  MwrAbsorberPH3& SetModelRadtran() {
    model_name_ = "Radtran";
    return *this;
  }
  MwrAbsorberPH3& SetModelHoffman() {
    model_name_ = "Hoffman";
    return *this;
  }

  Real GetAttenuation(Real wave1, Real wave2, AirParcel const& var) const;
};

class MwrAbsorberH2O : public Absorber {
 public:
  //! \todo check Karpowics model
  MwrAbsorberH2O() : Absorber("H2O") {}

  MwrAbsorberH2O& SetModeldeBoer() {
    model_name_ = "deBoer";
    return *this;
  }
  MwrAbsorberH2O& SetModelWaters() {
    model_name_ = "Waters";
    return *this;
  }
  MwrAbsorberH2O& SetModelGoodman() {
    model_name_ = "Goodman";
    return *this;
  }
  MwrAbsorberH2O& SetModelKarpowicz() {
    model_name_ = "Karpowicz";
    return *this;
  }

  void CheckFail() const override;
  Real GetAttenuation(Real wave1, Real wave2,
                      AirParcel const& var) const override;
};

class MwrAbsorberH2S : public Absorber {
 public:
  MwrAbsorberH2S() : Absorber("H2S") {}

  Real GetAttenuation(Real wave1, Real wave2, AirParcel const& var) const;
};

class MwrAbsorberElectron : public Absorber {
 public:
  MwrAbsorberElectron() : Absorber("Electron") {}

  MwrAbsorberElectron& SetModelAppletonHartree() {
    model_name_ = "AppletonHartree";
    return *this;
  }

  MwrAbsorberElectron& SetModelChengLi() {
    model_name_ = "ChengLi";
    return *this;
  }
  MwrAbsorberElectron& SetModelReference() {
    model_name_ = "Reference";
    return *this;
  }

  Real GetAttenuation(Real wave1, Real wave2, AirParcel const& var) const;
};

}  // namespace GiantPlanets

#endif  // SRC_OPACITY_MICROWAVE_MWR_ABSORBERS_H_
