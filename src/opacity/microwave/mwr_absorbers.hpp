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

#ifndef MWR_ABSORBERS_HPP_
#define MWR_ABSORBERS_HPP_

// C++
#include <string>
#include <vector>

// athena
//#include "../../athena.hpp"

// harp
#include <harp/absorber.hpp>

class MwrAbsorberCIA: public Absorber {
public:
  /** Construction function for CIA absorption
    * @param xHe molar mixing ratio of He
    * @param xCH4 molar mixing ratio of CH4
    * @param fequal mixing ratio equilibrium hydrogen
    */
  MwrAbsorberCIA(MeshBlock *pmb, ParameterInput *pin,
    Real xHe, Real xCH4, Real fequal = 0.);

  Real GetAttenuation(Real wave1, Real wave2,
      CellVariables const& var) const;

private:
  Real xHe_, xCH4_, fequal_;
};

class MwrAbsorberNH3: public Absorber {
public:
  /** Constructor for NH3 absorption (method == 1)
    * @param imols index of NH3 molecule
    * @param xHe molar mixing ratio of He
    * @param xH2O molar mixing ratio of H2O
    */
  MwrAbsorberNH3(MeshBlock *pmb, ParameterInput *pin,
    int imol, Real xHe, Real xH2O);

  /** Constructor for NH3 absorption (method == 2)
    * @param imols index of NH3 molecule
    * @param xHe molar mixing ratio of He
    * @param xH2O molar mixing ratio profile of H2O
    * @param pres pressure levels of xH2O
    */
  MwrAbsorberNH3(MeshBlock *pmb, ParameterInput *pin,
    int imol, Real xHe, Real *xH2O, Real *pres, int np);

  /** Constructor for NH3 absorption (method == 3)
    * @param imols an array of length 2, stores the index of NH3 and H2O
    * @param xHe molar mixing ratio of He
    */
  MwrAbsorberNH3(MeshBlock *pmb, ParameterInput *pin,
    std::vector<int> imols, Real xHe, Real power = 0.);

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

  Real GetAttenuation(Real wave1, Real wave2,
      CellVariables const& var) const;

private:
  std::string model_name_;
  int method_;
  Real xHe_, xH2O_;
  std::vector<Real> ref_xh2o_, ref_pres_;
  Real power_;
};

class MwrAbsorberPH3: public Absorber {
public:
  /** Constructor for PH3 absorption (method == 1)
    * @param imol index of PH3 molecule
    * @param xHe molar mixing ratio of He
    */
  MwrAbsorberPH3(MeshBlock *pmb, ParameterInput *pin,
    int imol, Real xHe);

  /** Constructor for PH3 absorption (method == 2)
    * @param xHe molar mixing ratio of He
    * @param xPH3 molar mixing ratio profle of PH3
    * @param pres reference pressure level
    * @param np number of elements in the array
    */
  MwrAbsorberPH3(MeshBlock *pmb, ParameterInput *pin,
    Real xHe, Real *xPH3, Real *pres, int np);

  MwrAbsorberPH3& SetModelRadtran() {
    model_name_ = "Radtran";
    return *this;
  }
  MwrAbsorberPH3& SetModelHoffman() {
    model_name_ = "Hoffman";
    return *this;
  }

  Real GetAttenuation(Real wave1, Real wave2,
      CellVariables const& var) const;
private:
  std::string model_name_;
  int method_;
  Real xHe_;
  std::vector<Real> ref_xph3_, ref_pres_;
};

class MwrAbsorberH2O: public Absorber {
public:
  /** Constructor for H2O absorption
    * @param xHe molar mixing ratio of He
    * @param imol index of H2O molecule
    */
  MwrAbsorberH2O(MeshBlock *pmb, ParameterInput *pin,
    int imol, Real xHe, Real scale = 0.);

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

  Real GetAttenuation(Real wave1, Real wave2,
      CellVariables const& var) const;
private:
  std::string model_name_;
  Real xHe_, scale_;
};

class MwrAbsorberH2S: public Absorber {
public:
  /** Constructor for H2O absorption (method == 1)
    * @param imol index of PH3 molecule
    * @param xHe molar mixing ratio of He
    */
  MwrAbsorberH2S(MeshBlock *pmb, ParameterInput *pin,
    int imol, Real xHe);

  /** Constructor for PH3 absorption (method == 2)
    * @param xHe molar mixing ratio of He
    * @param xH2S molar mixing ratio profle of PH3
    * @param pres reference pressure level
    * @param np number of elements in the array
    */
  MwrAbsorberH2S(MeshBlock *pmb, ParameterInput *pin,
    Real xHe, Real *xH2S, Real *pres, int np);

  Real GetAttenuation(Real wave1, Real wave2,
      CellVariables const& var) const;

private:
  int method_;
  Real xHe_;
  std::vector<Real> ref_xh2s_, ref_pres_;
};

class MwrAbsorberElectron: public Absorber {
public:
  MwrAbsorberElectron(MeshBlock *pmb, ParameterInput *pin, int ion);

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

  Real GetAttenuation(Real wave1, Real wave2,
      CellVariables const& var) const;

private:
  std::string model_name_;
};

#endif
