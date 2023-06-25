/**@file
 * @brief This file contains declaration of Absorber
 *
 * **Author** : Cheng Li, California Institute of Technology <br>
 * **Contact** : cli@gps.caltech.edu <br>
 * **Revision history** :
 * - June 21 2016, start documenting this file
 * - July 28 2016, merge scatter into absorber
 * - June 24 2017, adapt to Athena++ framework
 * - April 03 2019, merge to snap
 * - July 27 2019, add multiple dependent molecules
 */

#ifndef SRC_HARP_ABSORBER_HPP_
#define SRC_HARP_ABSORBER_HPP_

// C/C++
#include <string>
#include <vector>

// athena
#include <athena/athena.hpp>

class CellVariables;

class Absorber {
 public:
  Absorber(std::string name);

  std::string GetName() const { return name_; }

  virtual ~Absorber();

  virtual void LoadCoefficient(std::string fname, size_t bid);

  virtual Real GetAttenuation(Real wave1, Real wave2,
                              CellVariables const& var) const;

  virtual Real GetSingleScatteringAlbedo(Real wave1, Real wave2,
                                         CellVariables const& var) const;

  virtual void GetPhaseMomentum(Real* pp, Real wave1, Real wave2,
                                CellVariables const& var, int np) const;

 protected:
  // data
  std::string name_;

  //! id of dependent molecules
  std::vector<int> imols_;

  //! mixr of dependent molecules
  std::vector<Real> mixrs_;
};

#endif  // SRC_HARP_ABSORBER_HPP_
