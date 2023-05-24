#ifndef SRC_HARP_ABSORBER_HPP_
#define SRC_HARP_ABSORBER_HPP_

// C++ headers
#include <athena/athena.hpp>
#include <configure.hpp>
#include <string>
#include <vector>

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

class MeshBlock;
class ParameterInput;
class Thermodynamics;
class CellVariables;

class Absorber {
 public:
  Absorber(MeshBlock* pmb, ParameterInput* pin, std::string bname,
           std::string name);

  virtual ~Absorber();

  virtual void loadCoefficient(std::string fname, int bid);

  virtual Real getAttenuation(Real wave1, Real wave2,
                              CellVariables const& var) const;

  virtual Real getSingleScatteringAlbedo(Real wave1, Real wave2,
                                         CellVariables const& var) const;

  virtual void getPhaseMomentum(Real* pp, Real wave1, Real wave2,
                                CellVariables const& var, int np) const;

  Absorber* use(Thermodynamics const* p) {
    pthermo_ = p;
    return this;
  }

 protected:
  // data
  std::string name_;

  /**< id of dependent molecules */
  std::vector<int> imols_;

  /**< mixr of dependent molecules */
  std::vector<Real> mixrs_;

  // connections
  Thermodynamics const* pthermo_;
};

#endif  // SRC_HARP_ABSORBER_HPP_
