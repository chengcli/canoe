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
#include <map>
#include <memory>
#include <string>
#include <vector>

// athena
#include <athena/athena.hpp>

// utils
#include <utils/parameter_map.hpp>

class AirParcel;

class Absorber {
 public:
  explicit Absorber(std::string name);

  Absorber(std::string name, std::vector<std::string> const& names,
           ParameterMap params);

  std::string GetName() const { return name_; }

  std::string GetCategory() const { return category_; }

  void SetModel(std::string name) { model_name_ = name; }

  virtual ~Absorber();

  virtual void LoadCoefficient(std::string fname, size_t bid);

  virtual Real GetAttenuation(Real wave1, Real wave2,
                              AirParcel const& var) const;

  virtual Real GetSingleScatteringAlbedo(Real wave1, Real wave2,
                                         AirParcel const& var) const;

  virtual void GetPhaseMomentum(Real* pp, Real wave1, Real wave2,
                                AirParcel const& var, int np) const;

 protected:
  //! absorber name
  std::string name_;

  //! absorber category
  std::string category_;

  //! absorption model model
  std::string model_name_;

  //! id of dependent molecules
  std::vector<int> imols_;

  //! parameters
  ParameterMap params_;
};

using AbsorberPtr = std::shared_ptr<Absorber>;

#endif  // SRC_HARP_ABSORBER_HPP_
