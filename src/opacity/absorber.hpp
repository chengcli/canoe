#ifndef SRC_OPACITY_ABSORBER_HPP_
#define SRC_OPACITY_ABSORBER_HPP_

// C/C++
#include <map>
#include <memory>
#include <string>
#include <vector>

// external
#include <yaml-cpp/yaml.h>

// athena
#include <athena/athena.hpp>

// canoe
#include <virtual_groups.hpp>

class AirParcel;

//! \brief base class of all absorbers
class Absorber : public NamedGroup, public ParameterGroup {
 public:  // constructor and destructor
  Absorber(std::string name,
           std::vector<std::string> const& names) virtual ~Absorber();

 public:  // functions
  //! Set absorption model
  void SetModel(std::string name) { model_name_ = name; }

  //! Load absorption coefficient from file
  virtual void LoadCoefficient(std::string fname, size_t bid) {}

  //! Get attenuation coefficient [1/m]
  virtual Real GetAttenuation(Real wave1, Real wave2,
                              AirParcel const& var) const {
    return 0.;
  }

  //! Get single scattering albedo [1]
  virtual Real GetSingleScatteringAlbedo(Real wave1, Real wave2,
                                         AirParcel const& var) const {
    return 0.;
  }

  //! Get phase function [1]
  virtual void GetPhaseMomentum(Real* pp, Real wave1, Real wave2,
                                AirParcel const& var, int np) const {}

 protected:
  //! absorption model model
  std::string model_name_;

  //! id of dependent molecules
  std::vector<int> imols_;
};

using AbsorberPtr = std::shared_ptr<Absorber>;

class AbsorberFactory {
 public:
  //! Create an absorber from YAML node
  //!
  //! \param[in] my YAML node containing the current absorber
  static AbsorberPtr CreateFrom(YAML::Node& my);

  //! \brief Create absorbers from YAML node
  //!
  //! \param[in] absorber_names names of absorbers
  //! \param[in] rad YAML node containing the whole radiation input file
  static std::vector<std::shared_ptr<Absorber>> CreateFrom(
      std::vector<std::string> const& absorber_names, YAML::Node& rad);
};

#endif  // SRC_HARP_ABSORBER_HPP_
