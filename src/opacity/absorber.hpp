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
class Absorber : public NamedGroup,
                 public ParameterGroup,
                 public SpeciesIndexGroup,
                 public StringReprGroup {
 public:  // constructor and destructor
  Absorber(std::string name);
  virtual ~Absorber();

 public:  // member functions
  //! Set absorption model
  void SetModel(std::string name) { model_name_ = name; }

  //! Combines SetOpacityFile() and LoadOpacity()
  void LoadOpacityFromFile(std::string filename);

  //! Set opacity filename to internal variable, does not load opacity
  void SetOpacityFile(std::string filename);

  //! Load opacity from internal variable
  void LoadOpacity(int bid);

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

  virtual void CheckFail() const {}

 public:  // StringRepr
  std::string ToString() const override;

 protected:
  //! absorption model model
  std::string model_name_;

  //! opacity filename
  std::string opacity_filename_;
};

using AbsorberPtr = std::shared_ptr<Absorber>;
using AbsorberContainer = std::vector<AbsorberPtr>;

class AbsorberFactory {
 public:
  //! \todo make it a static member of Absorber
  //! search path for radiation input file
  std::string search_path;

  //! \brief Create an absorber from YAML node
  //!
  //! \param[in] my YAML node containing the current absorber
  //! \param[in] band_name name of the radiation band
  static AbsorberPtr CreateFrom(YAML::Node const& my, std::string band_name);

  //! \brief Create absorbers from YAML node
  //!
  //! \param[in] names names of absorbers
  //! \param[in] band_name name of the radiation band
  //! \param[in] rad YAML node containing the radiation input file
  static AbsorberContainer CreateFrom(std::vector<std::string> const& names,
                                      std::string band_name,
                                      YAML::Node const& rad);

 protected:
  //! \brief Only create an absorber based on its name and class
  //!
  //! \param[in] name name of the absorber
  //! \param[in] type class identifier of the absorber
  static AbsorberPtr createAbsorberPartial(std::string name, std::string type);
};

#endif  // SRC_OPACITY_ABSORBER_HPP_
