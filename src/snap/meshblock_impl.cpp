// canoe
#include <configure.hpp>

// athena
#include <athena/athena.hpp>
#include <athena/parameter_input.hpp>

// harp
#include <harp/radiation.hpp>
#include <harp/radiation_band.hpp>

// application
#include <application/exceptions.hpp>

// inversion
#include <inversion/inversion.hpp>
#include <inversion/inversion_helper.hpp>

// utils
#include <utils/vectorize.hpp>

// snap
#include "decomposition/decomposition.hpp"
#include "implicit/implicit_solver.hpp"
#include "meshblock_impl.hpp"
#include "reconstruct/face_reconstruct.hpp"
#include "thermodynamics/thermodynamics.hpp"

MeshBlock::IndexMap::IndexMap(MeshBlock *pmb, ParameterInput *pin)
    : pmy_block_(pmb) {
  // species id
  std::string str = pin->GetOrAddString("species", "vapor", "");
  std::vector<std::string> names = Vectorize<std::string>(str.c_str(), " ,");

  for (size_t i = 0; i < names.size(); ++i) {
    vapor_index_map_[names[i]] = 1 + i;
  }

  str = pin->GetOrAddString("species", "tracer", "");
  names = Vectorize<std::string>(str.c_str(), " ,");
  for (size_t i = 0; i < names.size(); ++i) {
    tracer_index_map_[names[i]] = i;
  }

  str = pin->GetOrAddString("species", "cloud", "");
  names = Vectorize<std::string>(str.c_str(), " ,");
  for (size_t i = 0; i < names.size(); ++i) {
    cloud_index_map_[names[i]] = i;
  }

  str = pin->GetOrAddString("species", "chemistry", "");
  names = Vectorize<std::string>(str.c_str(), " ,");
  for (size_t i = 0; i < names.size(); ++i) {
    chemistry_index_map_[names[i]] = i;
  }

  str = pin->GetOrAddString("species", "particle", "");
  names = Vectorize<std::string>(str.c_str(), " ,");
  for (size_t i = 0; i < names.size(); ++i) {
    particle_index_map_[names[i]] = i;
  }

  str = pin->GetOrAddString("species", "static", "");
  names = Vectorize<std::string>(str.c_str(), " ,");
  for (size_t i = 0; i < names.size(); ++i) {
    static_index_map_[names[i]] = i;
  }
}

size_t MeshBlock::IndexMap::GetSpeciesId(std::string category_name) const {
  std::string delimiter = ".";

  // Find the position of the delimiter
  size_t delimiter_pos = category_name.find(delimiter);

  if (delimiter_pos == std::string::npos) {
    throw NotFoundError("GetSpeciesId",
                        "Delimiter '" + delimiter + "' in " + category_name);
  }

  // Extract the substrings
  std::string category = category_name.substr(0, delimiter_pos);
  std::string name = category_name.substr(delimiter_pos + delimiter.length());

  if (category == "vapor") {
    return GetVaporId(name);
  } else if (category == "tracer") {
    return NHYDRO + GetTracerId(name);
  } else if (category == "cloud") {
    return NHYDRO + NSCALARS + GetCloudId(name);
  } else if (category == "chemistry") {
    return NHYDRO + NSCALARS + NCLOUDS + GetChemistryId(name);
  } else {
    throw NotFoundError("GetSpeciesId", "Category " + category);
  }
}

MeshBlock::Impl::Impl(MeshBlock *pmb, ParameterInput *pin) : pmy_block_(pmb) {
  du.NewAthenaArray(NHYDRO, pmb->ncells3, pmb->ncells2, pmb->ncells1);

  // thermodynamics
  pthermo = new Thermodynamics(pmb, pin);

  // decomposition
  pdec = new Decomposition(pmb);

  // implicit methodsphydro
  phevi = new ImplicitSolver(pmb, pin);

  // reconstruction
  precon = new FaceReconstruct(pmb, pin);

  // radiation
  prad = new Radiation(pmb, pin);

  // inversion queue
  fitq = create_inversion_queue(pmb, pin);

  // reference pressure
#ifdef HYDROSTATIC
  reference_pressure_ = pin->GetReal("mesh", "ReferencePressure");

  // pressure scale height
  pressure_scale_height_ = pin->GetReal("mesh", "PressureScaleHeight");
#else
  reference_pressure_ = 1.0;
  pressure_scale_height_ = 1.0;
#endif  // HYDROSTATIC
}

// Athena++ demands destruct pbval AFTER all boundary values
// But in this mod, boundary values are destructed BEFORE pbval
// TODO(cli) check if this is OK
MeshBlock::Impl::~Impl() {
  // std::cout << "Impl desctructor" << std::endl;
  delete pthermo;
  delete pdec;
  delete phevi;
  delete precon;

  delete prad;
  for (auto q : fitq) {
    delete q;
  }
}
