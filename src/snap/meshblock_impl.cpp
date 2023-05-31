// athena
#include <athena/athena.hpp>
#include <athena/parameter_input.hpp>

// canoe
#include <configure.hpp>

// harp
#include <harp/radiation.hpp>
#include <harp/radiation_band.hpp>

// inversion
#include <inversion/inversion.hpp>
#include <inversion/inversion_helper.hpp>

// snap
#include "decomposition/decomposition.hpp"
#include "implicit/implicit_solver.hpp"
#include "meshblock_impl.hpp"
#include "reconstruct/face_reconstruct.hpp"
#include "thermodynamics/thermodynamics.hpp"

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
  reference_pressure_ = pin->GetReal("mesh", "ReferencePressure");

  // pressure scale height
  pressure_scale_height_ = pin->GetReal("mesh", "PressureScaleHeight");
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
