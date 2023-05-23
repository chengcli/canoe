// Athena++ headers
#include "meshblock_impl.hpp"

#include <athena/athena.hpp>
#include <athena/parameter_input.hpp>
#include <configure.hpp>
#include <harp/radiation.hpp>
#include <harp/radiation_band.hpp>
#include <inversion/inversion.hpp>

#include "../hydro/decomposition/decomposition.hpp"
#include "../hydro/implicit/implicit_solver.hpp"
#include "../reconstruct/face_reconstruct.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "block_index.hpp"

MeshBlock::Impl::Impl(MeshBlock *pmb, ParameterInput *pin) : pmy_block_(pmb) {
  du.NewAthenaArray(NumHydros, pmb->ncells3, pmb->ncells2, pmb->ncells1);

  // block index
  pblock = new BlockIndex(pmb);

  // thermodynamics
  pthermo = new Thermodynamics(pmb, pin);

  // hydro decomposition
  pdec = new Decomposition(pmb->phydro);

  // implicit methods
  phevi = new ImplicitSolver(pmb, pin);

  // reconstruction
  precon = new FaceReconstruct(pmb, pin);

  // radiation
  prad = new Radiation(pmb, pin);

  // inversion queue
  new_inversion_queue(fitq, pmb, pin);
}

// Athena++ demands destruct pbval AFTER all boundary values
// But in this mod, boundary values are destructed BEFORE pbval
// TODO, check if this is OK
MeshBlock::Impl::~Impl() {
  // std::cout << "Impl desctructor" << std::endl;
  delete pblock;
  delete pthermo;
  delete pdec;
  delete phevi;
  delete precon;

  delete prad;
  for (auto q : fitq) {
    delete q;
  }
}
