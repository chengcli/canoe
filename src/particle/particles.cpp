/** @file particles.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Sunday May 30, 2021 13:10:12 PDT
 * @bug No known bugs.
 */

// C++ headers
#include <ctime>
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include "../coordinates/coordinates.hpp"
#include "../debugger/debugger.hpp"
#include "../mesh/mesh.hpp"
#include "particle_buffer.hpp"
#include "particles.hpp"

Particles::Particles(MeshBlock *pmb, ParameterInput *pin)
    : pmy_block(pmb),
      myname("HEAD"),
      prev(nullptr),
      next(nullptr),
      seeds_per_cell_(1),
      nmax_per_cell_(1),
      density_floor_(0) {
  pmb->pdebug->Enter("Particle List");
  std::stringstream &msg = pmb->pdebug->msg;

  ppb = new ParticleBuffer(this);

  char particle_names[1024], *p;
  std::string str = pin->GetOrAddString("particles", "particles", "");
  std::strcpy(particle_names, str.c_str());
  p = std::strtok(particle_names, " ,");

  int npart = 0;
  while (p != NULL) {
    Particles *pnew;
    std::string name;
    char *c = std::strchr(p, '.');
    if (c != NULL)
      name = c + 1;
    else
      name = p;
    if (std::strncmp(p, "2pcp", 4) == 0) {
      pnew = AddParticles(TwoPhaseCloudParticles(pmb, pin, name));
    } else if (std::strncmp(p, "scp", 3) == 0) {
      pnew = AddParticles(SimpleCloudParticles(pmb, pin, name));
    } else {
      msg << "### FATAL ERROR in function Particles::Particles" << std::endl
          << "Particles '" << p << "' "
          << "unrecognized." << std::endl;
      ATHENA_ERROR(msg);
    }
    p = std::strtok(NULL, " ,");
    npart += pnew->u.GetDim4();
  }

  pmb->pdebug->Leave();
}

// constructor, initializes data structure and parameters
Particles::Particles(MeshBlock *pmb, ParameterInput *pin, std::string name,
                     int nct)
    : pmy_block(pmb), myname(name), prev(nullptr), next(nullptr) {
  pmb->pdebug->Enter("Basic Particles");
  std::stringstream &msg = pmb->pdebug->msg;
  msg << "- " << name << " particle categories = " << nct << std::endl;
  ppb = new ParticleBuffer(this);
  int nc1 = pmb->ncells1, nc2 = pmb->ncells2, nc3 = pmb->ncells3;

  xface_.resize(nc3 + nc2 + nc1 + 3);
  for (int k = 0; k <= nc3; ++k) xface_[k] = pmb->pcoord->x3f(k);
  for (int j = 0; j <= nc2; ++j) xface_[nc3 + 1 + j] = pmb->pcoord->x2f(j);
  for (int i = 0; i <= nc1; ++i)
    xface_[nc3 + nc2 + 2 + i] = pmb->pcoord->x1f(i);

  xcenter_.resize(nc3 + nc2 + nc1);
  for (int k = 0; k < nc3; ++k) xcenter_[k] = pmb->pcoord->x3v(k);
  for (int j = 0; j < nc2; ++j) xcenter_[nc3 + j] = pmb->pcoord->x2v(j);
  for (int i = 0; i < nc1; ++i) xcenter_[nc3 + nc2 + i] = pmb->pcoord->x1v(i);

  dims_.resize(3);
  dims_[0] = nc3;
  dims_[1] = nc2;
  dims_[2] = nc1;

  u.NewAthenaArray(nct, nc3, nc2, nc1);
  u.ZeroClear();
  u1_.NewAthenaArray(nct, nc3, nc2, nc1);
  u1_.ZeroClear();
  pcell_.NewAthenaArray(nct, nc3, nc2, nc1);

  seeds_per_cell_ =
      pin->GetOrAddInteger("particles", name + ".seeds_per_cell", 1);
  nmax_per_cell_ = pin->GetOrAddInteger("particles", name + ".nmax_per_cell",
                                        5 * seeds_per_cell_);
  density_floor_ = pin->GetOrAddReal("particles", name + ".dfloor", 1.E-10);

  if (nmax_per_cell_ < seeds_per_cell_) {
    msg << "### FATAL ERROR in Particles::Particles"
        << "Maximum particles per cell: " << nmax_per_cell_ << " is less than "
        << "seed particles per cell: " << seeds_per_cell_ << std::endl;
    ATHENA_ERROR(msg);
  }
  pmb->pdebug->Leave();
}

// destructor
Particles::~Particles() {
  if (prev != nullptr) prev->next = next;
  if (next != nullptr) next->prev = prev;
  delete ppb;
}

Particles::Particles(Particles const &other)
    : u(other.u), u1_(other.u1_), pcell_(other.pcell_) {
  if (this == &other) return;
  *this = other;
  ppb = new ParticleBuffer(this);
}

// functions
Particles *Particles::FindParticle(std::string name) {
  std::stringstream msg;
  Particles *p = this;

  while ((p != nullptr) && (p->myname != name)) p = p->next;
  if (p == nullptr) {
    msg << "### FATAL ERROR in Particles::FindParticles"
        << "Particles " << name << " not found" << std::endl;
    ATHENA_ERROR(msg);
  }

  return p;
}

void Particles::Initialize() {
  Particles *p = this;
  while (p != nullptr) {
    p->Particulate(p->mp, p->u);
    p = p->next;
  }
}

void Particles::TimeIntegrate(std::vector<MaterialPoint> &mp, Real time,
                              Real dt) {
  pmy_block->pdebug->Call("Particles::TimeIntegrate-" + myname);
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    for (std::vector<MaterialPoint>::iterator it = mp.begin(); it != mp.end();
         ++it) {
      it->x1 += it->v1 * dt;
      it->x2 += it->v2 * dt;
      it->x3 += it->v3 * dt;
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    for (std::vector<MaterialPoint>::iterator it = mp.begin(); it != mp.end();
         ++it) {
      it->x1 += it->v1 * dt;
      it->x2 += it->v2 * dt / it->x1;
      it->x3 += it->v3 * dt / (it->x1 * sin(it->x2));
    }
  }
  pmy_block->pdebug->Leave();
}

void Particles::WeightedAverage(std::vector<MaterialPoint> &mp_out,
                                std::vector<MaterialPoint> const &mp_in,
                                Real ave_wghts[]) {
  pmy_block->pdebug->Call("Particles::WeightedAverage-" + myname);
  size_t psize = mp_out.size();
  for (size_t i = 0; i < psize; ++i) {
    mp_out[i].x1 = ave_wghts[0] * mp_out[i].x1 + ave_wghts[1] * mp_in[i].x1;
    mp_out[i].x2 = ave_wghts[0] * mp_out[i].x2 + ave_wghts[1] * mp_in[i].x2;
    mp_out[i].x3 = ave_wghts[0] * mp_out[i].x3 + ave_wghts[1] * mp_in[i].x3;

    mp_out[i].v1 = ave_wghts[0] * mp_out[i].v1 + ave_wghts[1] * mp_in[i].v1;
    mp_out[i].v2 = ave_wghts[0] * mp_out[i].v2 + ave_wghts[1] * mp_in[i].v2;
    mp_out[i].v3 = ave_wghts[0] * mp_out[i].v3 + ave_wghts[1] * mp_in[i].v3;
  }
  pmy_block->pdebug->Leave();
}

size_t Particles::RestartDataSizeInBytes() {
  size_t size = 0;
  Particles *p = this;

  while (p != nullptr) {
    size += sizeof(int) + p->mp.size() * sizeof(MaterialPoint);
    p = p->next;
  }

  // gather maximum size
#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE, &size, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif

  return size;
}

size_t Particles::DumpRestartData(char *pdst) {
  Particles *p = this;

  while (p != nullptr) {
    int size = p->mp.size();
    std::memcpy(pdst, &size, sizeof(int));
    pdst += sizeof(int);
    std::memcpy(pdst, p->mp.data(), size * sizeof(MaterialPoint));
    pdst += size * sizeof(MaterialPoint);
    p = p->next;
  }
  return RestartDataSizeInBytes();
}

size_t Particles::LoadRestartData(char *psrc) {
  Particles *p = this;
  int size;

  while (p != nullptr) {
    std::memcpy(&size, psrc, sizeof(int));
    psrc += sizeof(int);
    p->mp.resize(size);
    std::memcpy(p->mp.data(), psrc, size * sizeof(MaterialPoint));
    psrc += size * sizeof(MaterialPoint);
    p = p->next;
  }

  return RestartDataSizeInBytes();
}
