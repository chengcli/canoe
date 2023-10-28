/** @file particulate.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Sunday Jun 13, 2021 12:16:32 PDT
 * @bug No known bugs.
 */

// C++ headers
#include <cassert>
#include <functional>  // hash
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include "../coordinates/coordinates.hpp"
#include "../debugger/debugger.hpp"
#include "../globals.hpp"
#include "particles.hpp"

// Particulate executes after chemistry, which is applied to aggreated
// quantities u It synchronized the particles in mp and the aggreated quantities
// after chemistry
void Particles::Particulate(std::vector<MaterialPoint> &mp,
                            AthenaArray<Real> const &u) {
  MeshBlock *pmb = pmy_block;
  pmb->pdebug->Call("Particles::Particulate-" + myname);

  Coordinates *pco = pmb->pcoord;
  // particle buffer
  std::vector<MaterialPoint> mpb;

  for (int t = 0; t < u.GetDim4(); ++t)
    for (int k = pmb->ks; k <= pmb->ke; ++k)
      for (int j = pmb->js; j <= pmb->je; ++j)
        for (int i = pmb->is; i <= pmb->ie; ++i) {
          // u1_ stores the total particle density before chemistry
          // u  stores the total particle density after chemistry
          // if delta_u is positive, new particles are added. Otherwise,
          // particle densities are adjusted to reflect the change of density
          Real delta_u = u(t, k, j, i) - u1_(t, k, j, i);
          // 1. count current number of particles in the cell
          int nparts = 0;
          MaterialPoint *pc = pcell_(t, k, j, i);
          while (pc != nullptr) {
            pc = pc->next;
            nparts++;
          }

          // 2. if the number of particles (n0) is greater than the maximum
          // allowed particles per cell (nmax_per_cell_), mark particles for
          // merging (destroy). Particle densities have been sorted from low to
          // high. particles to be removed (merged) are marked as inactive (id =
          // -1) and zero density (rho = 0). nparts reflects the number of
          // active particles in cell
          int n0 = nparts;
          Real drho = 0., sum = 0.;
          pc = pcell_(t, k, j, i);
          /* method 1
          if (pc != nullptr) {
            while (nparts > nmax_per_cell_ || pc->rho+drho < density_floor_) {
              if (nparts == 1) break;
              drho += (pc->rho+drho)/(nparts-1);
              pc->id = -1;
              pc->rho = 0;
              pc = pc->next;
              nparts--;
            }
          }*/

          // method 2
          for (int n = 0; n < n0 - nmax_per_cell_; ++n) {
            sum += pc->rho;
            pc->id = -1;
            pc->rho = 0;
            pc = pc->next;
            nparts--;
          }
          drho = sum / nparts;

          // Distribute the total density of inactive particles to active
          // particles.
          pcell_(t, k, j, i) = pc;
          while (pc != nullptr) {
            pc->rho += drho;
            pc = pc->next;
          }

          // 3. add new particles to the empty particle buffer mpb or increase
          // the density of existing particles
          if (delta_u > 0.) {
            // avg stores the mean density to be allocated to particles
            Real avg = delta_u / seeds_per_cell_;
            // num sotres the available number of particles to be added to cell
            int num = std::min(nmax_per_cell_ - nparts, seeds_per_cell_);
            // add new particles with density avg
            std::string str = myname + std::to_string(t) +
                              std::to_string(pmb->gid) + std::to_string(i) +
                              std::to_string(j) + std::to_string(k) + 't' +
                              std::to_string(pmb->pmy_mesh->time) + 'n';
            for (int n = 0; n < num; ++n) {
              MaterialPoint p;
              p.id = std::abs((int)std::hash<std::string>{}(
                         str + std::to_string(n))) %
                     (1 << 20);
              p.type = t;
              p.time = pmb->pmy_mesh->time;
              p.x1 = pco->x1f(i) + (1. * rand() / RAND_MAX) * pco->dx1f(i);
              p.x2 = pco->x2f(j) + (1. * rand() / RAND_MAX) * pco->dx2f(j);
              p.x3 = pco->x3f(k) + (1. * rand() / RAND_MAX) * pco->dx3f(k);
              p.v1 = 0.;
              p.v2 = 0.;
              p.v3 = 0.;
              p.rho = avg;
              mpb.push_back(p);
            }

            // add the remaining mass to existing particles
            pc = pcell_(t, k, j, i);
            for (int n = 0; n < seeds_per_cell_ - num; ++n) {
              pc->rho += avg;
              pc = pc->next;
            }
            // 4. removes existing particles
          } else if (delta_u < 0.) {
            Real avg = std::abs(delta_u) / nparts;
            pc = pcell_(t, k, j, i);
            // std::cout << "c =  " << u(t,k,j,i) << " c1 = " << u1_(t,k,j,i) <<
            // std::endl; std::cout << "[" << pco->x1f(i) << "," <<
            // pco->x1f(i+1) << "]" << std::endl; std::cout << "delta_u = " <<
            // delta_u << " nparts = " << nparts << std::endl;
            for (int n = 0; n < nparts; ++n) {
              // std::cout << pc->x1 << " " << pc->rho << " ";
              if (pc->rho > avg) {
                delta_u += avg;
                pc->rho -= avg;
              } else {
                delta_u += pc->rho;
                avg = std::abs(delta_u) / (nparts - n - 1);
                // available_ids_.push_back(pc->id);
                pc->id = -1;
                pc->rho = 0;
              }
              // std::cout << pc->rho << std::endl;
              // std::cout << delta_u << " " << avg << std::endl;
              pc = pc->next;
            }
            std::stringstream msg;
            /*if (std::abs(delta_u) > density_floor_) {
              msg << "### FATAL ERROR in Particles::Particulate:" << std::endl
                  << "(" << k << "," << j << "," << i << ") = " << nparts <<
            std::endl
                  << "u = " << u(t,k,j,i) << " u1 = " << u1_(t,k,j,i) <<
            std::endl
                  << "delta_u = " << delta_u << std::endl
                  << "density_floor_ = " << density_floor_ << std::endl;
              ATHENA_ERROR(msg);
            }*/
          }
        }

  // transfer from particle buffer (mpb) to storage (mp);
  // mp.reserve(mp.size() + mpb.size());
  // mp.insert(mp.end(), mpb.begin(), mpb.end());
  for (std::vector<MaterialPoint>::iterator it = mpb.begin(); it != mpb.end();
       ++it)
    mp.push_back(*it);

#if DEBUG_LEVEL > 2
  pmb->pdebug->CheckParticleConservation(cnames_, mp);
#endif

  pmb->pdebug->Leave();
}
