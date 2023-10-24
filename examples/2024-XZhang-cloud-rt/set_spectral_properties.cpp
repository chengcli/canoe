/** @file set_spectral_properties.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Friday Apr 15, 2022 10:41:21 EDT
 * @bug No known bugs.
 */

// Athena++ header
#include "radiation.hpp"
#include "absorber.hpp"
#include "../mesh/mesh.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "../coordinates/coordinates.hpp"
#include "../particles/particles.hpp"
#include "../hydro/hydro.hpp"
#include "../scalars/scalars.hpp"
#include "../reconstruct/interpolation.hpp"

// setting optical properties
void RadiationBand::setSpectralProperties(AthenaArray<Real> const& w,
  int k, int j, int il, int iu)
{
  MeshBlock *pmb = pmy_rad->pmy_block;
  int is = pmb->is, ie = pmb->ie;
  int nspec = spec.size();
  int npmom = bpmom.GetDim4() - 1;

  // set tau, ssalb, pmom, etc...
  int ncells1 = pmy_rad->pmy_block->ncells1;
  std::fill(*tau_, *tau_ + nspec*ncells1, 0.);
  std::fill(*ssa_, *ssa_ + nspec*ncells1, 0.);
  std::fill(**pmom_, **pmom_ + nspec*ncells1*(npmom+1), 0.);

  Absorber *a = pabs;
  Thermodynamics *pthermo = pmy_rad->pmy_block->pthermo;
  Coordinates *pcoord = pmy_rad->pmy_block->pcoord;
  Hydro *phydro = pmy_rad->pmy_block->phydro;
  PassiveScalars *pscalars = pmy_rad->pmy_block->pscalars;

  std::vector<Real> mypmom(1+npmom);
  CellVariables var(pmb->GetNumVariablesInCell());
  Particles *ppart;

  while (a != nullptr) {
    for (int i = il; i <= iu; ++i) {
      for (int n = 0; n < NSCALARS; ++n) var.s[n] = pscalars->s(n,k,j,i);
      pthermo->PrimitiveToChemical(var.q, w.at(k,j,i));
      //! \todo do we need it?
      // molar concentration to molar mixing ratio
      Real nmols = var.q[IPR]/(Thermodynamics::Rgas*var.q[IDN]);
      for (int n = 1; n <= NVAPOR; ++n) var.q[n] /= nmols;

      // molar density of clouds, mol/m^3, xiz changed to kg/m^3
      ppart = pmb->ppart;
      int ip = 0;
      while (ppart != nullptr) {
        for (int n = 0; n < ppart->u.GetDim4(); ++n)
          //var.c[ip++] = ppart->u(n,k,j,i)/ppart->GetMolecularWeight(n);
          var.c[ip++] = ppart->u(n,k,j,i);
        ppart = ppart->next;
      }

      tem_[i] = var.q[IDN];
      //std::cout << i << " " << tem_[i] << std::endl;
      for (int m = 0; m < nspec; ++m) {
        Real kcoeff = a->getAttenuation(spec[m].wav1, spec[m].wav2, var);  // 1/m
        Real dssalb = a->getSingleScatteringAlbedo(spec[m].wav1, spec[m].wav2, var)*kcoeff;

        // tau 
        tau_[m][i] += kcoeff;
        // ssalb
        ssa_[m][i] += dssalb;

        // pmom
        a->getPhaseMomentum(mypmom.data(), spec[m].wav1, spec[m].wav2, var, npmom);
        for (int p = 0; p <= npmom; ++p)
          pmom_[m][i][p] += mypmom[p]*dssalb;
      }
    }
    a = a->next;
  }

  /* set temperature at cell interface*/

  temf_[il] = 3.*tem_[il] - 2.*tem_[il+1];
  temf_[il+1] = (tem_[il] + tem_[il+1])/2.;
  for (int i = il+2; i <= iu-1; ++i)
    temf_[i] = interp_cp4(tem_[i-2], tem_[i-1], tem_[i], tem_[i+1]);
  temf_[iu] = (tem_[iu] + tem_[iu-1])/2.;
  temf_[iu+1] = 3.*tem_[iu] - 2.*tem_[iu-1];

/* 
  Real r = pmb->pmy_mesh->mesh_size.x1rat;
  temf_[il] = (r*(2 + r)*tem_[il] - tem_[il+1])/(1 + r);
  temf_[il+1] = (r*r*tem_[il] + tem_[il+1])/(1 + r);
  for (int i = il+2; i <= iu-1; ++i)
    temf_[i] = (-pow(r,8)*tem_[i-2] + pow(r,4)*(1 + 2*r*(1 + r + r*r))*tem_[i-1] + 
                r*(2 + r*(2 + r*(2 + r)))*tem_[i] - tem_[i+1])
                /(r*(1 + r)*(1 + r*r)*(1 + r + r*r));
  temf_[iu+1] = (r*r*tem_[iu] + tem_[iu-1])/(1 + r);
  temf_[iu] = (r*(2 + r)*tem_[iu] - tem_[iu-1])/(1 + r);
*/

  // absorption coefficiunts -> optical thickness
  for (int m = 0; m < nspec; ++m) {
    for (int i = il; i <= iu; ++i) {
      if (tau_[m][i] > 1e-6 && ssa_[m][i] > 1e-6) {  // has scattering
        for (int p = 0; p <= npmom; ++p)
          pmom_[m][i][p] /= ssa_[m][i];
        ssa_[m][i] /= tau_[m][i];
      } else {
        ssa_[m][i] = 0.;
        pmom_[m][i][0] = 1.;
        for (int p = 1; p <= npmom; ++p)
          pmom_[m][i][p] = 0.;
      }
      if (HYDROSTATIC) {
        Real grav = -phydro->hsrc.GetG1();
        // \delta z = \delt Z * P/(\rho g H), cli, TODO
        tau_[m][i] *= pcoord->dx1f(i)*w(IPR,k,j,i)/(w(IDN,k,j,i)*grav*phydro->scale_height);
      } else {
        tau_[m][i] *= pcoord->dx1f(i);
      }
    }
  }

  // aggregated band properties
  for (int i = il; i <= iu; ++i) {
    btau(k,j,i) = 0; bssa(k,j,i) = 0;
    for (int p = 0; p <= npmom; ++p) bpmom(p,k,j,i) = 0.;

    for (int m = 0; m < nspec; ++m) {
      btau(k,j,i) += tau_[m][i];
      bssa(k,j,i) += ssa_[m][i]*tau_[m][i];
      for (int p = 0; p <= npmom; ++p)
        bpmom(p,k,j,i) += pmom_[m][i][p]*ssa_[m][i];
    }

    for (int p = 0; p <= npmom; ++p)
      bpmom(p,k,j,i) /= bssa(k,j,i);
    bssa(k,j,i) /= btau(k,j,i);
    btau(k,j,i) /= nspec;
  }

//xiz 2022 tried the total optical depth output		
  for (int i = iu-1; i >= il; --i) {
    btau(k,j,i) += btau(k,j,i+1);
  }
}


