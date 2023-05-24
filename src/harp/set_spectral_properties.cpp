/** @file set_spectral_properties.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Friday Apr 15, 2022 10:41:21 EDT
 * @bug No known bugs.
 */

// Athena++ header
#include <athena/coordinates/coordinates.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/scalars/scalars.hpp>
#include <athena/stride_iterator.hpp>

// harp
#include <configure.hpp>
#include <snap/thermodynamics/thermodynamics.hpp>

#include "absorber.hpp"
#include "radiation.hpp"
#include "radiation_band.hpp"
// #include "../particles/particles.hpp"
#include <snap/cell_variables.hpp>
#include <snap/reconstruct/interpolation.hpp>

// setting optical properties
void RadiationBand::setSpectralProperties(int k, int j, int il, int iu) {
  int nspec = spec_.size();
  int npmom = bpmom.GetDim4() - 1;

  // set tau, ssalb, pmom, etc...
  tau_.ZeroClear();
  ssa_.ZeroClear();
  pmom_.ZeroClear();

  std::vector<Real> mypmom(1 + npmom);

  CellVariables var;
  // Particles *ppart;

  AthenaArray<Real> const& w = phydro_->w;

  for (auto a : absorbers) {
    for (int i = il; i <= iu; ++i) {
      for (int n = 0; n < NSCALARS; ++n) var.s[n] = pscalars_->s(n, k, j, i);
      pthermo_->PrimitiveToChemical(var.w, w.at(k, j, i));
      //! \todo do we need it?
      // molar concentration to molar mixing ratio
      Real nmols = var.w[IPR] / (Thermodynamics::Rgas * var.w[IDN]);
      for (int n = 1; n <= NVAPOR; ++n) var.w[n] /= nmols;

      /* molar density of clouds, mol/m^3
      ppart = pmb->ppart;
      int ip = 0;
      while (ppart != nullptr) {
        for (int n = 0; n < ppart->u.GetDim4(); ++n)
          var.c[ip++] = ppart->u(n,k,j,i)/ppart->GetMolecularWeight(n);
        ppart = ppart->next;
      }*/

      tem_(i) = var.w[IDN];
      // std::cout << i << " " << tem_(i] << std::endl;
      for (int m = 0; m < nspec; ++m) {
        Real kcoeff =
            a->getAttenuation(spec_[m].wav1, spec_[m].wav2, var);  // 1/m
        Real dssalb =
            a->getSingleScatteringAlbedo(spec_[m].wav1, spec_[m].wav2, var) *
            kcoeff;
        // tau
        tau_(m, i) += kcoeff;
        // ssalb
        ssa_(m, i) += dssalb;
        // pmom
        a->getPhaseMomentum(mypmom.data(), spec_[m].wav1, spec_[m].wav2, var,
                            npmom);
        for (int p = 0; p <= npmom; ++p) pmom_(m, i, p) += mypmom[p] * dssalb;
      }
    }
  }

  // set temperature at cell interface
  temf_(il) = 3. * tem_(il) - 2. * tem_(il + 1);
  temf_(il + 1) = (tem_(il) + tem_(il + 1)) / 2.;
  for (int i = il + 2; i <= iu - 1; ++i)
    temf_(i) = interp_cp4(tem_(i - 2), tem_(i - 1), tem_(i), tem_(i + 1));
  temf_(iu) = (tem_(iu) + tem_(iu - 1)) / 2.;
  temf_(iu + 1) = 3. * tem_(iu) - 2. * tem_(iu - 1);

  // absorption coefficiunts -> optical thickness
  for (int m = 0; m < nspec; ++m) {
    for (int i = il; i <= iu; ++i) {
      if (tau_(m, i) > 1e-6 && ssa_(m, i) > 1e-6) {  // has scattering
        for (int p = 0; p <= npmom; ++p) pmom_(m, i, p) /= ssa_(m, i);
        ssa_(m, i) /= tau_(m, i);
      } else {
        ssa_(m, i) = 0.;
        pmom_(m, i, 0) = 1.;
        for (int p = 1; p <= npmom; ++p) pmom_(m, i, p) = 0.;
      }
#if HYDROSTATIC
      Real grav = -phydro_->hsrc.GetG1();
      // \delta z = \delt Z * P/(\rho g H), cli, TODO
      tau_(m, i) *= pcoord_->dx1f(i) * w(IPR, k, j, i) /
                    (w(IDN, k, j, i) * grav * Constants::PressureScaleHeight);
#else
      tau_(m, i) *= pcoord_->dx1f(i);
#endif
    }
  }

  // aggregated band properties
  for (int i = il; i <= iu; ++i) {
    btau(k, j, i) = 0;
    bssa(k, j, i) = 0;
    for (int p = 0; p <= npmom; ++p) bpmom(p, k, j, i) = 0.;

    for (int m = 0; m < nspec; ++m) {
      btau(k, j, i) += tau_(m, i);
      bssa(k, j, i) += ssa_(m, i) * tau_(m, i);
      for (int p = 0; p <= npmom; ++p)
        bpmom(p, k, j, i) += pmom_(m, i, p) * ssa_(m, i);
    }

    for (int p = 0; p <= npmom; ++p) bpmom(p, k, j, i) /= bssa(k, j, i);
    bssa(k, j, i) /= btau(k, j, i);
    btau(k, j, i) /= nspec;
  }
}
