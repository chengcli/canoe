// cnaoe
#include <configure.hpp>

// athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/reconstruct/interpolation.hpp>
#include <athena/scalars/scalars.hpp>
#include <athena/stride_iterator.hpp>

// canoe
#include <air_parcel.hpp>
#include <constants.hpp>
#include <impl.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// harp
#include "absorber.hpp"
#include "radiation.hpp"
#include "radiation_band.hpp"
// #include "../particles/particles.hpp"

// setting optical properties
void RadiationBand::SetSpectralProperties(int k, int j, int il, int iu) {
  int nspec = spec_.size();
  int npmom = bpmom.GetDim4() - 1;

  // set tau, ssalb, pmom, etc...
  tau_.ZeroClear();
  ssa_.ZeroClear();
  pmom_.ZeroClear();

  std::vector<Real> mypmom(1 + npmom);

  AirParcel var(AirParcel::Type::MoleFrac);
  // Particles *ppart;

  MeshBlock* pmb = pmy_block_;
  AthenaArray<Real> const& w = pmb->phydro->w;

  for (int i = il; i <= iu; ++i) {
    pmb->pimpl->GatherFromPrimitive(&var, k, j, i);
    tem_(i) = var.w[IDN];

    for (auto& a : absorbers_) {
      for (int m = 0; m < nspec; ++m) {
        Real kcoeff =
            a->GetAttenuation(spec_[m].wav1, spec_[m].wav2, var);  // 1/m
        Real dssalb =
            a->GetSingleScatteringAlbedo(spec_[m].wav1, spec_[m].wav2, var) *
            kcoeff;
        // tau
        tau_(m, i) += kcoeff;
        // ssalb
        ssa_(m, i) += dssalb;
        // pmom
        a->GetPhaseMomentum(mypmom.data(), spec_[m].wav1, spec_[m].wav2, var,
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
#ifdef HYDROSTATIC
      Real grav = -pmb->phydro->hsrc.GetG1();
      Real H0 = pmb->pimpl->GetPressureScaleHeight();
      // TODO(cli) check this
      // \delta z = \delt Z * P/(\rho g H)
      tau_(m, i) *= pmb->pcoord->dx1f(i) * w(IPR, k, j, i) /
                    (w(IDN, k, j, i) * grav * H0);
#else
      tau_(m, i) *= pmb->pcoord->dx1f(i);
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
