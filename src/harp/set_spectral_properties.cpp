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

// opacity
#include <opacity/absorber.hpp>

// harp
#include "radiation.hpp"
#include "radiation_band.hpp"

// setting optical properties
void RadiationBand::SetSpectralProperties(AirColumn& ac, Real const* x1f,
                                          Real gH0, int k, int j) {
  int nspec = GetNumSpecGrids();
  int npmom = GetNumPhaseMoments();

  // set tau, ssalb, pmom, etc...
  tau_.ZeroClear();
  ssa_.ZeroClear();
  pmom_.ZeroClear();

  std::vector<Real> mypmom(1 + npmom);

  for (int i = 0; i < ac.size(); ++i) {
    auto& air = ac[i];
    air.ToMoleFraction();
    tem_[i] = air.w[IDN];

    for (auto& a : absorbers_) {
      for (int m = 0; m < nspec; ++m) {
        auto& spec = pgrid_->spec[m];
        Real kcoeff = a->GetAttenuation(spec.wav1, spec.wav2, air);  // 1/m
        Real dssalb =
            a->GetSingleScatteringAlbedo(spec.wav1, spec.wav2, air) * kcoeff;
        // tau
        tau_(m, i) += kcoeff;
        // ssalb
        ssa_(m, i) += dssalb;
        // pmom
        a->GetPhaseMomentum(mypmom.data(), spec.wav1, spec.wav2, air, npmom);
        for (int p = 0; p <= npmom; ++p) pmom_(m, i, p) += mypmom[p] * dssalb;
        std::cout << "absorber = " << a->GetName() << std::endl;
        std::cout << "tau[" << i << "] = " << kcoeff << std::endl;
        std::cout << "ssa[" << i << "] = " << dssalb << std::endl;
        std::cout << "pmom[" << i << "] = ";
        for (int p = 0; p <= npmom; ++p) std::cout << mypmom[p] << " ";
        std::cout << std::endl << std::endl;
      }
    }
  }

  // set temperature at cell interface
  int il = 0, iu = ac.size() - 1;
  temf_[il] = 3. * tem_[il] - 2. * tem_[il + 1];
  temf_[il + 1] = (tem_[il] + tem_[il + 1]) / 2.;
  for (int i = il + 2; i <= iu - 1; ++i)
    temf_[i] = interp_cp4(tem_[i - 2], tem_[i - 1], tem_[i], tem_[i + 1]);
  temf_[iu] = (tem_[iu] + tem_[iu - 1]) / 2.;
  temf_[iu + 1] = 3. * tem_[iu] - 2. * tem_[iu - 1];

  // absorption coefficients -> optical thickness
  for (int m = 0; m < nspec; ++m) {
    for (int i = 0; i < ac.size(); ++i) {
      if (tau_(m, i) > 1e-6 && ssa_(m, i) > 1e-6) {  // has scattering
        for (int p = 0; p <= npmom; ++p) pmom_(m, i, p) /= ssa_(m, i);
        ssa_(m, i) /= tau_(m, i);
      } else {
        ssa_(m, i) = 0.;
        pmom_(m, i, 0) = 1.;
        for (int p = 1; p <= npmom; ++p) pmom_(m, i, p) = 0.;
      }
#ifdef HYDROSTATIC
      auto pthermo = Thermodynamics::GetInstance();
      Real Rgas = pthermo->RovRd(ac[i]) * pthermo->GetRd();
      // TODO(cli) check this
      // \delta z = \delt Z * (R T)/(g H0)
      tau_(m, i) *= (x1f[i + 1] - x1f[i]) * (Rgas * tem_[i]) / (gH0);
#else
      tau_(m, i) *= x1f[i + 1] - x1f[i];
#endif
    }
  }

  // aggregated band properties
  for (int i = 0; i < ac.size(); ++i) {
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
