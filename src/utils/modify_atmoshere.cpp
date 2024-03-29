// C/C++
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <memory>
#include <vector>

// athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/eos/eos.hpp>
#include <athena/field/field.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>

// canoe
#include <air_parcel.hpp>
#include <constants.hpp>
#include <impl.hpp>
#include <index_map.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// helper functions, will be moved in the future
int find_pressure_level_lesser_pybind(Real pres, AthenaArray<Real> const &w,
                                      int k, int j, int is, int ie) {
  for (int i = is; i <= ie; ++i)
    if (w(IPR, k, j, i) < pres) return i;

  return ie + 1;
}

// modify atmoshere with adlnTdlnP
void modify_atmoshere_adlnTdlnP(MeshBlock *pmb, Real adlnTdlnP, Real pmin,
                                Real pmax) {
  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  Hydro *phydro = pmb->phydro;
  auto pthermo = Thermodynamics::GetInstance();
  auto pcoord = pmb->pcoord;
  Real H0 = pcoord->GetPressureScaleHeight();
  Real dlnp = pcoord->dx1f(is) / H0;

  // loop over all aircolumns
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      int ibegin =
          find_pressure_level_lesser_pybind(pmax, phydro->w, k, j, is, ie);
      int iend =
          find_pressure_level_lesser_pybind(pmin, phydro->w, k, j, is, ie);

      auto &&air = AirParcelHelper::gather_from_primitive(pmb, k, j, ibegin);
      air.ToMoleFraction();

      for (int i = ibegin; i < iend; ++i) {
        pthermo->Extrapolate(&air, -dlnp, "dry", 0., adlnTdlnP);
        AirParcelHelper::distribute_to_primitive(pmb, k, j, i + 1, air);
      }
    }
  }
};

// modify atmoshere with adlnNH3dlnP
void modify_atmoshere_adlnNH3dlnP(MeshBlock *pmb, Real adlnNH3dlnP, Real pmin,
                                  Real pmax) {
  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  Hydro *phydro = pmb->phydro;
  auto pthermo = Thermodynamics::GetInstance();
  auto pcoord = pmb->pcoord;
  Real H0 = pcoord->GetPressureScaleHeight();
  Real dlnp = pcoord->dx1f(is) / H0;

  // index
  auto pindex = IndexMap::GetInstance();
  int iNH3 = pindex->GetVaporId("NH3");

  // loop over all aircolumns
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      int ibegin =
          find_pressure_level_lesser_pybind(pmax, phydro->w, k, j, is, ie);
      int iend =
          find_pressure_level_lesser_pybind(pmin, phydro->w, k, j, is, ie);

      auto &&air = AirParcelHelper::gather_from_primitive(pmb, k, j, ibegin);
      air.ToMoleFraction();

      for (int i = ibegin; i < iend; ++i) {
        pthermo->Extrapolate(&air, -dlnp, "dry");
        air.w[iNH3] += adlnNH3dlnP * air.w[iNH3] * dlnp;
        auto rates = pthermo->TryEquilibriumTP_VaporCloud(air, iNH3);
        air.w[iNH3] += rates[0];
        AirParcelHelper::distribute_to_primitive(pmb, k, j, i + 1, air);
      }
    }
};
