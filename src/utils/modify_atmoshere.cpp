// C/C++
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <memory>
#include <vector>
#include <map>

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
void modify_adlnTdlnP(MeshBlock *pmb, Real adlnTdlnP, Real pmin,
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

// modify atmoshere with adlnXdlnP
void modify_adlnXdlnP(MeshBlock *pmb, Real adlnXdlnP, Real pmin,
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
        air.w[iNH3] += adlnXdlnP * air.w[iNH3] * dlnp;
        auto rates = pthermo->TryEquilibriumTP_VaporCloud(air, iNH3);
        air.w[iNH3] += rates[0];
        AirParcelHelper::distribute_to_primitive(pmb, k, j, i + 1, air);
      }
    }
};

void modify_atm(MeshBlock *pmb, 
                std::map<std::string, std::vector<Real>> const& atm)
{
  int is = pmb->is, ie = pmb->ie;

  if (atm.find("TEM") == atm.end() || atm.find("PRE") == atm.end()) {
    throw std::runtime_error("The atmosphere data is not complete.");
  }

  int nlayer = atm.at("TEM").size();

  if (nlayer != ie - is + 1) {
    throw std::runtime_error("The number of layers is not consistent.");
  }

  AirColumn ac(pmb->ncells1);

  // temperature and pressure
  for (int i = 0; i < nlayer; ++i) {
    ac[is + i].SetType(AirParcel::Type::MoleFrac);
    ac[is + i].w[IDN] = atm.at("TEM")[i];
    ac[is + i].w[IPR] = atm.at("PRE")[i];
  }

  auto pindex = IndexMap::GetInstance();

  for (auto const &pair : atm) {
    if (pair.first == "HGT" || pair.first == "TEM" ||
        pair.first == "PRE" || pair.first == "IDX") {
      continue;
    }

    if (pindex->HasVapor(pair.first))
      for (int i = 0; i < nlayer; ++i) {
        ac[is + i].w[pindex->GetVaporId(pair.first)] = pair.second[i];
      }

    if (pindex->HasCloud(pair.first))
      for (int i = 0; i < nlayer; ++i) {
        ac[is + i].c[pindex->GetCloudId(pair.first)] = pair.second[i];
      }

    if (pindex->HasChemistry(pair.first))
      for (int i = 0; i < nlayer; ++i) {
        ac[is + i].q[pindex->GetChemistryId(pair.first)] = pair.second[i];
      }

    if (pindex->HasTracer(pair.first))
      for (int i = 0; i < nlayer; ++i) {
        ac[is + i].x[pindex->GetTracerId(pair.first)] = pair.second[i];
      }
  }

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) {
      AirParcelHelper::distribute_to_primitive(pmb, k, j, ac);
    }
};
