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
#include <dirty.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// modify atmoshere with adlnTdlnP
void modify_atmoshere_adlnTdlnP(MeshBlock *pmb, Real adlnTdlnP, Real pmin, Real pmax) {

    Hydro* phydro=pmb->phydro;

    // loop over all aircolumns
    for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
            int ibegin = find_pressure_level_lesser(pmax, phydro->w, k, j, is, ie);
            int iend = find_pressure_level_lesser(pmin, phydro->w, k, j, is, ie);

            auto &&air = AirParcelHelper::gather_from_primitive(this, k, j, ibegin);
            air.ToMoleFraction();

            for (int i = ibegin; i < iend; ++i) {
            pthermo->Extrapolate(&air, -dlnp, Thermodynamics::Method::DryAdiabat,
                                0., adlnTdlnP);
            AirParcelHelper::distribute_to_primitive(this, k, j, i + 1, air);
            }
        }
    }
};

// modify atmoshere with adlnNH3dlnP
void modify_atmoshere_adlnNH3dlnP(MeshBlock *pmb, Real adlnNH3dlnP, Real pmin, Real pmax){

    Hydro* phydro=pmb->phydro;

    // loop over all aircolumns
    for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
            int ibegin = find_pressure_level_lesser(pmax, phydro->w, k, j, is, ie);
            int iend = find_pressure_level_lesser(pmin, phydro->w, k, j, is, ie);

            auto &&air = AirParcelHelper::gather_from_primitive(this, k, j, ibegin);
            air.ToMoleFraction();

            for (int i = ibegin; i < iend; ++i) {
            pthermo->Extrapolate(&air, -dlnp, Thermodynamics::Method::DryAdiabat,
                                0., adlnNH3dlnP);
            AirParcelHelper::distribute_to_primitive(this, k, j, i + 1, air);
            }
        }
    }
}




