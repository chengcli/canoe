// athena
#include <athena/athena.hpp>
#include <athena/athena_arrays.hpp>
#include <athena/bvals/bvals.hpp>
#include <athena/coordinates/coordinates.hpp>
#include <athena/eos/eos.hpp>
#include <athena/field/field.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// canoe
#include <impl.hpp>
#include <index_map.hpp>

// climath
#include <climath/core.h>
#include <climath/interpolation.h>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// harp
#include <harp/radiation.hpp>

// special includes
#include <special/giants_enroll_vapor_functions_v1.hpp>

// utils
#include <utils/fileio.hpp>  // read_data_vector

Real grav, P0, T0;
int iSiO;

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(4 + NVAPOR);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
  SetUserOutputVariableName(2, "thetav");
  SetUserOutputVariableName(3, "mse");
  for (int n = 1; n <= NVAPOR; ++n) {
    std::string name = "rh" + std::to_string(n);
    SetUserOutputVariableName(3 + n, name.c_str());
  }
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  auto pthermo = Thermodynamics::GetInstance();

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        user_out_var(0, k, j, i) = pthermo->GetTemp(this, k, j, i);
        user_out_var(1, k, j, i) = pthermo->PotentialTemp(this, P0, k, j, i);
        // theta_v (potential virtual temperature)
        user_out_var(2, k, j, i) =
            user_out_var(1, k, j, i) * pthermo->RovRd(this, k, j, i);    // RovRd: gas constant of humid air over gas constant of dry air
        // mse
        user_out_var(3, k, j, i) =
            pthermo->MoistStaticEnergy(this, grav * pcoord->x1v(i), k, j, i);
        for (int n = 1; n <= NVAPOR; ++n)
          user_out_var(3 + n, k, j, i) =
              pthermo->RelativeHumidity(this, n, k, j, i);
      }
}

void Forcing(MeshBlock *pmb, Real const time, Real const dt,
             AthenaArray<Real> const &w, AthenaArray<Real> const &r,
             AthenaArray<Real> const &bcc, AthenaArray<Real> &du,
             AthenaArray<Real> &s) {
  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;
  auto pthermo = Thermodynamics::GetInstance();

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        auto &&air = AirParcelHelper::gather_from_primitive(pmb, k, j, i);

        Real cv = pthermo->GetCvMass(air, 0);

        // if (w(IPR, k, j, i) < prad) {
        //   du(IEN, k, j, i) += dt * hrate * w(IDN, k, j, i) * cv *
        //                       (1. + 1.E-4 * sin(2. * M_PI * rand() /
        //                       RAND_MAX));
        // }

        // if (air.w[IDN] < Tmin) {
        //   u(IEN,k,j,i) += w(IDN,k,j,i)*cv*(Tmin - temp)/sponge_tau*dt;
        // }
      }
}

void Mesh::InitUserMeshData(ParameterInput *pin) {
  grav = -pin->GetReal("hydro", "grav_acc1");

  P0 = pin->GetReal("problem", "P0");
  //T0 = pin->GetReal("problem", "T0");

  // index
  auto pindex = IndexMap::GetInstance();
  iSiO = pindex->GetVaporId("SiO");
  EnrollUserExplicitSourceFunction(Forcing);
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  srand(Globals::my_rank + time(0));

  Application::Logger app("main");
  app->Log("ProblemGenerator: lava_cir_rt");

  auto pthermo = Thermodynamics::GetInstance();

  // mesh limits
  Real x1min = pmy_mesh->mesh_size.x1min;
  Real x1max = pmy_mesh->mesh_size.x1max;

  // request temperature and pressure
  //app->Log("request T", T0);
  app->Log("request P", P0);

  // thermodynamic constants
  Real gamma = pin->GetReal("hydro", "gamma");
  Real Rd = pthermo->GetRd();
  Real cp = gamma / (gamma - 1.) * Rd;
  
  // read initial atmosphere conditions
  std::string input_atm_path = "/home/linfel/canoe_1/examples/2023-Linf-lava/balance_atm.txt";
  DataVector atm = read_data_vector(input_atm_path);
  size_t atm_size = atm["HGT"].size();

  AirParcel air(AirParcel::Type::MoleFrac);
  
  // convert km -> m, hPa -> Pa, ppm -> fraction
  double *hgtptr = atm["HGT"].data(); 
  double *preptr = atm["PRE"].data(); 
  double *sioptr = atm["SiO"].data(); 
  for (int i = 0; i < atm_size; ++i) {
    hgtptr[i] *= 1000.0;
    preptr[i] *= 100.0;
    sioptr[i] /= 1.E6;
  }


  // construct atmosphere from bottom up
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      
      // half a grid to cell center
      pthermo->Extrapolate(&air, pcoord->dx1f(is) / 2.,
                           Thermodynamics::Method::ReversibleAdiabat, grav);
      
      for (int i = is; i <= ie; ++i) {
        // interpolation to coord grid
        air.SetZero();
	air.w[IPR] = interp1(pcoord->x1v(i), atm["PRE"].data(),
                             atm["HGT"].data(), atm_size);
        air.w[IDN] = interp1(pcoord->x1v(i), atm["TEM"].data(),
                             atm["HGT"].data(), atm_size);
        air.w[iSiO] = interp1(pcoord->x1v(i), atm["SiO"].data(),
                              atm["HGT"].data(), atm_size);
        //air.w[IVX] = 1. * sin(2. * M_PI * rand() / RAND_MAX);
  	
	
	//std::cout << " i  " << i << " j " << j << " k " << k << "  coord  " << pcoord->x1v(i) << std::endl;
        //std::cout << "  coord"  << pcoord->x1v(i) << " IPR  " << air.w[IPR] << " IDN  " <<  air.w[IDN] << "  isio  " << air.w[iSiO] << std::endl;
	
	AirParcelHelper::distribute_to_conserved(this, k, j, i, air);
        pthermo->Extrapolate(&air, pcoord->dx1f(i),
                             Thermodynamics::Method::PseudoAdiabat, grav,
                             1.e-5);
      }

      peos->ConservedToPrimitive(phydro->u, phydro->w, pfield->b, phydro->w,
                                 pfield->bcc, pcoord, is, ie, j, j, k, k);

      pimpl->prad->CalFlux(this, k, j, is, ie + 1);
      //std::cout << " =======================================================" << std::endl;
      
    }
}

