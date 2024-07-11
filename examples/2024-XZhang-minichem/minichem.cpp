// C++ headers
#include <cmath>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>

// athena
#include <athena/eos/eos.hpp>
#include <athena/field/field.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>
#include <athena/stride_iterator.hpp>
#include <athena/scalars/scalars.hpp>

// canoe
#include <air_parcel.hpp>
#include <configure.hpp>
#include <impl.hpp>
#include <athena/coordinates/coordinates.hpp>

// climath
#include <climath/core.h>

// exo3
#include <exo3/cubed_sphere.hpp>
#include <exo3/cubed_sphere_utility.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// astro
#include <astro/celestrial_body.hpp>

// harp
#include <harp/radiation.hpp>

// minichem
#include <minichem/mini_chem.hpp>

Real Ps, Ts, Omega, grav, sponge_tau, bflux, gammad, Tmin;
int sponge_layer;
// species ['OH','H2','H2O','H','CO','CO2','O','CH4','C2H2','NH3','N2','HCN']
std::vector<double> vmass = {17.01, 2.02, 18.02, 1.01, 28.01, 44.01, 16., 16.05, 26.04, 17.04, 28.02, 27.03};
Real mmass = 2.238; //mean molecular mass in amu

MiniChem *mc;
/*
void vmr_from_prim_scalar(std::vector<Real> &vmr, 
                          AthenaArray<Real> const& prim_scalar, 
                          AthenaArray<Real> const& w,
                          int k, int j, int i)
{
  if (NCHEMISTRY > 0) {
     for (int n=0; n<NCHEMISTRY; ++n) {
        vmr[n] = prim_scalar(NCLOUD + n, k, j, i)/vmass[n]*mmass; //change mmr to vmr
      }
    }
}

void vmr_to_cons_scalar(AthenaArray<Real> &cons_scalar, 
                        std::vector<Real> const& vmr, 
                        AthenaArray<Real> const& w,
                        int k, int j, int i)
{
  if (NCHEMISTRY > 0) {
     for (int n=0; n<NCHEMISTRY; ++n) {
        cons_scalar(NCLOUD + n, k, j, i) = w(IDN, k, j, i)*vmr[n]*vmass[n]/mmass; //change vmr to den
      }
    }
}
*/
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(7);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
  SetUserOutputVariableName(2, "lat");
  SetUserOutputVariableName(3, "lon");
  SetUserOutputVariableName(4, "vlat");
  SetUserOutputVariableName(5, "vlon");
  SetUserOutputVariableName(6, "zenith");
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  auto pexo3 = pimpl->pexo3;
  auto pthermo = Thermodynamics::GetInstance();

  Real lat, lon;
  Real U, V;
  Direction ray = pimpl->prad->GetRayInput(0);
  Real zenith;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        user_out_var(0, k, j, i) = pthermo->GetTemp(this, k, j, i);
        user_out_var(1, k, j, i) = pthermo->PotentialTemp(this, Ps, k, j, i);

        pexo3->GetLatLon(&lat, &lon, k, j, i);
        pexo3->GetUV(&U, &V, phydro->w(IVY, k, j, i), phydro->w(IVZ, k, j, i),
                     k, j, i);
        user_out_var(2, k, j, i) = lat;
        user_out_var(3, k, j, i) = lon;
        user_out_var(4, k, j, i) = U;
        user_out_var(5, k, j, i) = V;

        ray = pimpl->planet->ParentZenithAngle(pmy_mesh->time, M_PI / 2. - lat,
                                               lon);
        zenith = std::acos(ray.mu) / M_PI * 180.0;
        user_out_var(6, k, j, i) = zenith;
      }
    }
}

void Forcing(MeshBlock *pmb, Real const time, Real const dt,
             AthenaArray<Real> const &w, const AthenaArray<Real> &prim_scalar,
             AthenaArray<Real> const &bcc, AthenaArray<Real> &du,
             AthenaArray<Real> &cons_scalar) {
  auto pexo3 = pmb->pimpl->pexo3;
  auto pthermo = Thermodynamics::GetInstance();
  auto prad = pmb->pimpl->prad;
  auto phydro = pmb->phydro;

  std::vector<Real> vmr(NCHEMISTRY);

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) 
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real lat, lon;
        pexo3->GetLatLon(&lat, &lon, k, j, i);

        // coriolis force
        Real f = 2. * Omega * sin(lat);
        Real f2 = 2. * Omega * cos(lat);
        Real U, V;

        pexo3->GetUV(&U, &V, w(IVY, k, j, i), w(IVZ, k, j, i), k, j, i);

        Real m1 = w(IDN, k, j, i) * w(IVX, k, j, i);
        Real m2 = w(IDN, k, j, i) * U;
        Real m3 = w(IDN, k, j, i) * V;

        Real ll_acc_U = f * m3;
        Real ll_acc_V = -f * m2;
        Real acc2, acc3;
        pexo3->GetVyVz(&acc2, &acc3, ll_acc_U, ll_acc_V, k, j, i);
        pexo3->ContravariantVectorToCovariant(j, k, acc2, acc3, &acc2, &acc3);
        du(IM2, k, j, i) += dt * acc2;
        du(IM3, k, j, i) += dt * acc3;

        // minichem
        Real temp = pthermo->GetTemp(w.at(k, j , i));
        Real pres = w(IPR, k, j, i);

        //change mmr to vmr
        for (int n=0; n<NCHEMISTRY; ++n) vmr[n] = prim_scalar(NCLOUD + n, k, j, i)/vmass[n]*mmass; 
//       vmr_from_prim_scalar(vmr, prim_scalar, w, k, j, i);
        // call minichem        
        mc->Run(temp, pres, dt, vmr.data(), "NCHO");

       //normalize scale VMR to 1 and change to density
        Real sumVMR = std::accumulate(vmr.begin(), vmr.end(), static_cast<Real>(0));
        for (int n=0; n<NCHEMISTRY; ++n) {
          cons_scalar(NCLOUD + n, k, j, i) = phydro->w(IDN, k, j, i)*vmr[n]/sumVMR*vmass[n]/mmass; 
        }
//        vmr_to_cons_scalar(cons_scalar, vmr, w, k, j, i);

/*        
        Real rho = w(IDN, k, j, i);
        if (i > pmb->is+10) {
          for (int n = 0; n < NCHEMISTRY; ++n) {
             cons_scalar(NCLOUD + n, k, j, i) -= rho * prim_scalar(NCLOUD + n, k, j, i) * n * dt / 1.E20;
          }
        } else {
          for (int n = 0; n < NCHEMISTRY; ++n) {
            cons_scalar(NCLOUD + n, k, j, i) = rho*1.E-6;
          }
        }
*/

      }
}

Real TimeStep(MeshBlock *pmb) {
  auto prad = pmb->pimpl->prad;
  Real time = pmb->pmy_mesh->time;

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) {
      prad->CalTimeStep(pmb, k, j, pmb->is, pmb->ie);
    }

  return prad->GetTimeStep() * (1. + time / prad->GetRelaxTime());
}

void Mesh::InitUserMeshData(ParameterInput *pin) {

  // forcing parameters
  Omega = pin->GetReal("problem", "Omega");
  grav = -pin->GetReal("hydro", "grav_acc1");
  Ts = pin->GetReal("problem", "Ts");
  Ps = pin->GetReal("problem", "Ps");
  gammad = pin->GetReal("hydro", "gamma");
  Tmin = pin->GetReal("problem", "Tmin");

  // forcing function
  EnrollUserExplicitSourceFunction(Forcing);
  // EnrollUserTimeStepFunction(TimeStep);
  
  // minichem
  mc = new MiniChem();
  mc->SetDataFile("chem_data/mini_chem_data_NCHO.txt");
  mc->SetSpeciesFile("chem_data/mini_chem_sp_NCHO.txt");
  mc->SetNetworkDir("chem_data/1x/");
  mc->SetMetallicityStr("1x");
  mc->Initialize();

}

//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  srand(Globals::my_rank + time(0));
  auto pexo3 = pimpl->pexo3;
  auto pthermo = Thermodynamics::GetInstance();

  AirParcel air(AirParcel::Type::MoleFrac);
  srand(Globals::my_rank + time(0));

  // construct atmosphere from bottom up
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      air.SetZero();
      air.w[IPR] = Ps;
      air.w[IDN] = Ts;

      // half a grid to cell center
      pthermo->Extrapolate(&air, pcoord->dx1f(is) / 2., "reversible", grav);

      int i = is;
      for (; i <= ie; ++i) {
        air.w[IVX] = 0.001 * (1. * rand() / RAND_MAX - 0.5);
        if (air.w[IDN] < Tmin) break;
        AirParcelHelper::distribute_to_primitive(this, k, j, i, air);
        AirParcelHelper::distribute_to_conserved(this, k, j, i, air);
        pthermo->Extrapolate(&air, pcoord->dx1f(i), "pseudo", grav);
      }

      // Replace adiabatic atmosphere with isothermal atmosphere if temperature
      // is too low
      for (; i <= ie; ++i) {
        air.w[IVX] = 0.001 * (1. * rand() / RAND_MAX - 0.5);
        AirParcelHelper::distribute_to_primitive(this, k, j, i, air);
        AirParcelHelper::distribute_to_conserved(this, k, j, i, air);
        pthermo->Extrapolate(&air, pcoord->dx1f(i), "isothermal", grav);
      }
    }

  // initialize conserved variables
  Real Rd = pthermo->GetRd();
  int n_sp = 13;
  // interpolate from ce table
  std::vector<double> vmr_ic(n_sp);
  std::string ic_file = "chem_data/IC/mini_chem_IC_FastChem_1x.txt";

  for (int k=ks; k<=ke; k++) 
    for (int j=js; j<=je; j++) 
      for (int i=is; i<=ie; i++) {
         double mu;
         double T_in = pthermo->GetTemp(this, k, j, i);
         double P_in = phydro->w(IPR, k, j, i);
//         std::cout<<"interp  "<<P_in<<" "<<T_in<<std::endl;
         interp_ce_table(n_sp, T_in, P_in, vmr_ic.data(), &mu, ic_file);

        //normalize scale VMR to 1
        Real sumVMR = std::accumulate(vmr_ic.begin(), vmr_ic.end(), static_cast<Real>(0));
        for (auto& value : vmr_ic) {value /= sumVMR;}
   
         if (NCHEMISTRY > 0) {
           for (int n=0; n<NCHEMISTRY; ++n) {
//           std::cout<<"mini  "<<rho<<" "<<n<<" "<<vmr_ic[n]<<" "<<vmass[n]<<" "<<mu<<std::endl;
            pscalars->s(NCLOUD+n,k,j,i) = phydro->w(IDN, k, j, i)*vmr_ic[n]*vmass[n]/mmass;
           }
         }
      }

}
