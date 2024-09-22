/* -------------------------------------------------------------------------------------
 * SNAP Example Program
 *
 * Contributer:
 * Cheng Li, University of Michigan
 *
 * Year: 2023
 * Contact: chengcli@umich.edu
 * Reference: Test Jupiter CRM
 * -------------------------------------------------------------------------------------
 */
// C++ headers
#include <cmath>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>

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
#include <athena/scalars/scalars.hpp>

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
#include <snap/stride_iterator.hpp>
#include <snap/thermodynamics/atm_thermodynamics.hpp>

// special includes
#include <special/giants_enroll_vapor_functions_v1.hpp>

// minichem
#include <minichem/mini_chem.hpp>

// global parameters
Real grav, P0, T0, Z0, Tmin, radius, omega, nu_scalar_iso, tau_scale_factor;
bool use_polar_beta;
int s_ind;

/*
The following are the "primed" coefficients.
CO example for low T --> [CO(gas) - {C(graphite) + 0.5*O2(gas)}].
Primed arrays are defined to limit calls to the Gibbs functionality:
a_prime_CO_low = a_CO_low - (a_C_gr_low + 0.5*a_O2_low)
*/
// NASA-9 polynomials to be treated as globals for CO (for deltaG(T)):
const std::vector<double> a_prime_CO_low = {-8.12674056e+04, 1.44584303e+03, -8.48882019e+00,
                                            3.60377848e-02, -8.72224514e-05, 9.84235413e-08,
                                            -4.22141668e-11, -2.02794511e+04, 5.58505087e+01};
const std::vector<double> a_prime_CO_hig = {6.45288795e+05, -5.20591633e+02, -1.94199375e+00,
                                            2.28448402e-03, -1.59490759e-06, 4.77375195e-10, 
                                            -5.61346232e-14, -8.00533100e+03, 2.22041168e+01};

const std::vector<double> a_prime_CH4_low = {-3.71537240e+05,  6.36843991e+03, -4.21090244e+01,
                                            1.10931447e-01, -1.73395947e-04,  1.52615054e-07, 
                                            -5.64354811e-11, -3.76219727e+04, 2.22877252e+02};
const std::vector<double> a_prime_CH4_hig = {2.27281672e+06, -9.56418553e+03,  7.59149994e+00, 
                                            -9.81636909e-04,-6.23317903e-07,  3.49499945e-10, 
                                            -4.86686647e-14,  5.06568957e+04, -7.27351090e+01};

const std::vector<double> a_prime_H2O_low = {-6.31350233e+04,  1.13414166e+03, -7.84242484e+00,  
                                            1.77729128e-02, -2.45367931e-05,  1.79953325e-08, 
                                            -5.22454675e-12, -3.40265003e+04, 3.84314488e+01};
const std::vector<double> a_prime_H2O_hig = {9.93128806e+05, -2.74796323e+03,  7.60880232e-01,
                                            4.05825396e-04, -2.00208030e-07,  2.46298394e-11, 
                                            -8.05713078e-16, -1.07376349e+04, -1.44689563e+01};

// NASA7 Polynomials from BURCAT for PH3 Equilibrium (unprimed - use carefully; for G(T)):
const std::vector<double> H3PO4_nasa_coeff_low = {-1.82480022E+00, 6.92983740E-02, -1.01448758E-04, 
                                                  6.91910155E-08, -1.77690951E-11, -1.36317149E+05, 
                                                  3.30354162E+01};
const std::vector<double> H3PO4_nasa_coeff_hig = {1.38520189E+01, 4.93643885E-03, -1.55822515E-06, 
                                                  2.29340638E-10, -1.28366239E-14, -1.39420678E+05, 
                                                  -4.22913846E+01};
const std::vector<double> PH3_nasa_coeff_low   = {4.17009763E+00, -5.06487157E-03, 2.86027846E-05, 
                                                  -3.13123782E-08, 1.13447768E-11, 2.03144445E+02, 
                                                  2.02004617E+00};
const std::vector<double> PH3_nasa_coeff_hig   = {3.71229298E+00, 5.85959002E-03, -2.16607791E-06, 
                                                  3.56195511E-10, -2.15913467E-14, -1.88863997E+02, 
                                                  1.92781913E+00};
const std::vector<double> H2O_nasa_coeff_low   = {4.198640560E+00,  -2.036434100E-03, 6.520402110E-06,  
                                                  -5.487970620E-09,   1.771978170E-12, -3.029372670E+04, 
                                                   -8.490322080E-01};
const std::vector<double> H2O_nasa_coeff_hig   = { 3.033992490E+00,   2.176918040E-03, -1.640725180E-07,  
                                                  -9.704198700E-11,  1.682009920E-14, -3.000429710E+04,
                                                   4.966770100E+00};
const std::vector<double> H2_nasa_coeff_low   = {2.344331120E+00,   7.980520750E-03, -1.947815100E-05, 
                                                2.015720940E-08,  -7.376117610E-12, -9.179351730E+02, 
                                                6.830102380E-01};
const std::vector<double> H2_nasa_coeff_hig   = {3.337279200E+00,  -4.940247310E-05, 4.994567780E-07, 
                                                -1.795663940E-10,   2.002553760E-14, -9.501589220E+02, 
                                                -3.205023310E+00};

// Additional NASA7 polynomials from BURCAT for PH3 timescale
const std::vector<double> H_nasa_coeff_low   = {2.500000000E+00,   7.053328190E-13, -1.995919640E-15, 
                                               2.300816320E-18,  -9.277323320E-22, 2.547365990E+04,  
                                               -4.466828530E-01};
const std::vector<double> H_nasa_coeff_hig   = {2.500000010E+00,  -2.308429730E-11, 1.615619480E-14,
                                               -4.735152350E-18,   4.981973570E-22, 2.547365990E+04, 
                                               -4.466829140E-01};
const std::vector<double> PO2_nasa_coeff_low   = {3.70338274E+00, 1.99856635E-03, 8.62444075E-06, 
                                                -1.34479846E-08, 5.58470721E-12, -3.51050019E+04, 
                                                8.53959173E+00};
const std::vector<double> PO2_nasa_coeff_hig   = {5.27587261E+00, 1.80638038E-03, -7.50028350E-07, 
                                                1.40062873E-10, -9.17506743E-15, -3.57348111E+04,
                                                -5.74241509E-01};
const std::vector<double> HOPO2_nasa_coeff_low = {1.63070941E+00, 2.49990293E-02, -2.23292246E-05, 
                                                  5.83782168E-09, 1.06880806E-12, -8.72608520E+04, 
                                                  1.81667382E+01};
const std::vector<double> HOPO2_nasa_coeff_hig = {9.07543004E+00, 3.08252827E-03, -1.12339846E-06, 
                                                  1.83697534E-10, -1.11135239E-14, -8.91864155E+04, 
                                                  -1.97912936E+01 };


//  GeH4 Equilibrium using Fegley & Lodders, 1994:
const std::vector<double> xgeh4FL_powers = {-41.451744864332106, -41.0785482543876, -40.7053516444431, -40.3321550344986, 
                                            -39.9589584245541, -39.39916350963733, -38.65277028974833, -38.09297537483157, 
                                            -37.34658215494257, -36.41359063008131, -35.6671974101923, -34.92080419030329, 
                                            -33.987812665442036, -33.05482114058078, -32.12182961571952, -31.18883809085826, 
                                            -30.255846565997004, -29.509453346107996, -28.389863516274485, -27.45687199141323, 
                                            -26.52388046655197, -25.59088894169071, -24.4712991118572, -23.165110977051437, 
                                            -21.858922842245676, -20.552734707439914, -19.784051579218, -19.437355378438777, 
                                            -18.964677728774003, -18.491802306034053, -17.987480964341515, -17.420267784743782, 
                                            -16.947392362003832, -16.474912485414233, -16.065522219804997, -15.498309040207264, 
                                            -14.994185471589896, -14.490259676047701, -13.986136107430331, -13.51365623084073,
                                            -13.072622273203724, -12.663825326820005, -12.19233431560627, -11.595257401657337, 
                                            -11.091924925340663, -10.526293930344321, -10.055198465280935, -9.678440757075332, 
                                            -9.333128967822319, -9.145837865632974, -9.086110396930561, -9.027371793604019, 
                                            -8.9684354172023, -8.878448667998338, -8.78707750726816, -8.632814508632794, 
                                            -8.445918952593795, -8.259023396554792, -8.10377153254356, -7.980361133635268};

const std::vector<double> T_FL  = {86.25592417061611, 90.99526066350711, 95.73459715639808, 100.47393364928911, 
                                  105.21327014218011, 112.32227488151659, 121.80094786729859, 128.9099526066351, 
                                  138.3886255924171, 150.23696682464453, 159.71563981042652, 169.19431279620852, 
                                  181.042654028436, 192.8909952606635, 204.73933649289103, 216.5876777251185, 
                                  228.43601895734594, 237.91469194312793, 252.13270142180096, 263.9810426540284, 
                                  275.8293838862559, 287.6777251184834, 301.8957345971564, 318.4834123222749, 
                                  335.07109004739334, 351.65876777251185, 361.42050703802516, 365.82330682911646, 
                                  372.2979900709245, 378.33388827990257, 385.2697909562978, 393.92782573467537, 
                                  400.6918180431376, 409.25628396945325, 417.36123005301306, 427.54075073088126, 
                                  437.31391368838473, 448.48514414580217, 459.251353898205, 470.5373262111708, 
                                  481.2933992661092, 495.9770467280551, 515.3153520416222, 543.0998139112735, 
                                  564.8881825722897, 596.7596786233797, 628.662368906931, 657.9654510556622, 
                                  687.8875174733606, 718.1198760123704, 764.4873166288049, 833.2145490798682, 
                                  911.6076137075584, 1020.7736100456253, 1117.1366234913971, 1248.1609479140955, 
                                  1368.731643273825, 1515.0870464147667, 1696.360723661587, 1944.3007796352736};

// NASA7 polynomials for GeH4 timescale:     
const std::vector<double> GeH4_calc_coeff_low = {2.54992789E+00, 7.13885765E-03, 1.43758337E-05,
                                                -2.33592977E-08, 9.65676013E-12, 9.69756465E+03, 9.02678812E+00};
const std::vector<double> GeH4_calc_coeff_hig = {5.41474159E+00, 7.24155154E-03, -2.71818301E-06,
                                                4.51535021E-10, -2.75635275E-14, 8.46356611E+03, -7.83419271E+00};                                              
const std::vector<double> GeH2_calc_coeff_low = {3.7766671605834596, 0.0009625978066155319, 3.939252214813524e-06, 
                                                -3.8051748899537e-09, 1.0224338840428087e-12, 28811.902151701965, 
                                                  5.704596034306216};
// GeH2 coefficients have more error in the high T limit (as demonstrated by CH2 and SiH2)
const std::vector<double> GeH2_calc_coeff_hig = GeH2_calc_coeff_low;

// Chemical parameter initializations:
double kBoltz    = 1.38e-16; // erg/K (cgs)
double R_uni     = 8.3145; // J/mol/K
double X_H2      = 0.864;
double X_CH4     = 2.37e-3 * X_H2; // (Wong+2004)
double R_dry_J   = R_uni / 2.22 / 1e-3; // J/kg/K using conventional mean_mu for a dry Jupiter

double X_H2O     = 9.61e-4;        // default to use in no vapor case.
double X_P       = 7.0e-7 * X_H2;  // enriched by 1.4
double X_Ge      = 3.93e-8 * X_H2; // enriched by 4.4
double X_H2S     = 8.9e-5 * X_H2;  // enriched by 3.38
double E_Ge      = 4.4;

// Molecular weights for conversion between mass and mole fractions:
double mu_H2O    = 18.01528; // molecular weight of H2O.
double mu_CO     = 28.01;    // molecular weight of CO.
double mu_PH3    = 33.99758; // molecular weight of PH3.
double mu_GeH4   = 76.62;    // molecular weight of GeH4.

double temperature;
// double temp_vir;
double mean_mu;
double n_particle_density;

// X_H2O_05sol  = 4.805e-4  // (0.5 x Solar)
// X_H2O_1sol   = 9.61e-4  // (1.0 x Solar)
// X_H2O_25sol  = 2.4e-3   // (2.5 x Solar)
// X_H2O_5sol   = 4.805e-3 // (5.0 x Solar)

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  if (NVAPOR > 0) {
    AllocateUserOutputVariables(5 + NVAPOR);
    SetUserOutputVariableName(0, "temp");
    SetUserOutputVariableName(1, "theta");
    SetUserOutputVariableName(2, "pres");
    SetUserOutputVariableName(3, "thetav");
    SetUserOutputVariableName(4, "mse");
    for (int n = 1; n <= NVAPOR; ++n) {
      std::string name = "rh" + std::to_string(n);
      SetUserOutputVariableName(4 + n, name.c_str());
    }
  } else {
    AllocateUserOutputVariables(3);
    SetUserOutputVariableName(0, "temp");
    SetUserOutputVariableName(1, "theta");
    SetUserOutputVariableName(2, "pres");
  }
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  auto pthermo = Thermodynamics::GetInstance();

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        user_out_var(0, k, j, i) = pthermo->GetTemp(this, k, j, i);
        user_out_var(1, k, j, i) = pthermo->PotentialTemp(this, P0, k, j, i);
        user_out_var(2, k, j, i) = phydro->w(IPR, k, j, i);

        if (NVAPOR > 0) {
          // theta_v
          user_out_var(3, k, j, i) =
              user_out_var(1, k, j, i) * pthermo->RovRd(this, k, j, i);
          // mse
          user_out_var(4, k, j, i) =
              pthermo->MoistStaticEnergy(this, grav * pcoord->x1v(i), k, j, i);
          // relative humidity
          for (int n = 1; n <= NVAPOR; ++n)
            user_out_var(4 + n, k, j, i) =
                pthermo->RelativeHumidity(this, n, k, j, i);
        }
      }
}

// Function to perform linear interpolation.
// Given a value x and arrays X and T, it finds the appropriate interval in X and interpolates the corresponding value in T.
double linearInterpolate(const std::vector<double>& X, const std::vector<double>& T, double x) {
    // Check if x is out of bounds of the X array
    if (x <= X.front()) return T.front();
    if (x >= X.back()) return T.back();
    // Find the right place in the array (lower bound of the interval)
    auto lower = std::lower_bound(X.begin(), X.end(), x);
    int index = lower - X.begin();
    // Edge case handling: if x matches the last element, interpolate with the second last
    if (index == X.size()) {
        index--;
    }
    // Calculate the interpolation
    double x1 = X[index - 1], x2 = X[index];
    double t1 = T[index - 1], t2 = T[index];
    // The interpolation formula
    double t = t1 + (x - x1) * (t2 - t1) / (x2 - x1);
    return t;
}


double EnthalpyPerRT(const std::vector<double> &a, double T){
        
  double enthalpy_per_RT = (
                      -a[0]/pow(T,2) + a[1]*log(T)/T + a[2] + a[3]*T/2
                      + (a[4]*pow(T,2))/3 + (a[5]*pow(T,3))/4 + (a[6]*pow(T,4))/5
                      + a[7]/T
                        );
  return enthalpy_per_RT;
}

double EntropyPerR(const std::vector<double> &a, double T){
    
  double entropy_per_R = (
                      -(a[0]/pow(T,2))/2 - a[1]/T + a[2]*log(T) + a[3]*T
                      + (a[4]*pow(T,2))/2 + (a[5]*pow(T,3))/3 + (a[6]*pow(T,4))/4
                      + a[8]
                        );
  return entropy_per_R;
}


double EnthalpyPerRT_NASA7(const std::vector<double> &a, double T){
        
  double enthalpy_per_RT = (
                            a[0] + a[1]*T/2 + (a[2]*pow(T,2))/3 + (a[3]*pow(T,3))/4
                            + (a[4]*pow(T,4))/5 + a[5]/T
                             );
  return enthalpy_per_RT;
}

double EntropyPerR_NASA7(const std::vector<double> &a, double T){
    
  double entropy_per_R = (
                            a[0]*log(T) + a[1]*T + (a[2]*pow(T,2))/2 + 
                            (a[3]*pow(T,3))/3 + (a[4]*pow(T,4))/4 + a[6]
                             );
  return entropy_per_R;
}

double GibbsPerRT(const std::vector<double> &a, double T){
    
  return EnthalpyPerRT(a, T) - EntropyPerR(a, T); // unitless
}

double GibbsPerRT_NASA7(const std::vector<double> &a, double T){
    
  return EnthalpyPerRT_NASA7(a, T) - EntropyPerR_NASA7(a, T); // unitless
}

double KeqCO(double T){

  double delta_Gibbs_f_CO;
  double delta_Gibbs_f_CH4;
  double delta_Gibbs_f_H2O;
  if (T < 1000){
    delta_Gibbs_f_CO  = GibbsPerRT(a_prime_CO_low, T);
    delta_Gibbs_f_CH4 = GibbsPerRT(a_prime_CH4_low, T);
    delta_Gibbs_f_H2O = GibbsPerRT(a_prime_H2O_low, T);
  }
  else{
    delta_Gibbs_f_CO  = GibbsPerRT(a_prime_CO_hig, T);
    delta_Gibbs_f_CH4 = GibbsPerRT(a_prime_CH4_hig, T);
    delta_Gibbs_f_H2O = GibbsPerRT(a_prime_H2O_hig, T);
  }
  return exp(-(delta_Gibbs_f_CO - delta_Gibbs_f_CH4 - delta_Gibbs_f_H2O));
}

double KeqPH3(double T){

  double Gibbs_f_H3PO4;
  double Gibbs_f_H2;
  double Gibbs_f_PH3;
  double Gibbs_f_H2O;
  if (T < 1000){
    Gibbs_f_H3PO4  = GibbsPerRT_NASA7(H3PO4_nasa_coeff_low, T);
    Gibbs_f_H2     = GibbsPerRT_NASA7(H2_nasa_coeff_low, T);
    Gibbs_f_PH3    = GibbsPerRT_NASA7(PH3_nasa_coeff_low, T);
    Gibbs_f_H2O    = GibbsPerRT_NASA7(H2O_nasa_coeff_low, T);
  }
  else{
    Gibbs_f_H3PO4  = GibbsPerRT_NASA7(H3PO4_nasa_coeff_hig, T);
    Gibbs_f_H2     = GibbsPerRT_NASA7(H2_nasa_coeff_hig, T);
    Gibbs_f_PH3    = GibbsPerRT_NASA7(PH3_nasa_coeff_hig, T);
    Gibbs_f_H2O    = GibbsPerRT_NASA7(H2O_nasa_coeff_hig, T);
  }
  return exp(-(Gibbs_f_H3PO4 + 4*Gibbs_f_H2 - Gibbs_f_PH3 - 4*Gibbs_f_H2O));
}

double KeqPO2(double T){

  double Gibbs_f_PO2;
  double Gibbs_f_H;
  double Gibbs_f_HOPO2;
  double Gibbs_f_H2O;
  if (T < 1000){
    Gibbs_f_PO2   = GibbsPerRT_NASA7(PO2_nasa_coeff_low, T);
    Gibbs_f_H     = GibbsPerRT_NASA7(H_nasa_coeff_low, T);
    Gibbs_f_HOPO2 = GibbsPerRT_NASA7(HOPO2_nasa_coeff_low, T);
    Gibbs_f_H2O   = GibbsPerRT_NASA7(H2O_nasa_coeff_low, T);
  }
  else{
    Gibbs_f_PO2   = GibbsPerRT_NASA7(PO2_nasa_coeff_hig, T);
    Gibbs_f_H     = GibbsPerRT_NASA7(H_nasa_coeff_hig, T);
    Gibbs_f_HOPO2 = GibbsPerRT_NASA7(HOPO2_nasa_coeff_hig, T);
    Gibbs_f_H2O   = GibbsPerRT_NASA7(H2O_nasa_coeff_hig, T);
  }
  return exp(-(Gibbs_f_HOPO2 + Gibbs_f_H - Gibbs_f_PO2 - Gibbs_f_H2O));
}

double KeqGeH2(double T){
  double Gibbs_f_GeH2;
  double Gibbs_f_H2;
  double Gibbs_f_GeH4;
  if (T < 1000){
    Gibbs_f_GeH2 = GibbsPerRT_NASA7(GeH2_calc_coeff_low, T);
    Gibbs_f_H2   = GibbsPerRT_NASA7(H2_nasa_coeff_low, T);
    Gibbs_f_GeH4 = GibbsPerRT_NASA7(GeH4_calc_coeff_low, T);
  }
  else{
    Gibbs_f_GeH2 = GibbsPerRT_NASA7(GeH2_calc_coeff_hig, T);
    Gibbs_f_H2   = GibbsPerRT_NASA7(H2_nasa_coeff_hig, T);
    Gibbs_f_GeH4 = GibbsPerRT_NASA7(GeH4_calc_coeff_hig, T);
  }
  return exp(-(Gibbs_f_GeH2 + Gibbs_f_H2 - Gibbs_f_GeH4));
}

double X_CO_Eq(double press, double Temp, double X_CH4, double X_H2O, double X_H2){

  double K_eq_CO = KeqCO(Temp);
  return X_CH4 * X_H2O * K_eq_CO / pow(X_H2, 3) / pow((press/1E5), 2);

}

double X_PH3_Eq(double Temp, double X_H2O, double X_H2){

  double K_eq_PH3 = KeqPH3(Temp);
  return X_P / (1 + (K_eq_PH3 * pow(X_H2O/X_H2, 4)));

}

void Forcing(MeshBlock *pmb, Real const time, Real const dt,
             AthenaArray<Real> const &w, const AthenaArray<Real> &prim_scalar,
             AthenaArray<Real> const &bcc, AthenaArray<Real> &du,
             AthenaArray<Real> &cons_scalar) {
  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;
  auto pthermo = Thermodynamics::GetInstance();
  auto phydro = pmb->phydro;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        if (use_polar_beta) {
          Real x2 = pmb->pcoord->x2v(j);
          Real x3 = pmb->pcoord->x3v(k);
          Real dist = sqrt(x2*x2 + x3*x3);

          Real fcor = -omega*sqr(dist/radius);
          du(IM2,k,j,i) += dt*fcor*w(IDN,k,j,i)*w(IM3,k,j,i);
          du(IM3,k,j,i) -= dt*fcor*w(IDN,k,j,i)*w(IM2,k,j,i);
        }
        
        temperature  = pthermo->GetTemp(w.at(k,j,i));
        
        // mean molecular weight in g/mol:
        mean_mu      = 1E3 * w(IDN,k,j,i) * R_uni * temperature / w(IPR,k,j,i); // technically uses virtual temperature
        // particle number density in molecules/cm3:
        n_particle_density =  10.0 * w(IDN,k,j,i) * R_dry_J / kBoltz; // --> kg/m3 * (J/kg/K) / (ergs/K) * 1e7 * 1e-6
        // Chemical Tracer method using linear relaxation scheme:
        if (NVAPOR > 0){
          // if moist case is used, set X_H2O to model's evolving water content.
          X_H2O = w(1,k,j,i) * (mean_mu/mu_H2O);
        }

        if (NSCALARS > 0){
          if (NSCALARS == 1){
            for (int n=0; n<NSCALARS; ++n){
              // Passive Scalar CO
              if (s_ind == 1) {
                double tau = tau_scale_factor * 5E-6 * exp(2.8E4/temperature) * pow(w(IPR,k,j,i)/1E5, -1.29);
//                double C_CO_equi = X_CO_Eq(w(IPR,k,j,i), temperature, 
//                                    X_CH4, X_H2O, X_H2) * w(IDN,k,j,i)
//                                    * (mu_CO/mean_mu); // conversion factor to Mass fraction;
                //pmb->pscalars->s(n,k,j,i) -= dt  * (pmb->pscalars->s(n,k,j,i) - C_CO_equi) / tau;

                double C_CO_equi = X_CO_Eq(w(IPR,k,j,i), temperature,
                                    X_CH4, X_H2O, X_H2)
                                    * (mu_CO/mean_mu); // conversion factor to Mass fraction;
//                pmb->pscalars->s(n,k,j,i) -= dt  * (pmb->pscalars->r(n,k,j,i) - C_CO_equi) / tau * w(IDN,k,j,i);
               cons_scalar(NCLOUD + n, k, j, i) -= dt  * (prim_scalar(NCLOUD + n, k, j, i) - C_CO_equi) / tau * w(IDN,k,j,i);
              }

              // Passive Scalar PH3
              if (s_ind == 2) {
                double khf   = 2.35e-12 * exp(-1.067e4/temperature);
                double K_PO2 = KeqPO2(temperature);
                
                double tau =  pow(w(IPR,k,j,i)/1E5, 1.5) * pow(X_H2, 3.5) / 
                              (khf * K_PO2 * n_particle_density * pow(X_H2O, 3));
                double C_PH3_equi = X_PH3_Eq(temperature, 
                                    X_H2O, X_H2) * w(IDN,k,j,i)
                                    * (mu_PH3/mean_mu); // conversion factor to Mass fraction;
//                pmb->pscalars->s(n,k,j,i) -= dt  * (pmb->pscalars->s(n,k,j,i) - C_PH3_equi) / tau;
              }

              // Passive Scalar for GeH4
              if (s_ind == 3){
                double k4b = 3.57e-14 * pow(temperature, 0.7)*exp(-4956./temperature);
                double K_GeH2 = KeqGeH2(temperature);
                
                // No discontinuity occurs with 1E3. We might need to add a minimum timescale to limit numerical instability
                double tau = fmax(10., (w(IPR,k,j,i)/1E5) * X_H2 / (k4b * K_GeH2 * n_particle_density * X_H2S));
                double C_GeH4_equi = pow(10, linearInterpolate(T_FL, xgeh4FL_powers, temperature)) * 
                                      (mu_GeH4/mean_mu) * w(IDN,k,j,i) * X_H2 * E_Ge/6.0;
                
//                pmb->pscalars->s(n,k,j,i) -= dt  * (pmb->pscalars->s(n,k,j,i) - C_GeH4_equi) / tau;
              }
            }
          }
          if (NSCALARS == 2){
            for (int n=0; n<NSCALARS; ++n){
              // Passive scalar PH3 with NSCALAR=2 configuration
              if (n == 0) {
                double khf   = 2.35e-12 * exp(-1.067e4/temperature);
                double K_PO2 = KeqPO2(temperature);
                
                double tau =  pow(w(IPR,k,j,i)/1E5, 1.5) * pow(X_H2, 3.5) / 
                              (khf * K_PO2 * n_particle_density * pow(X_H2O, 3));
                double C_PH3_equi = X_PH3_Eq(temperature, 
                                    X_H2O, X_H2) * w(IDN,k,j,i)
                                    * (mu_PH3/mean_mu); // conversion factor to Mass fraction;
//                pmb->pscalars->s(n,k,j,i) -= dt  * (pmb->pscalars->s(n,k,j,i) - C_PH3_equi) / tau;
              }
              // Passive scalar GeH4 with NSCALAR=2 configuration
              if (n == 1) {
                double k4b = 3.57e-14 * pow(temperature, 0.7)*exp(-4956./temperature);
                double K_GeH2 = KeqGeH2(temperature);
                
                // No discontinuity occurs with 1E3. We might need to add a minimum timescale to limit numerical instability
                double tau = fmax(10., (w(IPR,k,j,i)/1E5) * X_H2 / (k4b * K_GeH2 * n_particle_density * X_H2S));
                double C_GeH4_equi = pow(10, linearInterpolate(T_FL, xgeh4FL_powers, temperature)) * 
                                      (mu_GeH4/mean_mu) * w(IDN,k,j,i) * X_H2 * E_Ge/6.0;
                
//                pmb->pscalars->s(n,k,j,i) -= dt  * (pmb->pscalars->s(n,k,j,i) - C_GeH4_equi) / tau;
              }
            }
          }
        }
      }

}


void Mesh::InitUserMeshData(ParameterInput *pin)
{
  grav = - pin->GetReal("hydro", "grav_acc1");
  P0 = pin->GetReal("problem", "P0");
  T0 = pin->GetReal("problem", "T0");
  tau_scale_factor = pin->GetOrAddReal("problem", "tau_scale_factor", 1.);
  
  // Scalar indexing: CO=1, PH3=2, GeH4=3; defaults to CO by design
  s_ind = pin->GetOrAddInteger("hydro", "s_ind", 1);
  Z0 = pin->GetOrAddReal("problem", "Z0", 0.);
  nu_scalar_iso = pin->GetOrAddReal("problem", "nu_scalar_iso", 0.25);
  use_polar_beta = pin->GetOrAddBoolean("problem", "use_polar_beta", false);
  if (use_polar_beta) {
    radius = pin->GetReal("problem", "radius");
    omega = pin->GetReal("hydro", "OmegaZ");
  }
  Tmin = pin->GetReal("hydro", "min_tem");

  EnrollUserExplicitSourceFunction(Forcing);
  // if (NSCALARS > 0) {
  //   EnrollUserExplicitSourceFunction(ChemTracer);
  // }
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  //ATHENA_LOG("convection");
  Real gamma = pin->GetReal("hydro", "gamma");


  // construct a 1D pseudo-moist adiabat with given relative humidity
  Real x1min = pmy_mesh->mesh_size.x1min;
  Real x1max = pmy_mesh->mesh_size.x1max;
  Real x1rat = pmy_mesh->mesh_size.x1rat;

  Real dz, **w1, *z1, *p1, *t1;
  if (x1rat != 1.0) {
    dz = (x1max - x1min)*(x1rat - 1.)/(pow(x1rat, pmy_mesh->mesh_size.nx1) - 1.);
    dz = std::min(dz, dz*pow(x1rat, pmy_mesh->mesh_size.nx1))/2.;
  } else {
    dz = (x1max - x1min)/pmy_mesh->mesh_size.nx1/2.;
  }
  int nx1 = (int)((x1max - x1min)/dz);
  NewCArray(w1, nx1, NHYDRO+2*NVAPOR);
  z1 = new Real [nx1];
  p1 = new Real [nx1];
  t1 = new Real [nx1];

  // 1.1 estimate surface temperature and pressure
  Real Rd = pthermo->GetRd();
  Real cp = gamma/(gamma - 1.)*Rd;
  Real Ts = T0 - grav/cp*(x1min - Z0);
  Real Ps = P0*pow(Ts/T0, cp/Rd);
  int max_iter = 200, iter = 0;


  for (int n = 1; n <= NVAPOR; ++n) {
    Real qv = pin->GetReal("problem", "qvapor" + std::to_string(n))/1.E3;
    w1[0][n] = qv;
  }
  z1[0] = x1min;
  for (int i = 1; i < nx1; ++i)
    z1[i] = z1[i-1] + dz;

  Real t0, p0;
  if (Globals::my_rank == 0)
    std::cout << "- request T = " << T0 << " P = " << P0 << " at Z = " << Z0 << std::endl;
  while (iter++ < max_iter) {
    pthermo->ConstructAdiabat(w1, Ts, Ps, grav, dz, nx1, Adiabat::pseudo);
//    pthermo->ConstructAtmosphere(w1, Ts, Ps, grav, dz, nx1, Adiabat::pseudo, 0.);

    // 1.2 replace adiabatic atmosphere with isothermal atmosphere if temperature is too low
    int ii = 0;
    for (; ii < nx1-1; ++ii)
      if (pthermo->Temp(w1[ii]) < Tmin) break;
    Real Tv = w1[ii][IPR]/(w1[ii][IDN]*Rd);
    for (int i = ii; i < nx1; ++i) {
      w1[i][IPR] = w1[ii][IPR]*exp(-grav*(z1[i] - z1[ii])/(Rd*Tv));
      w1[i][IDN] = w1[i][IPR]/(Rd*Tv);
      for (int n = 1; n <= NVAPOR; ++n)
        w1[i][n] = w1[ii][n];
    }

    // 1.3 find TP at z =  Z0
    for (int i = 0; i < nx1; ++i) {
      p1[i] = w1[i][IPR];
      t1[i] = pthermo->Temp(w1[i]);
    }
    p0 = interp1(Z0, p1, z1, nx1);
    t0 = interp1(Z0, t1, z1, nx1);

    Ts += T0 - t0;
    Ps *= P0/p0;
    if ((fabs(T0 - t0) < 0.01) && (fabs(P0/p0 - 1.) < 1.E-4)) break;
    if (Globals::my_rank == 0)
      std::cout << "- iteration #" << iter << ": " << "T = " << t0 << " P = " << p0 << std::endl;
  }

  if (iter > max_iter) {
    std::stringstream msg;
    msg << "### FATAL ERROR in problem generator"
        << std::endl << "maximum iteration reached."
        << std::endl << "T0 = " << t0
        << std::endl << "P0 = " << p0;
    ATHENA_ERROR(msg);
  }

  // setup initial condition
  srand(Globals::my_rank + time(0));
  for (int i = is; i <= ie; ++i) {
    Real buf[NHYDRO+2*NVAPOR];
    interpn(buf, &pcoord->x1v(i), *w1, z1, &nx1, 1, NHYDRO+2*NVAPOR);
    buf[IVX] = buf[IVY] = buf[IVZ] = 0.;

    // set gas concentration
    for (int n = 0; n < NHYDRO; ++n)
      for (int k = ks; k <= ke; ++k)
        for (int j = js; j <= je; ++j)
          phydro->w(n,k,j,i) = buf[n];

  
  if (NSCALARS > 0) {    
    if (NSCALARS == 1){ // If configured to single scalar.  Here, s_ind is important.
      for (int n=0; n<NSCALARS; ++n) {
        for (int i = is; i <= ie; ++i) {
          for (int k = ks; k <= ke; ++k){
            for (int j = js; j <= je; ++j){
              temperature = pthermo->Temp(phydro->w.at(k,j,i));
              mean_mu = 1E3 * phydro->w(IDN,k,j,i) * R_uni * temperature / phydro->w(IPR,k,j,i);
              // Set water mole fraction appropriately:
              if (NVAPOR > 0){
                //if moist case is used, set X_H2O to the evolving water content.
                X_H2O = phydro->w(1,k,j,i) * (mean_mu/mu_H2O);
              }

              // s_ind = 1 for CO initial equilibrium condition 
              if (s_ind == 1) {
      //        pscalars->s(n,k,j,i) =  X_CO_Eq(phydro->w(IPR,k,j,i), temperature, 
      //                                X_CH4, X_H2O, X_H2) * (mu_CO/mean_mu) * phydro->w(IDN,k,j,i);
              pscalars->r(n,k,j,i) =  X_CO_Eq(phydro->w(IPR,k,j,i), temperature,
                                      X_CH4, X_H2O, X_H2) * (mu_CO/mean_mu);
              }
              // s_ind = 2 for PH3 initial equilibrium condition
              if (s_ind == 2) {
              pscalars->s(n,k,j,i) =  X_PH3_Eq(temperature, 
                                      X_H2O, X_H2) * (mu_PH3/mean_mu) * phydro->w(IDN,k,j,i);
              }
              // s_ind = 3 for GeH4 initial equilibrium condition
              if (s_ind == 3) {
              // Dividing out 6.0 due to the enrichment of Ge used in Fegley & Lodders, 1994.
              pscalars->s(n,k,j,i) =  pow(10, linearInterpolate(T_FL, xgeh4FL_powers, temperature)) * 
                                      (mu_GeH4/mean_mu) * phydro->w(IDN,k,j,i) * X_H2 * E_Ge/6.0;
              }
            }
          }
        }
      }
    }
    if (NSCALARS == 2){ // If configured to 2 scalars.  CO is ignored and PH3 and GeH4 are modeled.
      for (int n=0; n<NSCALARS; ++n) {
        for (int i = is; i <= ie; ++i) {
          for (int k = ks; k <= ke; ++k){
            for (int j = js; j <= je; ++j){
              temperature = pthermo->Temp(phydro->w.at(k,j,i));
              mean_mu = 1E3 * phydro->w(IDN,k,j,i) * R_uni * temperature / phydro->w(IPR,k,j,i);
              // Set water mole fraction appropriately:
              if (NVAPOR > 0){
                //if moist case is used, set X_H2O to the evolving water content.
                X_H2O = phydro->w(1,k,j,i) * (mean_mu/mu_H2O);
              }
              // n = 0 for PH3 initial equilibrium condition if NSCALARS == 2
              if (n == 0) {
              pscalars->s(n,k,j,i) =  X_PH3_Eq(temperature, 
                                      X_H2O, X_H2) * (mu_PH3/mean_mu) * phydro->w(IDN,k,j,i);
              }
              // n = 1 for GeH4 initial equilibrium condition if NSCALARS == 2
              if (n == 1) {
              // Dividing out 6.0 due to the enrichment of Ge used in Fegley & Lodders, 1994.
              pscalars->s(n,k,j,i) =  pow(10, linearInterpolate(T_FL, xgeh4FL_powers, temperature)) * 
                                      (mu_GeH4/mean_mu) * phydro->w(IDN,k,j,i) * X_H2 * E_Ge/6.0;
              }
            }
          }
        }
      }
    }
  }
  // std::cout << "Interpolated at 1050.0 K = " << pow(10, linearInterpolate(T_FL, xgeh4FL_powers, 1050.0)) << std::endl;
    // // Initialize some constant tracer 1/3 or the height in the box:
    // if (NSCALARS > 0) {
    //   for (int n=0; n<NSCALARS; ++n) {
    //     if (i==ie/3){
    //       for (int k = ks; k <= ke; ++k){
    //         for (int j = js; j <= je; ++j){
    //           pscalars->s(n,k,j,i) = X_CO_Eq(phydro->w(IPR,k,j,i), pthermo->Temp(phydro->w.at(k,j,i)), 
    //                                 1.0, 1.0, 1.0) * phydro->w(IDN,k,j,i);  
    //         }
    //       }
    //     }
    //     else{
    //       for (int k = ks; k <= ke; ++k){
    //         for (int j = js; j <= je; ++j){
    //           pscalars->s(n,k,j,i) = 0.0;
    //         }
    //       }
    //     }
    //   }
    // }

    /* set cloud concentration
    Particles *pp = ppart;
    for (int n = 0; n < NVAPOR; ++n) {
      if (pp == nullptr) {
        msg << "### FATAL ERROR in problem generator"
            << std::endl << "Vapor #" << n
            << "does not associate with any particle.";
        ATHENA_ERROR(msg);
      }
      for (int k = ks; k <= ke; ++k)
        for (int j = js; j <= je; ++j)
          pp->u(0,k,j,i) = buf[NHYDRO+n] + buf[NHYDRO+NVAPOR+n];
      pp = pp->next;
    }*/

    // add noise
    for (int k = ks; k <= ke; ++k){
      for (int j = js; j <= je; ++j){
        phydro->w(IV1,k,j,i) = 0.01*(1.*rand()/RAND_MAX - 0.5);
        }
      }
  }

  //if (ppart != nullptr)
  //ppart->Initialize();
  //pphy->Initialize(phydro->w);
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);
  peos->PassiveScalarPrimitiveToConserved(pscalars->r, phydro->u, pscalars->s, pcoord, is,
      ie, js, je, ks, ke);

  // // set chemical tracer
  // for (int k = ks; k <= ke; ++k)
  //   for (int j = js; j <= je; ++j)
  //     for (int i = is; i <= ie; ++i) {
  //       pscalars->s(0,k,j,i) = 1.;
  //       pscalars->r(0,k,j,i) = 1.;
  //     }

  FreeCArray(w1);
  delete[] z1;
  delete[] p1;
  delete[] t1;
}
