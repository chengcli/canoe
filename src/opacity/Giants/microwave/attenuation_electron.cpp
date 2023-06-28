// C/C++
#include <cmath>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// harp
#include <harp/radiation.hpp>

// cli, 190801
inline double SQUARE(double x) { return x * x; }
inline double CUBE(double x) { return x * x * x; }
// end

// From F. Oyafuso
double attenuation_freefree_Reference(double freq, double P, double T)
// this is very crude
{
  int const Z = 1;
  double const fine_struct = 1.0 / 137.036;
  double const mc2 = 0.511e6;           // eV
  double const kT = (T / 300) * .0256;  // eV
  double const hbarc = 1973e-8;         // eV-cm
  double sigma_T = 6.65e-25;            // cm^2
  // double g_ff = 1; // gaunt factor
  double hbar_omega =
      4.13566733e-15 * freq * 1e9;  // eV (.6GHz) -- check that freq is in GHz

  double debroglie_thermal = sqrt(2.0 * M_PI * SQUARE(hbarc) / (mc2 * kT));
  double E_ionization = 5.986;  // eV
  double n_Al = 1.6e-6 * (P * 1e5) / (kT * 1.609e-19) * 1e-6;

  double ni =
      sqrt(2 * n_Al / pow(debroglie_thermal, 3) * exp(-E_ionization / kT));

  double kludge = 0.5;
  double alpha_ff = kludge * sqrt(32 * CUBE(M_PI) / 3.0) *
                    (SQUARE(Z) * fine_struct) * sqrt(mc2 / kT) *
                    (ni * CUBE(hbarc / hbar_omega)) * ni * sigma_T *
                    (1.0 - exp(-hbar_omega / kT));
  return alpha_ff;
}

// from C. Li (in units of 1/cm)
double attenuation_freefree_Chengli(double freq, double P, double T) {
  double ion_opacity_scale = 1.0;
  double ion_opacity_cutoff = 1600.;

  double alpha_ff;
  if (T < ion_opacity_cutoff)
    alpha_ff = 0;
  else {
    alpha_ff = ion_opacity_scale *
               exp(-20.957 + 9.51112E-3 * (1200. + T - ion_opacity_cutoff));
    // convert to cm^-1
    alpha_ff *= 1e-2;
  }
  return alpha_ff;
}

double attenuation_appleton_hartree_nomag(double freq_GHz, double P_bar,
                                          double T, double ne) {
  ne *= 1.E-6;                    // m^{-3} -> cm^{-3}
  double freq = freq_GHz * 1.E9;  // GHz -> Hz
  double P = P_bar * 1E6;         // bar -> Ba
  double e_cgs = 4.8E-10;         // esu
  double me_cgs = 9.11E-28;       // g

  double num = P / (Thermodynamics::kBoltz_cgs * T);  // cm^{-3}
  double colli_radius = 1.2E-8;                       // cm
  double lambda = Radiation::cLight_cgs / freq;       // wavelength (cm)
  double omega = 2. * M_PI * freq;  // angular frequency (rad/s)
  // electron plasma frequency (rad/s)
  double omega_pe = sqrt(4. * M_PI * ne * e_cgs * e_cgs / me_cgs);
  // electron collision frequency (?)
  double nu_e = 5. / 3. * (M_PI * colli_radius * colli_radius * num) *
                sqrt(3. * Thermodynamics::kBoltz_cgs * T / me_cgs);
  double X = omega_pe * omega_pe / (omega * omega);
  double Z = nu_e / omega;
  double n2real = 1. - X / (1. + Z * Z);
  double n2imag = X * Z / (1. + Z * Z);
  double Qg = n2real / n2imag;
  double alpha = 2. * M_PI / lambda / Qg;
  return alpha;  // 1/cm
}
