// C/C++
#include <cassert>
#include <sstream>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// utils
#include <utils/fileio.hpp>

// climath
#include <climath/interpolation.h>

#include <climath/root.hpp>

// snap
#include "molecules.hpp"

#define PRECISION 1.E-8

std::ostream& operator<<(std::ostream& os, Molecule const& mol) {
  os << "name: " << mol.m_name << std::endl
     << "molecular weight: " << mol.m_mu << " kg/mol" << std::endl
     << "heat capacity (cp): " << mol.m_cp << " J/(mol K)" << std::endl
     << "vaporization heat at triple point: " << mol.m_latent << " kJ/mol"
     << std::endl
     << "vapor pressure coefficients: " << mol.m_beta << " " << mol.m_gamma
     << std::endl
     << "standard entropy: " << mol.m_entropy << " J/(mol K)" << std::endl
     << "standard enthalpy: " << mol.m_enthalpy << " kJ/mol" << std::endl;

  return os;
}

void Molecule::LoadThermodynamicFile(std::string chemfile) {
  auto app = Application::GetInstance();
  auto full_path = app->FindResource(chemfile);

  std::stringstream inp(DecommentFile(full_path));
  double junk;

  inp >> m_name >> m_mu >> m_entropy >> m_enthalpy >> m_gibbs >> m_cp >> m_tr >>
      m_pr >> m_tc >> m_pc;
  inp >> m_nshomate;
  for (int i = 0; i < m_nshomate; i++) {
    inp >> m_shomate_sp.at(i) >> junk;
    for (int j = 0; j < NSHOMATE; j++) inp >> m_shomate.at(i * NSHOMATE + j);
  }
  m_shomate_sp.at(m_nshomate) = junk;

  inp >> m_cliq >> m_enliq >> m_csld >> m_ensld;

  m_mu *= 1.E-3;  // g/mol -> kg/mol
}

double Molecule::cp(double T) const {
  int i = locate(m_shomate_sp.data(), T, m_nshomate + 1);
  if (i == m_nshomate) i = m_nshomate - 1;
  if (i < 0) i = 0;

  double result;
  T /= 1.E3;
  result = m_shomate[i * NSHOMATE] + m_shomate[i * NSHOMATE + 1] * T +
           m_shomate[i * NSHOMATE + 2] * T * T +
           m_shomate[i * NSHOMATE + 3] * T * T * T +
           m_shomate[i * NSHOMATE + 4] / (T * T);

  return result;
}

double Molecule::enthalpy(double T) const {
  int i = locate(m_shomate_sp.data(), T, m_nshomate + 1);
  if (i == m_nshomate) i = m_nshomate - 1;
  if (i < 0) i = 0;

  double result;
  T /= 1.E3;
  result = m_shomate[i * NSHOMATE] * T +
           0.5 * m_shomate[i * NSHOMATE + 1] * T * T +
           1. / 3. * m_shomate[i * NSHOMATE + 2] * T * T * T +
           1. / 4. * m_shomate[i * NSHOMATE + 3] * T * T * T * T +
           -m_shomate[i * NSHOMATE + 4] / T + m_shomate[i * NSHOMATE + 5];

  return result;
}

double Molecule::entropy(double T) const {
  int i = locate(m_shomate_sp.data(), T, m_nshomate + 1);
  if (i == m_nshomate) i = m_nshomate - 1;
  if (i < 0) i = 0;

  double result;
  T /= 1.E3;
  result = m_shomate[i * NSHOMATE] * log(T) + m_shomate[i * NSHOMATE + 1] * T +
           0.5 * m_shomate[i * NSHOMATE + 2] * T * T +
           1. / 3. * m_shomate[i * NSHOMATE + 3] * T * T * T +
           -m_shomate[i * NSHOMATE + 4] / (2. * T * T) +
           m_shomate[i * NSHOMATE + 6];

  return result;
}

double Molecule::latent(double T) const {
  double result;
  result = m_latent + enthalpy(T) - enthalpy(m_tr) - m_cp * (T - m_tr) / 1.E3;

  return result;
}

double Molecule::isat_vapor_p(double P) const {
  double Tsat, Tmin, Tmax;
  Tmin = 50.;
  Tmax = 2000.;
  int error = root(Tmin, Tmax, PRECISION, &Tsat,
                   [this, P](double T) { return this->sat_vapor_p(T) - P; });
  if (error) {
    std::stringstream msg;
    msg << "Inverting saturation vapor pressure is not sucessfull" << std::endl;
    msg << "Pressure = " << P << " Pa" << std::endl;
    msg << "Tmin = " << Tmin << std::endl;
    msg << "Tmax = " << Tmax << std::endl;
    msg << "Saturation vapor pressure at Tmin = " << sat_vapor_p(Tmin)
        << std::endl;
    msg << "Saturation vapor pressure at Tmax = " << sat_vapor_p(Tmax)
        << std::endl;
    msg << *this << std::endl;
    throw RuntimeError("Molecule", msg.str());
  }

  return Tsat;
}
