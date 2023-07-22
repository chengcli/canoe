#ifndef SRC_SNAP_THERMODYNAMICS_MOLECULES_HPP_
#define SRC_SNAP_THERMODYNAMICS_MOLECULES_HPP_

// C/C++
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

/** @brief Molecule stored the necessary information for calculating the
 * thermodynamic properties of an air parcel.
 *
 * This class has the following private members storing the properties:
 * - name, the name of this molecule, solids and liquids should add suffix "(s)"
 *   and "(l)" to the corresponding molecule.
 * - mu, molecular weight [g/mol]
 * - cp, solid or liquid heat capacity [J/(mol K)]
 * - latent, latent heat at triple point [kJ/mol]
 * - entropy, standard entropy [J/(mol K)]
 * - enthalpy, standard enthalpy [kJ/mol]
 * - gibbs, standard gibbs free energy [kJ/mol]
 * - tr, triple point temperature [K]
 * - pr, triple point pressure [pa]
 * - tc, critical point temperature [K]
 * - pc, critical point pressure [bar]
 * - cliq, liquid heat capacity [J/(mol K)]
 * - enliq, vaporization heat of liquids at triple point [kJ/mol]
 * - csld, solid heat capacity [J/(mol K)]
 * - ensld, sublimation heat of solids at triple point [kJ/mol]
 * - beta, slope of idealized saturation vapor pressure curve,
 *   beta = (L(T) - (cp - ci) * T) / (R * Tr)
 * - gamma, correction factor of saturation vapor pressure curve,
 *   gamma = (ci - cp) / R
 *
 * The following members store coefficients in thermodynamic expressions:
 * - shomate, 7 shomate expression coefficients
 * - shomate_sp, temperature separations in shomate expressions
 */

// number of Shomate coefficients
#define NSHOMATE 7

// maximum number of shomate expressions
#define MAXSHOMATE 3

class Molecule {
  friend std::ostream &operator<<(std::ostream &os, Molecule const &mol);

 protected:
  std::string m_name;

  double m_mu, m_cp, m_latent, m_entropy, m_enthalpy, m_gibbs, m_tr, m_pr, m_tc,
      m_pc, m_cliq, m_enliq, m_csld, m_ensld, m_beta, m_gamma;

  int m_nshomate;

  std::array<double, MAXSHOMATE * NSHOMATE> m_shomate;

  std::array<double, MAXSHOMATE + 1> m_shomate_sp;

 public:
  Molecule(std::string name = "")
      : m_name(name),
        m_mu(0),
        m_cp(0),
        m_latent(0),
        m_entropy(0),
        m_enthalpy(0),
        m_gibbs(0),
        m_tr(0),
        m_pr(0),
        m_tc(0),
        m_pc(0),
        m_cliq(0),
        m_enliq(0),
        m_csld(0),
        m_ensld(0),
        m_beta(0),
        m_gamma(0),
        m_nshomate(0) {
    std::fill(m_shomate.begin(), m_shomate.end(), 0.);
    std::fill(m_shomate_sp.begin(), m_shomate_sp.end(), 0.);
  }

  virtual ~Molecule() {}

  void LoadThermodynamicFile(std::string chemfile);

  virtual double cp(double T) const;

  virtual double enthalpy(double T) const;

  virtual double entropy(double T) const;

  virtual double latent(double T) const;

  double sat_vapor_p(double T) const {
    return m_pr * exp((1. - m_tr / T) * m_beta - m_gamma * log(T / m_tr));
  };

  double isat_vapor_p(double P) const;

  std::string name() const { return m_name; }

  double mu() const { return m_mu; }

  double cp() const { return m_cp; }

  double tr() const { return m_tr; }

  double pr() const { return m_pr; }

  double tc() const { return m_tc; }

  double pc() const { return m_pc; }

  double latent() const { return m_latent; }
};

class Hydrogen {
 public:
  Hydrogen();

  static double fpara_equil(double T);
  static double cp_para(double T);
  static double cp_ortho(double T);
  static double cp_norm(double T);
  static double cp_equil(double T);
  static double cp_nist(double T);

 private:
  static double FJ_[20];
  static double DJ_[20];
  static double nist_shomate1_[7];  // 298 ~ 1000 K
  static double nist_shomate2_[7];  // 1000 ~ 2500 K
  static double nist_shomate3_[7];  // 2500 ~ 6000 K
  static void get_cp_(double *cp_para, double *cp_equil, double *cp_norm,
                      double *cp_ortho, double T);
};

class Helium {
 public:
  static double cp_nist(double T);

 private:
  static double nist_shomate_[7];  // 300 ~ 600 K
};

class Methane {
 public:
  static double cp_nist(double T);

 private:
  static double nist_shomate1_[7];  // 298 ~ 1300 K
  static double nist_shomate2_[7];  // 1300 ~ 6000 K
};

#endif  // SRC_SNAP_THERMODYNAMICS_MOLECULES_HPP_
