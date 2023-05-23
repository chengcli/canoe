/** @file molecules.hpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Saturday Jul 10, 2021 14:50:00 PDT
 * @bug No known bugs.
 */

#ifndef MOLECULES_HPP
#define MOLECULES_HPP

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

#endif /* end of include guard HYDROGEN_HPP */
