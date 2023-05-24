/** @file test_thermo.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Tuesday Jul 13, 2021 14:20:17 PDT
 * @bug No known bugs.
 */

#include <iostream>

#include "../molecules.hpp"
using namespace std;

int main() {
  double T = 1200.;
  std::cout << Hydrogen::fpara_equil(T) << std::endl;
  std::cout << Hydrogen::cp_para(T) << std::endl;
  std::cout << Hydrogen::cp_ortho(T) << std::endl;
  std::cout << Hydrogen::cp_norm(T) << std::endl;
  std::cout << Hydrogen::cp_equil(T) << std::endl;
  std::cout << Hydrogen::cp_nist(T) << std::endl;
  std::cout << std::endl;
  std::cout << Helium::cp_nist(100.) << std::endl;
  std::cout << Helium::cp_nist(300.) << std::endl;
  std::cout << Helium::cp_nist(500.) << std::endl;
  std::cout << Helium::cp_nist(1000.) << std::endl;
  std::cout << std::endl;
  std::cout << Methane::cp_nist(300.) << std::endl;
  std::cout << Methane::cp_nist(1299.) << std::endl;
  std::cout << Methane::cp_nist(1301.) << std::endl;
}
