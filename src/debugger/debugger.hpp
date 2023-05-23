#ifndef DEBUGGER_HPP
#define DEBUGGER_HPP

// C/C++ headers
#include <cmath>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

// Athena++ header
// #include <athena.hpp>
// #include <globals.hpp>

// typedef int (*TestFunc_t)(Real);

// forward declaration
namespace Globals {
extern int my_rank, nranks;
}

class MaterialPoint;

// DEBUG_LEVEL = 0 : no debug output
//             = 1 : output important steps
//             = 2 : output conservation check
//             = 3 : output detailed particle transfer and bounds check
class Debugger {
 public:
  // data
  static std::string const cgreen;
  static std::string const cend;

  Debugger *prev, *next;
  std::stringstream msg;

  // functions
  Debugger(int depth = 0);
  ~Debugger();

  /*Debugger* StartTracking(std::string name);

  void Track3D(std::string name, TestFunc_t test, AthenaArray<Real>& var, int
  n);

  void Track1D(std::string name, TestFunc_t test, AthenaArray<Real>& var, int n,
  int k, int j); void DumpTracking(std::string name, int c1, int c2, int c3,
  char const* mode);*/
  // void Enter(char const *name);
  //
  void Enter(std::string name, std::string heil = "Initializing");

  void Call(std::string name) { Enter(name, "Calling"); }

  void Leave();

  /*void CheckConservation(std::string name, AthenaArray<Real> const& var,
      int is, int ie, int js, int je, int ks, int ke);

  void CheckParticleConservation(std::vector<std::string> const& cnames,
      std::vector<MaterialPoint> const& mp);*/

  Debugger* Message(std::string str);

  template <typename T>
  Debugger* Message(std::string str, T const& a);

  template <typename T>
  Debugger* Message(std::string str, T* a, int n);

  template <typename T>
  Debugger* Message(std::string str, std::vector<T>& a);

  static void Fatal(std::string where, std::string what);
  static void Fatal(std::string where, std::string str, std::string what);

  static void Print(std::string str);

  template <typename T>
  static void Print(std::string name, T const& value);

 protected:
  int depth_, current_depth_;
  std::string fname_;

  // AthenaArray<Real> data_;
  std::vector<std::string> vnames_;
  std::vector<std::string> sections_;
  std::vector<std::string> idstack_next_;
};

// global debugger
extern std::unique_ptr<Debugger> pdebug;

#include "debugger_impl.hpp"

void increment_id(std::string& str);

// small test functions
template <typename T>
int IsPositive(T v) {
  return v > 0 ? 1 : 0;
}

template <typename T>
int IsNumber(T v) {
  return !std::isnan(v);
}

#endif
