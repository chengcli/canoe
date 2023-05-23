// C/C++ headers
#include <iostream>

// debugger headers
#include "debugger.hpp"

template <typename T>
Debugger* Debugger::Message(std::string str, T const& a) {
  msg << "- " << str << " = " << a << std::endl;
  return this;
}

template <typename T>
Debugger* Debugger::Message(std::string str, T* a, int n) {
  msg << "- " << str << " = ";
  for (int i = 0; i < n; ++i) msg << a[i] << " ";
  msg << std::endl;
  return this;
}

template <typename T>
Debugger* Debugger::Message(std::string str, std::vector<T>& a) {
  msg << "- " << str << " = ";
  for (size_t i = 0; i < a.size(); ++i) msg << a[i] << " ";
  msg << std::endl;
  return this;
}

template <typename T>
void Debugger::Print(std::string name, T const& value) {
  if (Globals::my_rank == 0) std::cout << name << " = " << value << std::endl;
}
