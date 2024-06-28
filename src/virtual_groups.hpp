#ifndef SRC_VIRTUAL_GROUPS_HPP_
#define SRC_VIRTUAL_GROUPS_HPP_

// C/C++
#include <map>
#include <string>
#include <unordered_map>

// athena
#include <athena/athena.hpp>  // Real

// canoe
#include <configure.hpp>

class Mesh;
class MeshBlock;
class OutputType;
class OutputParameters;

namespace YAML {
class Node;
}

class NamedGroup {
 public:
  explicit NamedGroup(std::string name) : myname_(name) {}
  NamedGroup(std::string name, std::string long_name)
      : myname_(name), long_name_(long_name) {}

  virtual ~NamedGroup() {}

  std::string GetName() const { return myname_; }
  std::string GetLongName() const { return long_name_; }
  void SetName(std::string name) { myname_ = name; }
  void SetLongName(std::string long_name) { long_name_ = long_name; }

 private:
  std::string myname_;
  std::string long_name_;
};

class FlagGroup {
 public:
  explicit FlagGroup(uint64_t flags = 0LL) : myflags_(flags) {}
  virtual ~FlagGroup() {}

  int TestFlag(uint64_t flag) const { return myflags_ & flag; }
  void SetFlag(uint64_t flag) { myflags_ |= flag; }

 private:
  //! internal flags
  uint64_t myflags_;
};

class StringReprGroup {
 public:
  virtual ~StringReprGroup() {}
  virtual std::string ToString() const = 0;
};

class ParameterGroup {
 public:
  virtual ~ParameterGroup() {}

  void SetRealsFrom(YAML::Node const &node);

  void SetIntsFrom(YAML::Node const &node);

  //! Set real parameter
  void SetPar(std::string const &name, Real value) {
    params_real_[name] = value;
  }

  //! Set int parameter
  void SetPar(std::string const &name, int value) { params_int_[name] = value; }

  //! Set string parameter
  void SetPar(std::string const &name, std::string const &value) {
    params_str_[name] = value;
  }

  //! Get parameter
  template <typename T>
  T GetPar(std::string const &name) const;

  //! Check if a parameter exists
  bool HasPar(std::string const &name) const {
    if (params_real_.count(name) > 0) return true;
    if (params_int_.count(name) > 0) return true;
    if (params_str_.count(name) > 0) return true;
    return false;
  }

 private:
  //! real parameters
  std::unordered_map<std::string, Real> params_real_;

  //! int parameters
  std::unordered_map<std::string, int> params_int_;

  //! string parameters
  std::unordered_map<std::string, std::string> params_str_;
};

// specialization
template <>
inline int ParameterGroup::GetPar<int>(std::string const &name) const {
  return params_int_.at(name);
}

template <>
inline Real ParameterGroup::GetPar<Real>(std::string const &name) const {
  return params_real_.at(name);
}

template <>
inline std::string ParameterGroup::GetPar<std::string>(
    std::string const &name) const {
  return params_str_.at(name);
}

class CounterGroup {
 public:
  virtual ~CounterGroup() {}
  void ResetCounter() { current_ = cooldown_; }
  void DecrementCounter(Real dt) { current_ -= dt; }
  Real GetCounter() const { return current_; }
  void SetCooldownTime(Real cooldown) { cooldown_ = cooldown; }

 private:
  Real cooldown_, current_ = 0.;
};

class RestartGroup {
 public:
  virtual ~RestartGroup() {}
  virtual size_t RestartDataSizeInBytes(Mesh const *pm) const = 0;
  virtual void DumpRestartData(char *pdst) const = 0;
  virtual size_t LoadRestartData(char *psrt) = 0;
};

class ASCIIOutputGroup {
 public:
  virtual ~ASCIIOutputGroup() {}
  virtual void WriteAsciiHeader(OutputParameters const *) const = 0;
  virtual void WriteAsciiData(OutputParameters const *) const = 0;
};

class BinaryOutputGroup {
 public:
  virtual ~BinaryOutputGroup() {}
  virtual void WriteBinaryHeader(std::ofstream &os) const = 0;
  virtual void WriteBinaryData(std::ofstream &os) const = 0;
};

class MeshOutputGroup {
 public:
  virtual ~MeshOutputGroup() {}
  virtual bool ShouldMeshOutput(std::string variable_name) const = 0;
  virtual void LoadMeshOutputData(OutputType *out, int *num_vars) const = 0;
};

class FITSOutputGroup {
 public:
  virtual ~FITSOutputGroup() {}
  virtual bool ShouldFITSOutput(std::string variable_name) const = 0;
  virtual void LoadFITSOutputData(OutputType *out, int *num_vars) const = 0;
};

#endif  // SRC_VIRTUAL_GROUPS_HPP_
