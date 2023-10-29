#ifndef SRC_VIRTUAL_GROUPS_HPP_
#define SRC_VIRTUAL_GROUPS_HPP_

// C/C++
#include <string>

class OutputType;
class MeshBlock;

class NamedGroup {
 public:
  NamedGroup(std::string name) : myname_(name) {}
  virtual ~NamedGroup() {}
  virtual std::string GetName() const { return myname_; }

 protected:
  std::string myname_;
};

class RestartGroup {
 public:
  virtual ~RestartGroup() {}
  virtual size_t RestartDataSizeInBytes() const = 0;
  virtual size_t DumpRestartData(char *pdst) const = 0;
  virtual size_t LoadRestartData(char *psrt) = 0;
};

class ASCIIOutputGroup {
 public:
  virtual ~ASCIIOutputGroup() {}
  virtual void WriteAsciiHeader(std::ofstream &os) const = 0;
  virtual void WriteAsciiData(std::ofstream &os) const = 0;
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
