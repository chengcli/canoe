#ifndef SRC_OUTPUTS_USER_OUTPUTS_HPP_
#define SRC_OUTPUTS_USER_OUTPUTS_HPP_

// C/C++
#include <string>
#include <vector>

// athena
#include <athena/outputs/outputs.hpp>

// outputs
#include "output_utils.hpp"

using DiagnosticTable = std::vector<std::vector<std::string>>;

//! \class DebugOutput
//  \brief derived OutputType class for debug dumps

class DebugOutput : public OutputType {
 public:
  explicit DebugOutput(OutputParameters oparams) : OutputType(oparams) {}
  void WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag) override;
};

//! \class NetcdfOutput
//  \brief derived OutputType class for Netcdf dumps

class NetcdfOutput : public OutputType {
 public:
  explicit NetcdfOutput(OutputParameters oparams);
  ~NetcdfOutput() {}
  void WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag) override;
  void CombineBlocks() override;
};

//! \class PnetcdfOutput
//  \brief derived OutputType class for parallel Netcdf dumps

class PnetcdfOutput : public OutputType {
 public:
  explicit PnetcdfOutput(OutputParameters oparams);
  ~PnetcdfOutput() {}
  void WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag) override;
};

//! \class FITSOutput
//  \brief derived OutputType class for FITS dumps

class FITSOutput : public OutputType {
 public:
  explicit FITSOutput(OutputParameters oparams);
  ~FITSOutput() {}
  void WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag) override;
};

#endif  // SRC_OUTPUTS_USER_OUTPUTS_HPP_
