#ifndef USER_OUTPUTS_HPP
#define USER_OUTPUTS_HPP

// C/C++ headers
#include <vector>
#include <string>

// Athena++ headers
#include <outputs/outputs.hpp>

// canoe headers
#include "output_utils.hpp"

using DiagnosticTable = std::vector<std::vector<std::string>>;

//! \class DebugOutput
//  \brief derived OutputType class for debug dumps

class DebugOutput: public OutputType {
public:
  explicit DebugOutput(OutputParameters oparams) : OutputType(oparams) {} 
  void WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag) override;
};

//! \class NetcdfOutput
//  \brief derived OutputType class for Netcdf dumps

class NetcdfOutput : public OutputType {
public:
  NetcdfOutput(OutputParameters oparams);
  ~NetcdfOutput() {};
  void WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag);
  void CombineBlocks() override;
};

//! \class PnetcdfOutput
//  \brief derived OutputType class for parallel Netcdf dumps

class PnetcdfOutput : public OutputType {
public:
  PnetcdfOutput(OutputParameters oparams);
  ~PnetcdfOutput() {};
  void WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag);
};

//! \class FITSOutput
//  \brief derived OutputType class for FITS dumps

class FITSOutput : public OutputType {
public:
  FITSOutput(OutputParameters oparams);
  ~FITSOutput() {}
  void WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag);
};

#endif
