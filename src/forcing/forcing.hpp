#ifndef SRC_FORCING_FORCING_HPP_
#define SRC_FORCING_FORCING_HPP_

// C/C++ headers
#include <string>
#include <vector>

// Athena++ headers
#include <athena/athena.hpp>

// canoe
#include <virtual_groups.hpp>

// Forward declarations
class MeshBlock;
class ParameterInput;

//! \brief manages all physics package data and functions
class Forcing : public ParameterGroup {
 public:
  // functions
  Forcing() {}
  virtual ~Forcing() {}

  virtual void Initialize(MeshBlock *pmb) {}

  virtual size_t RestartDataSizeInBytes() const { return 0; }
  virtual size_t DumpRestartData(char *pdst) const { return 0; }
  virtual size_t LoadRestartData(char *psrc) { return 0; }

  virtual void Apply(AthenaArray<Real> &du, MeshBlock *pmb, Real time,
                     Real dt) = 0;
};

class BotForcing : public Forcing {
 public:
  BotForcing(MeshBlock *pmb, int nvar);

  virtual size_t RestartDataSizeInBytes() const override;
  virtual size_t DumpRestartData(char *pdst) const override;
  virtual size_t LoadRestartData(char *psrc) override;

 protected:
  AthenaArray<Real> bot_data_;
};

class RelaxBotTemp : public BotForcing {
 public:
  RelaxBotTemp(MeshBlock *pmb, ParameterInput *pin);

  void Initialize(MeshBlock *pmb) override;
  void Apply(AthenaArray<Real> &du, MeshBlock *pmb, Real time,
             Real dt) override;
};

class RelaxBotVelo : public BotForcing {
 public:
  RelaxBotVelo(MeshBlock *pmb, ParameterInput *pin);

  void Initialize(MeshBlock *pmb) override;
  void Apply(AthenaArray<Real> &du, MeshBlock *pmb, Real time,
             Real dt) override;
};

class RelaxBotComp : public BotForcing {
 public:
  RelaxBotComp(MeshBlock *pmb, ParameterInput *pin);

  void Initialize(MeshBlock *pmb) override;
  void Apply(AthenaArray<Real> &du, MeshBlock *pmb, Real time,
             Real dt) override;
};

class TopSpongeLyr : public Forcing {
 public:
  TopSpongeLyr(MeshBlock *pmb, ParameterInput *pin);

  void Apply(AthenaArray<Real> &du, MeshBlock *pmb, Real time,
             Real dt) override;
};

class BotSpongeLyr : public Forcing {
 public:
  BotSpongeLyr(MeshBlock *pmb, ParameterInput *pin);

  void Apply(AthenaArray<Real> &du, MeshBlock *pmb, Real time,
             Real dt) override;
};

class LeftSpongeLyr : public Forcing {
 public:
  LeftSpongeLyr(MeshBlock *pmb, ParameterInput *pin);

  void Apply(AthenaArray<Real> &du, MeshBlock *pmb, Real time,
             Real dt) override;
};

class RightSpongeLyr : public Forcing {
 public:
  RightSpongeLyr(MeshBlock *pmb, ParameterInput *pin);

  void Apply(AthenaArray<Real> &du, MeshBlock *pmb, Real time,
             Real dt) override;
};

class TopCooling : public Forcing {
 public:
  TopCooling(MeshBlock *pmb, ParameterInput *pin);

  void Apply(AthenaArray<Real> &du, MeshBlock *pmb, Real time,
             Real dt) override;
};

class BotHeating : public Forcing {
 public:
  BotHeating(MeshBlock *pmb, ParameterInput *pin);

  void Apply(AthenaArray<Real> &du, MeshBlock *pmb, Real time,
             Real dt) override;
};

class BodyHeating : public Forcing {
 public:
  BodyHeating(MeshBlock *pmb, ParameterInput *pin);

  void Apply(AthenaArray<Real> &du, MeshBlock *pmb, Real time,
             Real dt) override;
};

using ForcingPtr = std::shared_ptr<Forcing>;
using ForcingContainer = std::vector<ForcingPtr>;

class ForcingFactory {
 public:
  static ForcingContainer CreateFrom(MeshBlock *pmb, ParameterInput *pin);
};

using ForcingFunc = void (Forcing::*)(AthenaArray<Real> &, MeshBlock *, Real,
                                      Real);

#endif  // SRC_FORCING_FORCING_HPP_
