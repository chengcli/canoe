#ifndef SRC_FORCING_FORCING_HPP_
#define SRC_FORCING_FORCING_HPP_

// C/C++ headers
#include <string>
#include <vector>

// Athena++ headers
#include <athena/athena_array.hpp>

// canoe
#include <virtual_groups.hpp>

// Forward declarations
class MeshBlock;
class ParameterInput;

//! \brief manages all physics package data and functions
class Forcing : public ParameterGroup {
 public:
  // functions
  Forcing(MeshBlock *pmb, ParameterInput *pin);
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
  virtual ~BotForcing() {}

  virtual size_t RestartDataSizeInBytes() const override;
  virtual size_t DumpRestartData(char *pdst) const override;
  virtual size_t LoadRestartData(char *psrc) override;

  void Apply(AthenaArray<Real> &du, MeshBlock *pmb, Real time,
             Real dt) override;

 protected:
  AthenaArray<Real> data_;
};

class RelaxBotTemp : public BotForcing {
 public:
  void Apply(AthenaArray<Real> &du, MeshBlock *pmb, Real time,
             Real dt) override;
};

class RelaxBotVelo : public BotForcing {
 public:
  void Apply(AthenaArray<Real> &du, MeshBlock *pmb, Real time,
             Real dt) override;
};

class RelaxBotComp : public BotForcing {
 public:
  void Apply(AthenaArray<Real> &du, MeshBlock *pmb, Real time,
             Real dt) override;
};

class TopSpongeLyr : public Forcing {
 public:
  void Apply(AthenaArray<Real> &du, MeshBlock *pmb, Real time,
             Real dt) override;
};

class BotSpongeLyr : public Forcing {
 public:
  void Apply(AthenaArray<Real> &du, MeshBlock *pmb, Real time,
             Real dt) override;
};

class LeftSpongeLyr : public Forcing {
 public:
  void Apply(AthenaArray<Real> &du, MeshBlock *pmb, Real time,
             Real dt) override;
};

class RightSpongeLyr : public Forcing {
 public:
  void Apply(AthenaArray<Real> &du, MeshBlock *pmb, Real time,
             Real dt) override;
};

class TopCooling : public Forcing {
 public:
  void Apply(AthenaArray<Real> &du, MeshBlock *pmb, Real time,
             Real dt) override;

 protected:
  Real flux_top_; /**< top heating flux [W/m^2] */
};

class BotHeating : public Forcing {
 public:
  void Apply(AthenaArray<Real> &du, MeshBlock *pmb, Real time,
             Real dt) override;

 protected:
  Real flux_bot_; /**< bot heating flux [W/m^2] */
};

class BodyHeating : public Forcing {
 public:
  void Apply(AthenaArray<Real> &du, MeshBlock *pmb, Real time,
             Real dt) override;
};

using ForcingPtr = std::shared_ptr<Forcing>;
using ForcingContainer = std::vector<ForcingPtr>;

class ForcingFactory {
 public:
  static ForcingContainer CreateFrom(MeshBlock *pmb, ParameterInput *pin);
};
}
;

using ForcingFunc = void (Forcing::*)(AthenaArray<Real> &, MeshBlock *, Real,
                                      Real);

#endif  // SRC_FORCING_FORCING_HPP_
