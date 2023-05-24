// Athena++ headers
#include <athena/athena.hpp>

// Forward declarations
class MeshBlock;
class ParameterInput;

template <typename T>
class AthenaArray;

class FaceReconstruct {
 public:
  FaceReconstruct(MeshBlock *pmb, ParameterInput *pin);

  // weno reconstruction functions of various orders in each dimension
  void Weno3X1(const int k, const int j, const int il, const int iu,
               const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
               AthenaArray<Real> &wl, AthenaArray<Real> &wr);

  void Weno3X2(const int k, const int j, const int il, const int iu,
               const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
               AthenaArray<Real> &wl, AthenaArray<Real> &wr);

  void Weno3X3(const int k, const int j, const int il, const int iu,
               const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
               AthenaArray<Real> &wl, AthenaArray<Real> &wr);

  void Weno5X1(const int k, const int j, const int il, const int iu,
               const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
               AthenaArray<Real> &wl, AthenaArray<Real> &wr);

  void Weno5X2(const int k, const int j, const int il, const int iu,
               const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
               AthenaArray<Real> &wl, AthenaArray<Real> &wr);

  void Weno5X3(const int k, const int j, const int il, const int iu,
               const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
               AthenaArray<Real> &wl, AthenaArray<Real> &wr);

  void Weno5X1Simple(const int k, const int j, const int il, const int iu,
                     const AthenaArray<Real> &w, AthenaArray<Real> &wl,
                     AthenaArray<Real> &wr);

  void Weno5X2Simple(const int k, const int j, const int il, const int iu,
                     const AthenaArray<Real> &w, AthenaArray<Real> &wl,
                     AthenaArray<Real> &wr);

  void Weno5X3Simple(const int k, const int j, const int il, const int iu,
                     const AthenaArray<Real> &w, AthenaArray<Real> &wl,
                     AthenaArray<Real> &wr);

 protected:
  MeshBlock *pmy_block_;
};
