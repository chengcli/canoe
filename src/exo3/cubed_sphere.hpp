#ifndef SRC_EXO3_CUBED_SPHERE_HPP_
#define SRC_EXO3_CUBED_SPHERE_HPP_

// C/C++
#include <memory>

// athena
#include <athena/athena.hpp>
#include <athena/mesh/mesh.hpp>

// canoe
#include <configure.hpp>

class CubedSphere {
 public:
  CubedSphere(MeshBlock *pmb);
  ~CubedSphere() {}

  static void TransformOX(int *ox2, int *ox3, int *tox2, int *tox3,
                          LogicalLocation const &loc);

  void GetLatLon(Real *lat, Real *lon, int k, int j, int i) const;
  void GetLatLonFace2(Real *lat, Real *lon, int k, int j, int i) const;
  void GetLatLonFace3(Real *lat, Real *lon, int k, int j, int i) const;

  void GetUV(Real *U, Real *V, Real V2, Real V3, int k, int j, int i) const;
  void GetVyVz(Real *V2, Real *V3, Real U, Real V, int k, int j, int i) const;

  void CalculateCoriolisForce2(int i2, int i3, Real v2, Real v3, Real Omega,
                               Real den, Real *cF2, Real *cF3) const;
  void CalculateCoriolisForce3(int i2, int i3, Real v1, Real v2, Real v3,
                               Real Omega, Real den, Real *cF1, Real *cF2,
                               Real *cF3) const;
  void ContravariantVectorToCovariant(int i2, int i3, Real v2, Real v3,
                                      Real *v2c, Real *v3c) const;
  void CovariantVectorToContravariant(int i2, int i3, Real v2, Real v3,
                                      Real *v2c, Real *v3c) const;

  void TransformOX(int *ox2, int *ox3, int *tox2, int *tox3) const {
    return TransformOX(ox2, ox3, tox2, tox3, pmy_block_->loc);
  }

  void SaveLR3DValues(AthenaArray<Real> &L_in, AthenaArray<Real> &R_in,
                      int direction, int k, int j, int il, int iu);
  void LoadLR3DValues(AthenaArray<Real> &L_in, AthenaArray<Real> &R_in,
                      int direction, int k, int j, int il, int iu);

  void SynchronizeFluxesSend();
  void SynchronizeFluxesRecv();
  void SynchronizeFluxesWait();

 protected:
  void sendNeighborBlocks(int ox2, int ox3, int tg_rank, int tg_gid);
  void recvNeighborBlocks(int ox2, int ox3, int tg_rank, int tg_gid);

  AthenaArray<Real> L3DValues[3], R3DValues[3];
  std::unordered_map<Connection, Connection> LRDataBuffer[4];

#ifdef MPI_PARALLEL
  int send_flag_[4], recv_flag_[4];
  MPI_Request send_request_[4];
  MPI_Request recv_request_[4];
#endif

  std::vector<Real> LRDataBuffer[4];

  MeshBlock *pmy_block_;
};

using CubedSpherePtr = std::shared_ptr<CubedSphere>;

#endif  // SRC_EXO3_CUBED_SPHERE_HPP_
