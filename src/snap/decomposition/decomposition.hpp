#ifndef SRC_SNAP_DECOMPOSITION_DECOMPOSITION_HPP_
#define SRC_SNAP_DECOMPOSITION_DECOMPOSITION_HPP_

// C/C++
#include <memory>

// athena
#include <athena/athena.hpp>
#include <athena/bvals/bvals.hpp>

// exchanger
#include <exchanger/exchanger.hpp>

// MPI headers
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

class MeshBlock;
template <typename T>
class AthenaArray;

class Decomposition {
 public:  // public access data
  // data
  bool has_top_neighbor, has_bot_neighbor;
  NeighborBlock tblock, bblock;

 public:
  // functions
  explicit Decomposition(MeshBlock *pmb);
  ~Decomposition();
  void FindNeighbors();
  int CreateMPITag(int lid, int bufid, int phy);

  void RecvFromTop(AthenaArray<Real> &psf, int kl, int ku, int jl, int ju);
  void SendToBottom(AthenaArray<Real> const &psf, int kl, int ku, int jl,
                    int ju);
  void SyncNewVariables(AthenaArray<Real> const &w, int kl, int ku, int jl,
                        int ju);
  void WaitToFinishSync(AthenaArray<Real> &w, int kl, int ku, int jl, int ju);

  void ApplyHydroBoundary(AthenaArray<Real> &w, AthenaArray<Real> &psf, int kl,
                          int ku, int jl, int ju);

  void RecvBuffer(AthenaArray<Real> &psf, int kl, int ku, int jl, int ju,
                  int il, int iu, NeighborBlock nb);
  void SendBuffer(AthenaArray<Real> const &psf, int kl, int ku, int jl, int ju);
  void PopulateBotEntropy(AthenaArray<Real> const &w, int kl, int ku, int jl,
                          int ju);
  void WaitToFinishSend();

  void ChangeToPerturbation(AthenaArray<Real> &w, int kl, int ku, int jl,
                            int ju);
  void RestoreFromPerturbation(AthenaArray<Real> &w, AthenaArray<Real> &wl,
                               AthenaArray<Real> &wr, int k, int j, int il,
                               int iu);

  void ChangeToBuoyancy(AthenaArray<Real> &w, int kl, int ku, int jl, int ju);
  void RestoreFromBuoyancy(AthenaArray<Real> &w, AthenaArray<Real> &wl,
                           AthenaArray<Real> &wr, int k, int j, int il, int iu);

  void ChangeToEntropy(AthenaArray<Real> &w, int kl, int ku, int jl, int ju);
  void RestoreFromEntropy(AthenaArray<Real> &w, AthenaArray<Real> &wl,
                          AthenaArray<Real> &wr, int k, int j, int il, int iu);

  void ChangeToAnomaly(AthenaArray<Real> &w, int kl, int ku, int jl, int ju);
  void RestoreFromAnomaly(AthenaArray<Real> &w, AthenaArray<Real> &wl,
                          AthenaArray<Real> &wr, int k, int j, int il, int iu);

  void ChangeToTemperature(AthenaArray<Real> &w, int kl, int ku, int jl,
                           int ju);
  void RestoreFromTemperature(AthenaArray<Real> &w, AthenaArray<Real> &wl,
                              AthenaArray<Real> &wr, int k, int j, int il,
                              int iu);

 private:
  MeshBlock *pmy_block_;

  // pressure decomposition
  AthenaArray<Real> psf_, psv_;
  AthenaArray<Real> dsf_, dsv_;
  AthenaArray<Real> pres_, dens_;  // save of original w
                                   //
  Real *buffer_, *send_buffer_;    // MPI data buffer
  Real *wsend_top_, *wrecv_top_;   // MPI data buffer
  Real *wsend_bot_, *wrecv_bot_;   // MPI data buffer
  int *brank_, *color_;            // bottom block rank and color

  AthenaArray<Real> entropy_;  // adiabatic index and pseudo entropy
  AthenaArray<Real> gamma_;    // adiabatic index

#ifdef MPI_PARALLEL
  MPI_Request req_send_top_;
  MPI_Request req_send_bot_;
  MPI_Request req_send_sync_top_;
  MPI_Request req_send_sync_bot_;
  MPI_Request req_recv_sync_top_;
  MPI_Request req_recv_sync_bot_;
#endif
};

using DecompositionPtr = std::shared_ptr<Decomposition>;

#endif  // SRC_SNAP_DECOMPOSITION_DECOMPOSITION_HPP_
