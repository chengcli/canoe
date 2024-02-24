#ifndef SRC_SNAP_IMPLICIT_IMPLICIT_SOLVER_HPP_
#define SRC_SNAP_IMPLICIT_IMPLICIT_SOLVER_HPP_

// C/C++
#include <memory>
#include <string>
#include <vector>

// Eigen
#include <Eigen/Core>

// athena
#include <athena/athena.hpp>
// #include <bvals/bvals_interfaces.hpp>

class MeshBlock;
class ParameterInput;
class BlockIndex;
class Coordinates;
class Mesh;
class BoundaryValues;
class EquationOfState;
class Thermodynamics;

template <typename T>
class AthenaArray;

class ImplicitSolver {
  friend class Hydro;

 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  // data
  bool has_top_neighbor, has_bot_neighbor;
  bool first_block, last_block;
  bool periodic_boundary;
  bool pole_at_bot, pole_at_top;
  NeighborBlock tblock, bblock;

  // functions
  ImplicitSolver(MeshBlock *pmb, ParameterInput *pin);
  ~ImplicitSolver();
  void SetDirection(CoordinateDirection dir);
  int GetImplicitFlag() const { return implicit_flag_; }
  void SolveImplicit3D(AthenaArray<Real> &du, AthenaArray<Real> &w, Real dt);

  // utility functions
  void FindNeighbors();
  int CreateMPITag(int lid, int bufid, std::string phy);

  // correction methods
  void PartialCorrection(AthenaArray<Real> &du, AthenaArray<Real> const &w,
                         Real dt, int k, int j, int is, int ie);
  void FullCorrection(AthenaArray<Real> &du, AthenaArray<Real> const &w,
                      Real dt, int k, int j, int is, int ie);

  // tri-diagonal solver
  template <typename T1, typename T2>
  void ForwardSweep(std::vector<T1> &a, std::vector<T1> &b, std::vector<T1> &c,
                    std::vector<T2> &delta, std::vector<T2> &corr, Real dt,
                    int k, int j, int il, int iu);

  template <typename T1, typename T2>
  void BackwardSubstitution(std::vector<T1> &a, std::vector<T2> &delta, int k,
                            int j, int il, int iu);

  // periodic solver
  template <typename T1, typename T2>
  void PeriodicForwardSweep(std::vector<T1> &a, std::vector<T1> &b,
                            std::vector<T1> &c, std::vector<T2> &corr, Real dt,
                            int k, int j, int il, int iu);

  template <typename T1, typename T2>
  void PeriodicBackwardSubstitution(std::vector<T1> &a, std::vector<T1> &c,
                                    std::vector<T2> &delta, int k, int j,
                                    int il, int iu);

  // forcing jacobians
  template <typename T>
  void JacobianGravityCoriolis(T &jac, Real const prim[], int k, int j, int i);

  // communications
  /*void SynchronizeConserved(AthenaArray<Real> const& du,
    int kl, int ku, int jl, int ju, int is, int ie);
  void WaitToFinishSync(int kl, int ku, int jl, int ju, int is, int ie);*/

  template <typename T>
  void SendBuffer(T const &a, int k, int j, NeighborBlock nb);

  template <typename T1, typename T2>
  void SendBuffer(T1 const &a, T2 const &b, int k, int j, NeighborBlock nb);

  template <typename T1, typename T2, typename T3, typename T4, typename T5,
            typename T6>
  void SendBuffer(T1 const &a, T2 const &b, T3 const &c, T4 const &d,
                  T5 const &e, T6 const &f, int k, int j, NeighborBlock ntop);

  template <typename T1, typename T2, typename T3, typename T4, typename T5,
            typename T6, typename T7>
  void SendBuffer(T1 const &a, T2 const &b, T3 const &c, T4 const &d,
                  T5 const &e, T6 const &f, T7 const &g, int k, int j,
                  NeighborBlock ntop);

  template <typename T>
  void RecvBuffer(T &a, int k, int j, NeighborBlock nb);

  template <typename T1, typename T2>
  void RecvBuffer(T1 &a, T2 &b, int k, int j, NeighborBlock nb);

  template <typename T1, typename T2, typename T3, typename T4, typename T5,
            typename T6>
  void RecvBuffer(T1 &a, T2 &b, T3 &c, T4 &d, T5 &e, T6 &f, int k, int j,
                  NeighborBlock nb);

  template <typename T1, typename T2, typename T3, typename T4, typename T5,
            typename T6, typename T7>
  void RecvBuffer(T1 &a, T2 &b, T3 &c, T4 &d, T5 &e, T6 &f, T7 &g, int k, int j,
                  NeighborBlock nb);

  template <typename T1, typename T2>
  void SaveCoefficients(std::vector<T1> &a, std::vector<T2> &b, int k, int j,
                        int il, int iu);

  template <typename T1, typename T2, typename T3, typename T4>
  void SaveCoefficients(std::vector<T1> &a, std::vector<T2> &b,
                        std::vector<T3> &c, std::vector<T4> &d, int k, int j,
                        int il, int iu);

  template <typename T1, typename T2>
  void LoadCoefficients(std::vector<T1> &a, std::vector<T2> &b, int k, int j,
                        int il, int iu);

  template <typename T1, typename T2, typename T3, typename T4>
  void LoadCoefficients(std::vector<T1> &a, std::vector<T2> &b,
                        std::vector<T3> &c, std::vector<T4> &d, int k, int j,
                        int il, int iu);

  // template<typename T>
  // void LoadForcingJacobian(T &phi, int k, int j ,int i, CoordinateDirection
  // dir);

 private:
  int implicit_flag_;
  MeshBlock *pmy_block_;

  CoordinateDirection mydir_;
  // Real *usend_top_, *urecv_top_;      // MPI data buffer
  // Real *usend_bot_, *urecv_bot_;      // MPI data buffer

  Real ***buffer_;  // MPI data buffer
  // Real ****jacobian_;     // archive of forcing jacobian
  AthenaArray<Real> du_;  // stores implicit solution
  AthenaArray<Real>
      coefficients_;  // archive of coefficients in the tri-diagonal matrix

  Eigen::Matrix<Real, 5, 5> p2_, p3_;  // perturbation matrices

#ifdef MPI_PARALLEL
  MPI_Request **req_send_data1_;
  MPI_Request **req_send_data2_;
  MPI_Request **req_send_data6_;
  MPI_Request **req_send_data7_;
  // MPI_Request req_send_sync_top_;
  // MPI_Request req_send_sync_bot_;
  // MPI_Request req_recv_sync_top_;
  // MPI_Request req_recv_sync_bot_;
#endif
};

using ImplicitSolverPtr = std::shared_ptr<ImplicitSolver>;

#endif  //  SRC_SNAP_IMPLICIT_IMPLICIT_SOLVER_HPP_
