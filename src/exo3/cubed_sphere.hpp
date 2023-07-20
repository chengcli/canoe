#ifndef CUBED_SPHERE_HPP
#define CUBED_SPHERE_HPP
#include "athena.hpp"
#include "coordinates/coordinates.hpp"

int FindBlockID(LogicalLocation const &loc);

void TransformOxForCubedSphere(int *ox2, int *ox3, int *tox2, int *tox3,
                               LogicalLocation const &loc);

void PackDataCubedSphere(const AthenaArray<Real> &src, Real *buf, int sn,
                         int en, int si, int ei, int sj, int ej, int sk, int ek,
                         int &offset, int ox1, int ox2, int ox3,
                         LogicalLocation const &loc);

void ProjectLocalCartesianAffine(const AthenaArray<Real> &src,
                                 AthenaArray<Real> &tgt, Real affine_angle,
                                 int sn, int en, int si, int ei, int sj, int ej,
                                 int sk, int ek, int Dir);

void DeProjectLocalCartesianAffine(AthenaArray<Real> &flux, Real affine_angle,
                                   int sn, int en, int si, int ei, int sj,
                                   int ej, int sk, int ek, int Dir);

Real CubedSphereMeshGeneratorX2(Real x, LogicalLocation const &loc);
Real CubedSphereMeshGeneratorX3(Real x, LogicalLocation const &loc);

// Helper functions
void GetLatLon(Real *lat, Real *lon, Coordinates *pcoord, int k, int j, int i);
void GetLatLonFace2(Real *lat, Real *lon, Coordinates *pcoord, int k, int j,
                    int i);
void GetLatLonFace3(Real *lat, Real *lon, Coordinates *pcoord, int k, int j,
                    int i);

void GetUV(Real *U, Real *V, Coordinates *pcoord, Real V2, Real V3, int k,
           int j, int i);
void GetVyVz(Real *V2, Real *V3, Coordinates *pcoord, Real U, Real V, int k,
             int j, int i);
// Helper functions adapted from Paul
void VecTransABPFromRLL(Real X, Real Y, int blockID, Real U, Real V, Real *V2,
                        Real *V3);
void VecTransRLLFromABP(Real X, Real Y, int blockID, Real V2, Real V3, Real *U,
                        Real *V);
void RLLFromXYP(Real dX, Real dY, int nP, Real &lon, Real &lat);
void XYPFromRLL(Real lon, Real lat, Real &dX, Real &dY, int &nP);

class CubedSphereLR {
 public:
  void SetMeshBlock(MeshBlock *pmb_in);
  void InitializeSizes(int nc3, int nc2, int nc1);
  void SaveLR3DValues(AthenaArray<Real> &L_in, AthenaArray<Real> &R_in,
                      int direction, int k, int j, int il, int iu);
  void LoadLR3DValues(AthenaArray<Real> &L_in, AthenaArray<Real> &R_in,
                      int direction, int k, int j, int il, int iu);
  void SynchronizeFluxes();
  void SendNeighborBlocks(LogicalLocation const &loc, int ox2, int ox3,
                          int tg_rank, int tg_gid);
  void RecvNeighborBlocks(LogicalLocation const &loc, int ox2, int ox3,
                          int tg_rank, int tg_gid);

  AthenaArray<Real> L3DValues, R3DValues;
  MeshBlock *pmb;
};

#endif
