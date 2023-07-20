// Athena++ headers
#include <athena/athena.hpp>
#include <athena/athena_arrays.hpp>
#include <athena/coordinates/coordinates.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// exo3
#include "affine_coordinate.hpp"

//----------------------------------------------------------------------------------------
//! Cartesian coordinates constructor

AffineCoordinate::AffineCoordinate(MeshBlock *pmb, ParameterInput *pin,
                                   bool flag)
    : Coordinates(pmb, pin, flag) {
  // Send something to confirm that we are using Affine
  std::cout << "===Note===: Affine coordinates activated" << std::endl;
  theta_ = PI / 3.;
  // theta_ = PI/2.;
  sin_theta_ = sin(theta_);
  cos_theta_ = cos(theta_);

  // initialize volume-averaged coordinates and spacing
  // x1-direction: x1v = dx/2
  for (int i = il - ng; i <= iu + ng; ++i) {
    x1v(i) = 0.5 * (x1f(i + 1) + x1f(i));
  }
  for (int i = il - ng; i <= iu + ng - 1; ++i) {
    if (pmb->block_size.x1rat != 1.0) {
      dx1v(i) = x1v(i + 1) - x1v(i);
    } else {
      // dx1v = dx1f constant for uniform mesh; may disagree with x1v(i+1) -
      // x1v(i)
      dx1v(i) = dx1f(i);
    }
  }

  // x2-direction: x2v = dy/2
  if (pmb->block_size.nx2 == 1) {
    x2v(jl) = 0.5 * (x2f(jl + 1) + x2f(jl));
    dx2v(jl) = dx2f(jl);
  } else {
    for (int j = jl - ng; j <= ju + ng; ++j) {
      x2v(j) = 0.5 * (x2f(j + 1) + x2f(j));
    }
    for (int j = jl - ng; j <= ju + ng - 1; ++j) {
      if (pmb->block_size.x2rat != 1.0) {
        dx2v(j) = x2v(j + 1) - x2v(j);
      } else {
        // dx2v = dx2f constant for uniform mesh; may disagree with x2v(j+1) -
        // x2v(j)
        dx2v(j) = dx2f(j);
      }
    }
  }

  // x3-direction: x3v = dz/2
  if (pmb->block_size.nx3 == 1) {
    x3v(kl) = 0.5 * (x3f(kl + 1) + x3f(kl));
    dx3v(kl) = dx3f(kl);
  } else {
    for (int k = kl - ng; k <= ku + ng; ++k) {
      x3v(k) = 0.5 * (x3f(k + 1) + x3f(k));
    }
    for (int k = kl - ng; k <= ku + ng - 1; ++k) {
      if (pmb->block_size.x3rat != 1.0) {
        dx3v(k) = x3v(k + 1) - x3v(k);
      } else {
        // dxkv = dx3f constant for uniform mesh; may disagree with x3v(k+1) -
        // x3v(k)
        dx3v(k) = dx3f(k);
      }
    }
  }
  // initialize geometry coefficients
  // x1-direction
  for (int i = il - ng; i <= iu + ng; ++i) {
    h2v(i) = 1.0;
    h2f(i) = 1.0;
    h31v(i) = 1.0;
    h31f(i) = 1.0;
    dh2vd1(i) = 0.0;
    dh2fd1(i) = 0.0;
    dh31vd1(i) = 0.0;
    dh31fd1(i) = 0.0;
  }

  // x2-direction
  if (pmb->block_size.nx2 == 1) {
    h32v(jl) = 1.0;
    h32f(jl) = 1.0;
    dh32vd2(jl) = 0.0;
    dh32fd2(jl) = 0.0;
  } else {
    for (int j = jl - ng; j <= ju + ng; ++j) {
      h32v(j) = 1.0;
      h32f(j) = 1.0;
      dh32vd2(j) = 0.0;
      dh32fd2(j) = 0.0;
    }
  }

  // Initialize coordinate-transfer related variables
  g_.NewAthenaArray(NMETRIC, nc1 + 1);
  gi_.NewAthenaArray(NMETRIC, nc1 + 1);

  // Not needed here, but precalculation of the parts of metrics
  // may be possible in gnomonic equiangle...
}

// Put in the changes in face2area etc, similar to cylindrical.cpp
// In affine coordinates, the differences lie in face1area and face2area.
// face3area cancels the sqrt(g) term exactly.

//----------------------------------------------------------------------------------------
// FaceXArea functions: compute area of face with normal in X-dir as vector

void AffineCoordinate::Face1Area(const int k, const int j, const int il,
                                 const int iu, AthenaArray<Real> &area) {
#pragma omp simd
  for (int i = il; i <= iu; ++i) {
    Real &area_i = area(i);
    area_i = sin_theta_ * dx2f(j) * dx3f(k);
  }
  return;
}

void AffineCoordinate::Face2Area(const int k, const int j, const int il,
                                 const int iu, AthenaArray<Real> &area) {
#pragma omp simd
  for (int i = il; i <= iu; ++i) {
    Real &area_i = area(i);
    area_i = dx1f(i) * dx3f(k);
  }
  return;
}

void AffineCoordinate::Face3Area(const int k, const int j, const int il,
                                 const int iu, AthenaArray<Real> &area) {
#pragma omp simd
  for (int i = il; i <= iu; ++i) {
    Real &area_i = area(i);
    area_i = dx1f(i) * dx2f(j);
  }
  return;
}

//----------------------------------------------------------------------------------------
// GetFaceXArea functions: return area of face with normal in X-dir at (i,j,k)

Real AffineCoordinate::GetFace1Area(const int k, const int j, const int i) {
  return dx2f(j) * dx3f(k) * sin_theta_;
}

Real AffineCoordinate::GetFace2Area(const int k, const int j, const int i) {
  return dx1f(i) * dx3f(k);
}

Real AffineCoordinate::GetFace3Area(const int k, const int j, const int i) {
  return dx1f(i) * dx2f(j);
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
// VolCenterFaceXArea functions: compute area of face with normal in X-dir as
// vector where the faces are joined by cell centers (for non-ideal MHD)

void AffineCoordinate::VolCenterFace1Area(const int k, const int j,
                                          const int il, const int iu,
                                          AthenaArray<Real> &area) {
#pragma omp simd
  for (int i = il; i <= iu; ++i) {
    Real &area_i = area(i);
    area_i = dx2v(j) * dx3v(k) * sin_theta_;
  }
  return;
}

void AffineCoordinate::VolCenterFace2Area(const int k, const int j,
                                          const int il, const int iu,
                                          AthenaArray<Real> &area) {
#pragma omp simd
  for (int i = il; i <= iu; ++i) {
    Real &area_i = area(i);
    area_i = dx1v(i) * dx3v(k);
  }
  return;
}

void AffineCoordinate::VolCenterFace3Area(const int k, const int j,
                                          const int il, const int iu,
                                          AthenaArray<Real> &area) {
#pragma omp simd
  for (int i = il; i <= iu; ++i) {
    Real &area_i = area(i);
    area_i = dx1v(i) * dx2v(j);
  }
  return;
}

// Cell Volume function: compute volume of cell as vector

void AffineCoordinate::CellVolume(const int k, const int j, const int il,
                                  const int iu, AthenaArray<Real> &vol) {
#pragma omp simd
  for (int i = il; i <= iu; ++i) {
    vol(i) = dx1f(i) * dx2f(j) * dx3f(k) * sin_theta_;
  }
  return;
}

//----------------------------------------------------------------------------------------
// GetCellVolume: returns cell volume at (i,j,k)

Real AffineCoordinate::GetCellVolume(const int k, const int j, const int i) {
  return dx1f(i) * dx2f(j) * dx3f(k) * sin_theta_;
}

//----------------------------------------------------------------------------------------
// For Affine coordinate: the metrics are all the same for volume, faces
void AffineCoordinate::CellMetric(const int k, const int j, const int il,
                                  const int iu, AthenaArray<Real> &g,
                                  AthenaArray<Real> &g_inv) {
  // Go through 1D block of cells
#pragma omp simd
  for (int i = il; i <= iu; ++i) {
    // Extract metric terms
    Real &g11 = g(I11, i);
    Real &g22 = g(I22, i);
    Real &g33 = g(I33, i);

    Real &gi11 = g_inv(I11, i);
    Real &gi22 = g_inv(I22, i);
    Real &gi33 = g_inv(I33, i);

    Real &g12 = g(I12, i);
    Real &g13 = g(I13, i);
    Real &g23 = g(I23, i);

    // Set metric terms, we only use covariant g for all calculations
    g11 = 1.0;
    g22 = 1.0;
    g12 = 0.0;
    g13 = 0.0;
    g23 = cos_theta_;
    g33 = 1.0;

    gi11 = 1.0;
    gi22 = 1.0 / (sin_theta_ * sin_theta_);
    gi33 = 1.0 / (sin_theta_ * sin_theta_);
  }
  return;
}

void AffineCoordinate::Face1Metric(const int k, const int j, const int il,
                                   const int iu, AthenaArray<Real> &g,
                                   AthenaArray<Real> &g_inv) {
  // Go through 1D block of cells
#pragma omp simd
  for (int i = il; i <= iu; ++i) {
    // Extract remaining geometric quantities

    // Extract metric terms
    Real &g11 = g(I11, i);
    Real &g22 = g(I22, i);
    Real &g33 = g(I33, i);

    Real &gi11 = g_inv(I11, i);
    Real &gi22 = g_inv(I22, i);
    Real &gi33 = g_inv(I33, i);

    Real &g12 = g(I12, i);
    Real &g13 = g(I13, i);
    Real &g23 = g(I23, i);

    // Set metric terms, we only use covariant g for all calculations
    g11 = 1.0;
    g22 = 1.0;
    g12 = 0.0;
    g13 = 0.0;
    g23 = cos_theta_;
    g33 = 1.0;

    gi11 = 1.0;
    gi22 = 1.0 / (sin_theta_ * sin_theta_);
    gi33 = 1.0 / (sin_theta_ * sin_theta_);
  }
  return;
}

void AffineCoordinate::Face2Metric(const int k, const int j, const int il,
                                   const int iu, AthenaArray<Real> &g,
                                   AthenaArray<Real> &g_inv) {
  // Go through 1D block of cells
#pragma omp simd
  for (int i = il; i <= iu; ++i) {
    // Extract remaining geometric quantities

    // Extract metric terms
    Real &g11 = g(I11, i);
    Real &g22 = g(I22, i);
    Real &g33 = g(I33, i);

    Real &gi11 = g_inv(I11, i);
    Real &gi22 = g_inv(I22, i);
    Real &gi33 = g_inv(I33, i);

    Real &g12 = g(I12, i);
    Real &g13 = g(I13, i);
    Real &g23 = g(I23, i);

    // Set metric terms, we only use covariant g for all calculations
    g11 = 1.0;
    g22 = 1.0;
    g12 = 0.0;
    g13 = 0.0;
    g23 = cos_theta_;
    g33 = 1.0;

    gi11 = 1.0;
    gi22 = 1.0 / (sin_theta_ * sin_theta_);
    gi33 = 1.0 / (sin_theta_ * sin_theta_);
  }
  return;
}

void AffineCoordinate::Face3Metric(const int k, const int j, const int il,
                                   const int iu, AthenaArray<Real> &g,
                                   AthenaArray<Real> &g_inv) {
  // Go through 1D block of cells
#pragma omp simd
  for (int i = il; i <= iu; ++i) {
    // Extract remaining geometric quantities

    // Extract metric terms
    Real &g11 = g(I11, i);
    Real &g22 = g(I22, i);
    Real &g33 = g(I33, i);

    Real &gi11 = g_inv(I11, i);
    Real &gi22 = g_inv(I22, i);
    Real &gi33 = g_inv(I33, i);

    Real &g12 = g(I12, i);
    Real &g13 = g(I13, i);
    Real &g23 = g(I23, i);

    // Set metric terms, we only use covariant g for all calculations
    g11 = 1.0;
    g22 = 1.0;
    g12 = 0.0;
    g13 = 0.0;
    g23 = cos_theta_;
    g33 = 1.0;

    gi11 = 1.0;
    gi22 = 1.0 / (sin_theta_ * sin_theta_);
    gi33 = 1.0 / (sin_theta_ * sin_theta_);
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for transforming primitives to locally flat frame: r-interface
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
//   prim_l: 1D array of left primitives, using global coordinates
//   prim_r: 1D array of right primitives, using global coordinates
// Outputs:
//   prim_l: values overwritten in local coordinates
//   prim_r: values overwritten in local coordinates
// Notes: no notes

void AffineCoordinate::PrimToLocal2(const int k, const int j, const int il,
                                    const int iu,
                                    const AthenaArray<Real> &b1_vals,
                                    AthenaArray<Real> &prim_left,
                                    AthenaArray<Real> &prim_right,
                                    AthenaArray<Real> &bx) {
  // Calculate metric coefficients for projection
  // Face2Metric(k, j, il, iu, g_, gi_);

  // Go through 1D block of cells
#pragma omp simd
  for (int i = il; i <= iu; ++i) {
    // Real &g11 = g_(I11,i);
    // Real &g22 = g_(I22,i);
    // Real &g33 = g_(I33,i);

    // Real &gi11 = gi_(I11,i);
    // Real &gi22 = gi_(I22,i);
    // Real &gi33 = gi_(I33,i);

    // Real &g12 = g_(I12,i);
    // Real &g13 = g_(I13,i);
    // Real &g23 = g_(I23,i);

    // Extract global projected 4-velocities
    Real uu1_l = prim_left(IVX, i);
    Real uu2_l = prim_left(IVY, i);
    Real uu3_l = prim_left(IVZ, i);
    Real uu1_r = prim_right(IVX, i);
    Real uu2_r = prim_right(IVY, i);
    Real uu3_r = prim_right(IVZ, i);

    // // Calculate transformation matrix
    // Real T11 = sqrt(g11);
    // Real T12 = g12/sqrt(g11);
    // Real T13 = g13/sqrt(g11);
    // Real T21 = 0.0;
    // Real T22 = 1.0/sqrt(gi22);
    // Real T23 = 0.0;
    // Real T31 = g13/sqrt(g33);
    // Real T32 = g23/sqrt(g33);
    // Real T33 = sqrt(g33);

    // // Transform projected 4-velocities
    // // Differ from Schwartzchild here only...
    // Real ux_l = T11*uu1_l+T12*uu2_l+T13*uu3_l;
    // Real uy_l = T21*uu1_l+T22*uu2_l+T23*uu3_l;
    // Real uz_l = T31*uu1_l+T32*uu2_l+T33*uu3_l;
    // Real ux_r = T11*uu1_r+T12*uu2_r+T13*uu3_r;
    // Real uy_r = T11*uu1_r+T12*uu2_r+T13*uu3_r;
    // Real uz_r = T11*uu1_r+T12*uu2_r+T13*uu3_r;

    // affine: quick
    Real ux_l = uu1_l;
    Real uy_l = uu2_l * sin_theta_;
    Real uz_l = uu2_l * cos_theta_ + uu3_l;
    Real ux_r = uu1_r;
    Real uy_r = uu2_r * sin_theta_;
    Real uz_r = uu2_r * cos_theta_ + uu3_r;

    // Set local projected 4-velocities
    prim_left(IVX, i) = ux_l;
    prim_left(IVY, i) = uy_l;
    prim_left(IVZ, i) = uz_l;
    prim_right(IVX, i) = ux_r;
    prim_right(IVY, i) = uy_r;
    prim_right(IVZ, i) = uz_r;
  }
  return;
}

void AffineCoordinate::PrimToLocal3(const int k, const int j, const int il,
                                    const int iu,
                                    const AthenaArray<Real> &b1_vals,
                                    AthenaArray<Real> &prim_left,
                                    AthenaArray<Real> &prim_right,
                                    AthenaArray<Real> &bx) {
  // Calculate metric coefficients for projection
  // Face3Metric(k, j, il, iu, g_, gi_);

  // Go through 1D block of cells
#pragma omp simd
  for (int i = il; i <= iu; ++i) {
    // Real &g11 = g_(I11,i);
    // Real &g22 = g_(I22,i);
    // Real &g33 = g_(I33,i);

    // Real &gi11 = gi_(I11,i);
    // Real &gi22 = gi_(I22,i);
    // Real &gi33 = gi_(I33,i);

    // Real &g12 = g_(I12,i);
    // Real &g13 = g_(I13,i);
    // Real &g23 = g_(I23,i);

    // Extract global projected 4-velocities
    Real uu1_l = prim_left(IVX, i);
    Real uu2_l = prim_left(IVY, i);
    Real uu3_l = prim_left(IVZ, i);
    Real uu1_r = prim_right(IVX, i);
    Real uu2_r = prim_right(IVY, i);
    Real uu3_r = prim_right(IVZ, i);

    // Calculate transformation matrix
    // Real T11 = sqrt(g11);
    // Real T12 = g12/sqrt(g11);
    // Real T13 = g13/sqrt(g11);
    // Real T21 = g12/sqrt(g22);
    // Real T22 = sqrt(g22);
    // Real T23 = g23/sqrt(g22);
    // Real T31 = 0.0;
    // Real T32 = 0.0;
    // Real T33 = 1.0/sqrt(gi33);

    // Transform projected 4-velocities
    // Differ from Schwartzchild here only...
    // Real ux_l = T11*uu1_l+T12*uu2_l+T13*uu3_l;
    // Real uy_l = T21*uu1_l+T22*uu2_l+T23*uu3_l;
    // Real uz_l = T31*uu1_l+T32*uu2_l+T33*uu3_l;
    // Real ux_r = T11*uu1_r+T12*uu2_r+T13*uu3_r;
    // Real uy_r = T11*uu1_r+T12*uu2_r+T13*uu3_r;
    // Real uz_r = T11*uu1_r+T12*uu2_r+T13*uu3_r;

    // affine: simple form
    Real ux_l = uu1_l;
    Real uy_l = uu2_l + uu3_l * cos_theta_;
    Real uz_l = uu3_l * sin_theta_;
    Real ux_r = uu1_r;
    Real uy_r = uu2_r + uu3_r * cos_theta_;
    Real uz_r = uu3_r * sin_theta_;

    // Set local projected 4-velocities
    prim_left(IVX, i) = ux_l;
    prim_left(IVY, i) = uy_l;
    prim_left(IVZ, i) = uz_l;
    prim_right(IVX, i) = ux_r;
    prim_right(IVY, i) = uy_r;
    prim_right(IVZ, i) = uz_r;
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for transforming fluxes to global frame: r-interface
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
//   flux: 3D array of hydrodynamical fluxes, using local coordinates
//   ey,ez: 3D arrays of magnetic fluxes (electric fields), using local
//   coordinates
// Outputs:
//   flux: values overwritten in global coordinates
//   ey,ez: values overwritten in global coordinates
// Notes:
//   expects values and x-fluxes of Mx/My/Mz in IM1/IM2/IM3 slots
//   puts r-fluxes of M1/M2/M3 in IM1/IM2/IM3 slots

void AffineCoordinate::FluxToGlobal2(
    const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bbx,
    AthenaArray<Real> &flux, AthenaArray<Real> &ey, AthenaArray<Real> &ez) {
  // Extract metrics
  Face2Metric(k, j, il, iu, g_, gi_);

  // Go through 1D block of cells
#pragma omp simd
  for (int i = il; i <= iu; ++i) {
    // Extract metric tensors
    // Real &g11 = g_(I11,i);
    // Real &g22 = g_(I22,i);
    // Real &g33 = g_(I33,i);

    // Real &gi11 = gi_(I11,i);
    // Real &gi22 = gi_(I22,i);
    // Real &gi33 = gi_(I33,i);

    // Real &g12 = g_(I12,i);
    // Real &g13 = g_(I13,i);
    // Real &g23 = g_(I23,i);

    // Extract local conserved quantities and fluxes
    const Real txx = flux(IM1, k, j, i);
    const Real txy = flux(IM2, k, j, i);
    const Real txz = flux(IM3, k, j, i);

    // Calculate transformation matrix
    // const Real D = g11*g33-g13*g13;
    // Real T11 = g33*sqrt(g11)/D;
    // Real T12 = sqrt(gi22)*(g13*g23-g12*g33)/D;
    // Real T13 = -g13*sqrt(g33)/D;
    // Real T21 = 0.0;
    // Real T22 = sqrt(gi22);
    // Real T23 = 0.0;
    // Real T31 = -g13*sqrt(g11)/D;
    // Real T32 = sqrt(gi22)*(g12*g13-g23*g11)/D;
    // Real T33 = g11*sqrt(g33)/D;

    // Extract global fluxes
    Real &t1_1 = flux(IM1, k, j, i);
    Real &t1_2 = flux(IM2, k, j, i);
    Real &t1_3 = flux(IM3, k, j, i);

    // Set fluxes
    t1_1 = txx;
    t1_2 = txy / sin_theta_;
    t1_3 = -txy * cos_theta_ / sin_theta_ + txz;
  }
  return;
}

void AffineCoordinate::FluxToGlobal3(
    const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bbx,
    AthenaArray<Real> &flux, AthenaArray<Real> &ey, AthenaArray<Real> &ez) {
  // Extract metrics
  Face3Metric(k, j, il, iu, g_, gi_);

  // Go through 1D block of cells
#pragma omp simd
  for (int i = il; i <= iu; ++i) {
    // Extract metric tensors
    // Real &g11 = g_(I11,i);
    // Real &g22 = g_(I22,i);
    // Real &g33 = g_(I33,i);

    // Real &gi11 = gi_(I11,i);
    // Real &gi22 = gi_(I22,i);
    // Real &gi33 = gi_(I33,i);

    // Real &g12 = g_(I12,i);
    // Real &g13 = g_(I13,i);
    // Real &g23 = g_(I23,i);

    // Extract local conserved quantities and fluxes
    const Real txx = flux(IM1, k, j, i);
    const Real txy = flux(IM2, k, j, i);
    const Real txz = flux(IM3, k, j, i);

    // Calculate transformation matrix
    // const Real D = g11*g22-g12*g12;
    // Real T11 = g22*sqrt(g11)/D;
    // Real T12 = -g12*sqrt(g22)/D;
    // Real T13 = sqrt(gi33)*(g12*g23-g13*g22)/D;
    // Real T21 = -g12*sqrt(g11)/D;
    // Real T22 = g11*sqrt(g22)/D;
    // Real T23 = sqrt(gi33)*(g12*g13-g23*g11)/D;
    // Real T31 = 0.0;
    // Real T32 = 0.0;
    // Real T33 = sqrt(gi33);

    // Extract global fluxes
    Real &t1_1 = flux(IM1, k, j, i);
    Real &t1_2 = flux(IM2, k, j, i);
    Real &t1_3 = flux(IM3, k, j, i);

    // Set fluxes
    t1_1 = txx;
    t1_2 = txy - txz * cos_theta_ / sin_theta_;
    t1_3 = txz / sin_theta_;
  }
  return;
}
