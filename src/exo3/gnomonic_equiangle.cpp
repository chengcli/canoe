// C/C++
#include <cmath>

// Athena++ headers
#include <athena/athena.hpp>
#include <athena/athena_arrays.hpp>
#include <athena/coordinates/coordinates.hpp>
#include <athena/defs.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// application
#include <application/application.hpp>

// exo3
#include "gnomonic_equiangle.hpp"

//----------------------------------------------------------------------------------------
//! Cartesian coordinates constructor

GnomonicEquiangle::GnomonicEquiangle(MeshBlock *pmb, ParameterInput *pin,
                                     bool flag)
    : Coordinates(pmb, pin, flag) {
  Application::Logger app("exo3");
  app->Log("Gnomonic Equiangle Coordinate");

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

  // initialize area-averaged coordinates used with MHD AMR
  if ((pmb->pmy_mesh->multilevel) && MAGNETIC_FIELDS_ENABLED) {
    for (int i = il - ng; i <= iu + ng; ++i) {
      x1s2(i) = x1s3(i) = x1v(i);
    }
    if (pmb->block_size.nx2 == 1) {
      x2s1(jl) = x2s3(jl) = x2v(jl);
    } else {
      for (int j = jl - ng; j <= ju + ng; ++j) {
        x2s1(j) = x2s3(j) = x2v(j);
      }
    }
    if (pmb->block_size.nx3 == 1) {
      x3s1(kl) = x3s2(kl) = x3v(kl);
    } else {
      for (int k = kl - ng; k <= ku + ng; ++k) {
        x3s1(k) = x3s2(k) = x3v(k);
      }
    }
  }

  // Initialize coordinate-transfer related variables
  g_.NewAthenaArray(NMETRIC, nc1 + 1);
  gi_.NewAthenaArray(NMETRIC, nc1 + 1);
}

void GnomonicEquiangle::Face1Area(const int k, const int j, const int il,
                                  const int iu, AthenaArray<Real> &area) {
#pragma omp simd
  Real xt = tan(x2v(j));
  Real yt = tan(x3v(k));
  Real sin_theta =
      sqrt((1.0 + xt * xt + yt * yt) / (1.0 + xt * xt) / (1.0 + yt * yt));

  Real x1 = tan(x2f(j));
  Real x2 = tan(x2f(j + 1));
  Real y = tan(x3v(k));
  Real delta1 = sqrt(1.0 + x1 * x1 + y * y);
  Real delta2 = sqrt(1.0 + x2 * x2 + y * y);
  Real dx2_ang = acos(1 / (delta1 * delta2) * (1 + x1 * x2 + y * y));

  Real x = tan(x2v(j));
  Real y1 = tan(x3f(k));
  Real y2 = tan(x3f(k + 1));
  delta1 = sqrt(1.0 + x * x + y1 * y1);
  delta2 = sqrt(1.0 + x * x + y2 * y2);
  Real dx3_ang = acos(1 / (delta1 * delta2) * (1 + x * x + y1 * y2));

  for (int i = il; i <= iu; ++i) {
    area(i) = dx2_ang * dx3_ang * x1f(i) * x1f(i) * sin_theta;
  }
  return;
}

void GnomonicEquiangle::Face2Area(const int k, const int j, const int il,
                                  const int iu, AthenaArray<Real> &area) {
  Real x = tan(x2f(j));
  Real y1 = tan(x3f(k));
  Real y2 = tan(x3f(k + 1));
  Real delta1 = sqrt(1.0 + x * x + y1 * y1);
  Real delta2 = sqrt(1.0 + x * x + y2 * y2);
#pragma omp simd
  for (int i = il; i <= iu; ++i) {
    Real dx3_lin = x1v(i) * acos(1 / (delta1 * delta2) * (1 + x * x + y1 * y2));
    Real &area_i = area(i);
    area_i = dx3_lin * dx1f(i);
  }
  return;
}

void GnomonicEquiangle::Face3Area(const int k, const int j, const int il,
                                  const int iu, AthenaArray<Real> &area) {
  Real x1 = tan(x2f(j));
  Real x2 = tan(x2f(j + 1));
  Real y = tan(x3f(k));
  Real delta1 = sqrt(1.0 + x1 * x1 + y * y);
  Real delta2 = sqrt(1.0 + x2 * x2 + y * y);
#pragma omp simd
  for (int i = il; i <= iu; ++i) {
    Real dx2_lin = x1v(i) * acos(1 / (delta1 * delta2) * (1 + x1 * x2 + y * y));
    Real &area_i = area(i);
    area_i = dx2_lin * dx1f(i);
  }
  return;
}

//----------------------------------------------------------------------------------------
// GetFaceXArea functions: return area of face with normal in X-dir at (i,j,k)

Real GnomonicEquiangle::GetFace1Area(const int k, const int j, const int i) {
  Real xt = tan(x2v(j));
  Real yt = tan(x3v(k));
  Real sin_theta =
      sqrt((1.0 + xt * xt + yt * yt) / (1.0 + xt * xt) / (1.0 + yt * yt));

  Real x1 = tan(x2f(j));
  Real x2 = tan(x2f(j + 1));
  Real y = tan(x3v(k));
  Real delta1 = sqrt(1.0 + x1 * x1 + y * y);
  Real delta2 = sqrt(1.0 + x2 * x2 + y * y);
  Real dx2_ang = acos(1 / (delta1 * delta2) * (1 + x1 * x2 + y * y));

  Real x = tan(x2v(j));
  Real y1 = tan(x3f(k));
  Real y2 = tan(x3f(k + 1));
  delta1 = sqrt(1.0 + x * x + y1 * y1);
  delta2 = sqrt(1.0 + x * x + y2 * y2);
  Real dx3_ang = acos(1 / (delta1 * delta2) * (1 + x * x + y1 * y2));

  return dx2_ang * dx3_ang * x1f(i) * x1f(i) * sin_theta;
}

Real GnomonicEquiangle::GetFace2Area(const int k, const int j, const int i) {
  Real x = tan(x2f(j));
  Real y1 = tan(x3f(k));
  Real y2 = tan(x3f(k + 1));
  Real delta1 = sqrt(1.0 + x * x + y1 * y1);
  Real delta2 = sqrt(1.0 + x * x + y2 * y2);
  Real dx3_lin = x1v(i) * acos(1 / (delta1 * delta2) * (1 + x * x + y1 * y2));
  return dx1f(i) * dx3_lin;
}

Real GnomonicEquiangle::GetFace3Area(const int k, const int j, const int i) {
  Real x1 = tan(x2f(j));
  Real x2 = tan(x2f(j + 1));
  Real y = tan(x3f(k));
  Real delta1 = sqrt(1.0 + x1 * x1 + y * y);
  Real delta2 = sqrt(1.0 + x2 * x2 + y * y);
  Real dx2_lin = x1v(i) * acos(1 / (delta1 * delta2) * (1 + x1 * x2 + y * y));
  return dx1f(i) * dx2_lin;
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
// VolCenterFaceXArea functions: compute area of face with normal in X-dir as
// vector where the faces are joined by cell centers (for non-ideal MHD)

void GnomonicEquiangle::VolCenterFace1Area(const int k, const int j,
                                           const int il, const int iu,
                                           AthenaArray<Real> &area) {
  Real xt = tan(x2v(j));
  Real yt = tan(x3v(k));
  Real sin_theta =
      sqrt((1.0 + xt * xt + yt * yt) / (1.0 + xt * xt) / (1.0 + yt * yt));

  Real x1 = tan(x2f(j));
  Real x2 = tan(x2f(j + 1));
  Real y = tan(x3v(k));
  Real delta1 = sqrt(1.0 + x1 * x1 + y * y);
  Real delta2 = sqrt(1.0 + x2 * x2 + y * y);
  Real dx2_ang = acos(1 / (delta1 * delta2) * (1 + x1 * x2 + y * y));

  Real x = tan(x2v(j));
  Real y1 = tan(x3f(k));
  Real y2 = tan(x3f(k + 1));
  delta1 = sqrt(1.0 + x * x + y1 * y1);
  delta2 = sqrt(1.0 + x * x + y2 * y2);
  Real dx3_ang = acos(1 / (delta1 * delta2) * (1 + x * x + y1 * y2));

#pragma omp simd
  for (int i = il; i <= iu; ++i) {
    Real &area_i = area(i);
    area_i = dx2_ang * dx3_ang * x1v(i) * x1v(i) * sin_theta;
  }
  return;
}

void GnomonicEquiangle::VolCenterFace2Area(const int k, const int j,
                                           const int il, const int iu,
                                           AthenaArray<Real> &area) {
  Real x = tan(x2v(j));
  Real y1 = tan(x3f(k));
  Real y2 = tan(x3f(k + 1));
  Real delta1 = sqrt(1.0 + x * x + y1 * y1);
  Real delta2 = sqrt(1.0 + x * x + y2 * y2);
#pragma omp simd
  for (int i = il; i <= iu; ++i) {
    Real dx3_lin = x1v(i) * acos(1 / (delta1 * delta2) * (1 + x * x + y1 * y2));
    Real &area_i = area(i);
    area_i = dx1v(i) * dx3_lin;
  }
  return;
}

void GnomonicEquiangle::VolCenterFace3Area(const int k, const int j,
                                           const int il, const int iu,
                                           AthenaArray<Real> &area) {
  Real x1 = tan(x2f(j));
  Real x2 = tan(x2f(j + 1));
  Real y = tan(x3v(k));
  Real delta1 = sqrt(1.0 + x1 * x1 + y * y);
  Real delta2 = sqrt(1.0 + x2 * x2 + y * y);
#pragma omp simd
  for (int i = il; i <= iu; ++i) {
    Real dx2_lin = x1v(i) * acos(1 / (delta1 * delta2) * (1 + x1 * x2 + y * y));
    Real &area_i = area(i);
    area_i = dx1v(i) * dx2_lin;
  }
  return;
}

// Cell Volume function: compute volume of cell as vector

void GnomonicEquiangle::CellVolume(const int k, const int j, const int il,
                                   const int iu, AthenaArray<Real> &vol) {
  Real xt = tan(x2v(j));
  Real yt = tan(x3v(k));
  Real sin_theta =
      sqrt((1.0 + xt * xt + yt * yt) / (1.0 + xt * xt) / (1.0 + yt * yt));

  Real x1 = tan(x2f(j));
  Real x2 = tan(x2f(j + 1));
  Real y = tan(x3v(k));
  Real delta1 = sqrt(1.0 + x1 * x1 + y * y);
  Real delta2 = sqrt(1.0 + x2 * x2 + y * y);
  Real dx2_ang = acos(1 / (delta1 * delta2) * (1 + x1 * x2 + y * y));

  Real x = tan(x2v(j));
  Real y1 = tan(x3f(k));
  Real y2 = tan(x3f(k + 1));
  delta1 = sqrt(1.0 + x * x + y1 * y1);
  delta2 = sqrt(1.0 + x * x + y2 * y2);
  Real dx3_ang = acos(1 / (delta1 * delta2) * (1 + x * x + y1 * y2));

#pragma omp simd
  for (int i = il; i <= iu; ++i) {
    vol(i) = dx1f(i) * dx2_ang * dx3_ang * x1v(i) * x1v(i) * sin_theta;
  }
  return;
}

//----------------------------------------------------------------------------------------
// CenterWidthX functions: return physical width in X-dir at (i,j,k) cell-center

void GnomonicEquiangle::CenterWidth1(const int k, const int j, const int il,
                                     const int iu, AthenaArray<Real> &dx1) {
#pragma omp simd
  for (int i = il; i <= iu; ++i) {
    dx1(i) = dx1f(i);  // Same as cartesian
  }
  return;
}

void GnomonicEquiangle::CenterWidth2(const int k, const int j, const int il,
                                     const int iu, AthenaArray<Real> &dx2) {
  Real x = tan(x2v(j));
  Real y1 = tan(x3f(k));
  Real y2 = tan(x3f(k + 1));
  Real delta1 = sqrt(1.0 + x * x + y1 * y1);
  Real delta2 = sqrt(1.0 + x * x + y2 * y2);
#pragma omp simd
  for (int i = il; i <= iu; ++i) {
    Real r = x1v(i);
    dx2(i) = r * acos(1 / (delta1 * delta2) * (1 + x * x + y1 * y2));
  }
  return;
}

void GnomonicEquiangle::CenterWidth3(const int k, const int j, const int il,
                                     const int iu, AthenaArray<Real> &dx3) {
  Real x1 = tan(x2f(j));
  Real x2 = tan(x2f(j + 1));
  Real y = tan(x3v(k));
  Real delta1 = sqrt(1.0 + x1 * x1 + y * y);
  Real delta2 = sqrt(1.0 + x2 * x2 + y * y);
#pragma omp simd
  for (int i = il; i <= iu; ++i) {
    Real r = x1v(i);
    dx3(i) = r * acos(1 / (delta1 * delta2) * (1 + x1 * x2 + y * y));
  }
  return;
}

//----------------------------------------------------------------------------------------
// GetCellVolume: returns cell volume at (i,j,k)

Real GnomonicEquiangle::GetCellVolume(const int k, const int j, const int i) {
  Real xt = tan(x2v(j));
  Real yt = tan(x3v(k));
  Real sin_theta =
      sqrt(1.0 + xt * xt + yt * yt / (1.0 + xt * xt) / (1.0 + yt * yt));
  return dx1f(i) * dx2f(j) * dx3f(k) * sin_theta;
}

//----------------------------------------------------------------------------------------
// For Affine coordinate: the metrics are all the same for volume, faces
void GnomonicEquiangle::CellMetric(const int k, const int j, const int il,
                                   const int iu, AthenaArray<Real> &g,
                                   AthenaArray<Real> &g_inv) {
  Real xt = tan(x2v(j));
  Real yt = tan(x3v(k));
  Real cos_theta = xt * yt / sqrt((1.0 + xt * xt) * (1.0 + yt * yt));
  Real sin_theta =
      sqrt((1.0 + xt * xt + yt * yt) / (1.0 + xt * xt) / (1.0 + yt * yt));
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
    g23 = cos_theta;
    g33 = 1.0;

    gi11 = 1.0;
    gi22 = 1.0 / sin_theta / sin_theta;
    gi33 = 1.0 / sin_theta / sin_theta;
  }
  return;
}

void GnomonicEquiangle::Face1Metric(const int k, const int j, const int il,
                                    const int iu, AthenaArray<Real> &g,
                                    AthenaArray<Real> &g_inv) {
  Real xt = tan(x2v(j));
  Real yt = tan(x3v(k));
  Real cos_theta = -xt * yt / sqrt((1.0 + xt * xt) * (1.0 + yt * yt));
  Real sin_theta =
      sqrt((1.0 + xt * xt + yt * yt) / (1.0 + xt * xt) / (1.0 + yt * yt));
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

    // Set metric terms
    g11 = 1.0;
    g22 = 1.0;
    g12 = 0.0;
    g13 = 0.0;
    g23 = cos_theta;
    g33 = 1.0;

    gi11 = 1.0;
    gi22 = 1.0 / sin_theta / sin_theta;
    gi33 = 1.0 / sin_theta / sin_theta;
  }
  return;
}

void GnomonicEquiangle::Face2Metric(const int k, const int j, const int il,
                                    const int iu, AthenaArray<Real> &g,
                                    AthenaArray<Real> &g_inv) {
  // Go through 1D block of cells
  Real xt = tan(x2f(j));
  Real yt = tan(x3v(k));
  Real cos_theta = -xt * yt / sqrt((1.0 + xt * xt) * (1.0 + yt * yt));
  Real sin_theta =
      sqrt((1.0 + xt * xt + yt * yt) / (1.0 + xt * xt) / (1.0 + yt * yt));
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

    // Set metric terms
    g11 = 1.0;
    g22 = 1.0;
    g12 = 0.0;
    g13 = 0.0;
    g23 = cos_theta;
    g33 = 1.0;

    gi11 = 1.0;
    gi22 = 1.0 / sin_theta / sin_theta;
    gi33 = 1.0 / sin_theta / sin_theta;
  }
  return;
}

void GnomonicEquiangle::Face3Metric(const int k, const int j, const int il,
                                    const int iu, AthenaArray<Real> &g,
                                    AthenaArray<Real> &g_inv) {
  // Go through 1D block of cells
  Real xt = tan(x2v(j));
  Real yt = tan(x3f(k));
  Real cos_theta = -xt * yt / sqrt((1.0 + xt * xt) * (1.0 + yt * yt));
  Real sin_theta =
      sqrt((1.0 + xt * xt + yt * yt) / (1.0 + xt * xt) / (1.0 + yt * yt));
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

    // Set metric terms
    g11 = 1.0;
    g22 = 1.0;
    g12 = 0.0;
    g13 = 0.0;
    g23 = cos_theta;
    g33 = 1.0;

    gi11 = 1.0;
    gi22 = 1.0 / sin_theta / sin_theta;
    gi33 = 1.0 / sin_theta / sin_theta;
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

void GnomonicEquiangle::PrimToLocal2(const int k, const int j, const int il,
                                     const int iu,
                                     const AthenaArray<Real> &b1_vals,
                                     AthenaArray<Real> &prim_left,
                                     AthenaArray<Real> &prim_right,
                                     AthenaArray<Real> &bx) {
  // Calculate metric coefficients for projection
  Face2Metric(k, j, il, iu, g_, gi_);
  // Go through 1D block of cells
#pragma omp simd
  for (int i = il; i <= iu; ++i) {
    Real &g11 = g_(I11, i);
    Real &g22 = g_(I22, i);
    Real &g33 = g_(I33, i);

    Real &gi11 = gi_(I11, i);
    Real &gi22 = gi_(I22, i);
    Real &gi33 = gi_(I33, i);

    Real &g12 = g_(I12, i);
    Real &g13 = g_(I13, i);
    Real &g23 = g_(I23, i);

    // Extract global projected 4-velocities
    Real uu1_l = prim_left(IVX, i);
    Real uu2_l = prim_left(IVY, i);
    Real uu3_l = prim_left(IVZ, i);
    Real uu1_r = prim_right(IVX, i);
    Real uu2_r = prim_right(IVY, i);
    Real uu3_r = prim_right(IVZ, i);

    // // Calculate transformation matrix
    Real T11 = 1.0;
    Real T12 = 0.0;
    Real T13 = 0.0;
    Real T21 = 0.0;
    Real T22 = 1.0 / sqrt(gi22);
    Real T23 = 0.0;
    Real T31 = 0.0;
    Real T32 = g23 / sqrt(g33);
    Real T33 = sqrt(g33);

    // // Transform projected velocities
    Real ux_l = T11 * uu1_l + T12 * uu2_l + T13 * uu3_l;
    Real uy_l = T21 * uu1_l + T22 * uu2_l + T23 * uu3_l;
    Real uz_l = T31 * uu1_l + T32 * uu2_l + T33 * uu3_l;
    Real ux_r = T11 * uu1_r + T12 * uu2_r + T13 * uu3_r;
    Real uy_r = T21 * uu1_r + T22 * uu2_r + T23 * uu3_r;
    Real uz_r = T31 * uu1_r + T32 * uu2_r + T33 * uu3_r;

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

void GnomonicEquiangle::PrimToLocal3(const int k, const int j, const int il,
                                     const int iu,
                                     const AthenaArray<Real> &b1_vals,
                                     AthenaArray<Real> &prim_left,
                                     AthenaArray<Real> &prim_right,
                                     AthenaArray<Real> &bx) {
  // Calculate metric coefficients for projection
  Face3Metric(k, j, il, iu, g_, gi_);

  // Go through 1D block of cells
#pragma omp simd
  for (int i = il; i <= iu; ++i) {
    Real &g11 = g_(I11, i);
    Real &g22 = g_(I22, i);
    Real &g33 = g_(I33, i);

    Real &gi11 = gi_(I11, i);
    Real &gi22 = gi_(I22, i);
    Real &gi33 = gi_(I33, i);

    Real &g12 = g_(I12, i);
    Real &g13 = g_(I13, i);
    Real &g23 = g_(I23, i);

    // Extract global projected 4-velocities
    Real uu1_l = prim_left(IVX, i);
    Real uu2_l = prim_left(IVY, i);
    Real uu3_l = prim_left(IVZ, i);
    Real uu1_r = prim_right(IVX, i);
    Real uu2_r = prim_right(IVY, i);
    Real uu3_r = prim_right(IVZ, i);

    // // Calculate transformation matrix
    Real T11 = 1.0;
    Real T12 = 0.0;
    Real T13 = 0.0;
    Real T21 = 0.0;
    Real T22 = sqrt(g22);
    Real T23 = g23 / sqrt(g22);
    Real T31 = 0.0;
    Real T32 = 0.0;
    Real T33 = 1 / sqrt(gi33);

    // // Transform projected velocities
    Real ux_l = T11 * uu1_l + T12 * uu2_l + T13 * uu3_l;
    Real uy_l = T21 * uu1_l + T22 * uu2_l + T23 * uu3_l;
    Real uz_l = T31 * uu1_l + T32 * uu2_l + T33 * uu3_l;
    Real ux_r = T11 * uu1_r + T12 * uu2_r + T13 * uu3_r;
    Real uy_r = T21 * uu1_r + T22 * uu2_r + T23 * uu3_r;
    Real uz_r = T31 * uu1_r + T32 * uu2_r + T33 * uu3_r;

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

void GnomonicEquiangle::FluxToGlobal2(
    const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bbx,
    AthenaArray<Real> &flux, AthenaArray<Real> &ey, AthenaArray<Real> &ez) {
  // Extract metrics
  Face2Metric(k, j, il, iu, g_, gi_);

  // Go through 1D block of cells
#pragma omp simd
  for (int i = il; i <= iu; ++i) {
    // Extract metric tensors
    Real &g11 = g_(I11, i);
    Real &g22 = g_(I22, i);
    Real &g33 = g_(I33, i);

    Real &gi11 = gi_(I11, i);
    Real &gi22 = gi_(I22, i);
    Real &gi33 = gi_(I33, i);

    Real &g12 = g_(I12, i);
    Real &g13 = g_(I13, i);
    Real &g23 = g_(I23, i);

    // Extract local conserved quantities and fluxes
    const Real txx = flux(IM1, k, j, i);
    const Real txy = flux(IM2, k, j, i);
    const Real txz = flux(IM3, k, j, i);

    // Calculate transformation matrix
    Real T11 = 1.0;
    Real T12 = 0.0;
    Real T13 = 0.0;
    Real T21 = 0.0;
    Real T22 = sqrt(gi22);
    Real T23 = 0.0;
    Real T31 = 0.0;
    Real T32 = -sqrt(gi22) * g23 / g33;
    Real T33 = 1.0 / sqrt(g33);

    // Extract global fluxes
    Real &t1_1 = flux(IM1, k, j, i);
    Real &t1_2 = flux(IM2, k, j, i);
    Real &t1_3 = flux(IM3, k, j, i);

    // Set fluxes
    t1_1 = T11 * txx + T12 * txy + T13 * txz;
    t1_2 = T21 * txx + T22 * txy + T23 * txz;
    t1_3 = T31 * txx + T32 * txy + T33 * txz;
  }
  return;
}

void GnomonicEquiangle::FluxToGlobal3(
    const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bbx,
    AthenaArray<Real> &flux, AthenaArray<Real> &ey, AthenaArray<Real> &ez) {
  // Extract metrics
  Face3Metric(k, j, il, iu, g_, gi_);

  // Go through 1D block of cells
#pragma omp simd
  for (int i = il; i <= iu; ++i) {
    // Extract metric tensors
    Real &g11 = g_(I11, i);
    Real &g22 = g_(I22, i);
    Real &g33 = g_(I33, i);

    Real &gi11 = gi_(I11, i);
    Real &gi22 = gi_(I22, i);
    Real &gi33 = gi_(I33, i);

    Real &g12 = g_(I12, i);
    Real &g13 = g_(I13, i);
    Real &g23 = g_(I23, i);

    // Extract local conserved quantities and fluxes
    const Real txx = flux(IM1, k, j, i);
    const Real txy = flux(IM2, k, j, i);
    const Real txz = flux(IM3, k, j, i);

    // Calculate transformation matrix
    Real T11 = 1.0;
    Real T12 = 0.0;
    Real T13 = 0.0;
    Real T21 = 0.0;
    Real T22 = 1.0 / sqrt(g22);
    Real T23 = -g23 / g22 * sqrt(gi33);
    Real T31 = 0.0;
    Real T32 = 0.0;
    Real T33 = sqrt(gi33);

    // Extract global fluxes
    Real &t1_1 = flux(IM1, k, j, i);
    Real &t1_2 = flux(IM2, k, j, i);
    Real &t1_3 = flux(IM3, k, j, i);

    // Set fluxes
    t1_1 = T11 * txx + T12 * txy + T13 * txz;
    t1_2 = T21 * txx + T22 * txy + T23 * txz;
    t1_3 = T31 * txx + T32 * txy + T33 * txz;
  }
  return;
}

void GnomonicEquiangle::AddCoordTermsDivergence(const Real dt,
                                                const AthenaArray<Real> *flux,
                                                const AthenaArray<Real> &prim,
                                                const AthenaArray<Real> &bcc,
                                                AthenaArray<Real> &u) {
  for (int k = pmy_block->ks; k <= pmy_block->ke; ++k) {
    for (int j = pmy_block->js; j <= pmy_block->je; ++j) {
#pragma omp simd
      for (int i = pmy_block->is; i <= pmy_block->ie; ++i) {
        // General variables
        Real v1 = prim(IVX, k, j, i);
        Real v2 = prim(IVY, k, j, i);
        Real v3 = prim(IVZ, k, j, i);
        Real r = x1v(i);
        Real x = tan(x2v(j));
        Real y = tan(x3v(k));
        Real C = sqrt(1.0 + x * x);
        Real D = sqrt(1.0 + y * y);
        Real delta = 1.0 / (1.0 + x * x + y * y);
        Real pr;
        Real rho;
        if (strcmp(EQUATION_OF_STATE, "shallow_yz") == 0) {
          pr = 0.5 * prim(IDN, k, j, i) * prim(IDN, k, j, i);
          rho = prim(IDN, k, j, i);
        } else {
          pr = prim(IPR, k, j, i);
          rho = prim(IDN, k, j, i);
          // Update flux 1 (excluded from shallow water case)
          Real src1 = 2.0 * pr / r +
                      rho * (v2 * v2 + v3 * v3 - 2 * v2 * v3 * x * y / (C * D));
          u(IM1, k, j, i) += dt * src1;
        }
        // Update flux 2
        Real src2 = pr * y * y / r * x / D -
                    rho * v2 / r * (v1 - y * v3 * delta * delta / (C * D * D));
        u(IM2, k, j, i) += dt * src2;

        // Update flux 3
        Real src3 = pr * x * x / r * y / C -
                    rho * v3 / r * (v1 - x * v2 * delta * delta / (C * C * D));
        u(IM3, k, j, i) += dt * src3;
      }
    }
  }
  return;
}
