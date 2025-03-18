// Athena++ headers
#include <athena/athena.hpp>
#include <athena/athena_arrays.hpp>
#include <athena/coordinates/coordinates.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// exo3
#include "warped_coordinate.hpp"

//----------------------------------------------------------------------------------------
//! Warped coordinates constructor

WarpedCoordinate::WarpedCoordinate(MeshBlock *pmb, ParameterInput *pin,
                                   bool flag)
    : Coordinates(pmb, pin, flag) {
  // Send something to confirm that we are using Warped
  std::cout << "===Note===: Warped coordinates activated" << std::endl;

  // initialize volume-averaged coordinates and spacing
  // x1-direction: x1v = dx/2
  for (int i = il - ng; i <= iu + ng; ++i) {
    // x1v(i) = mean_w2x1(x1f(i+1), x1f(i)) / mean_w2(x1f(i+1), x1f(i));
    x1v(i) = 0.5 * (x1f(i + 1) + x1f(i));
  }
  for (int i = il - ng; i <= iu + ng - 1; ++i) {
      dx1v(i) = x1v(i + 1) - x1v(i);
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

}

// Put in the changes in face2area etc, similar to cylindrical.cpp
// In affine coordinates, the differences lie in face1area and face2area.
// face3area cancels the sqrt(g) term exactly.

//----------------------------------------------------------------------------------------
// FaceXArea functions: compute area of face with normal in X-dir as vector

void WarpedCoordinate::Face1Area(const int k, const int j, const int il,
                                 const int iu, AthenaArray<Real> &area) {
#pragma omp simd
  for (int i = il; i <= iu; ++i) {
    Real &area_i = area(i);
    area_i = w2(x1f(i)) * dx2f(j) * dx3f(k);
  }
  return;
}

void WarpedCoordinate::Face2Area(const int k, const int j, const int il,
                                 const int iu, AthenaArray<Real> &area) {
#pragma omp simd
  for (int i = il; i <= iu; ++i) {
    Real &area_i = area(i);
    area_i = dx1f(i) * dx3f(k);
  }
  return;
}

void WarpedCoordinate::Face3Area(const int k, const int j, const int il,
                                 const int iu, AthenaArray<Real> &area) {
#pragma omp simd
  for (int i = il; i <= iu; ++i) {
    Real &area_i = area(i);
    area_i = mean_w2(x1f(i+1), x1f(i)) * dx1f(i) * dx2f(j);
  }
  return;
}

//----------------------------------------------------------------------------------------
// GetFaceXArea functions: return area of face with normal in X-dir at (i,j,k)

Real WarpedCoordinate::GetFace1Area(const int k, const int j, const int i) {
  return dx2f(j) * dx3f(k) * w2(x1f(i));
}

Real WarpedCoordinate::GetFace2Area(const int k, const int j, const int i) {
  return dx1f(i) * dx3f(k);
}

Real WarpedCoordinate::GetFace3Area(const int k, const int j, const int i) {
  return dx1f(i) * dx2f(j) * mean_w2(x1f(i+1), x1f(i));
}

// Cell Volume function: compute volume of cell as vector

void WarpedCoordinate::CellVolume(const int k, const int j, const int il,
                                  const int iu, AthenaArray<Real> &vol) {
#pragma omp simd
  for (int i = il; i <= iu; ++i) {
    vol(i) = dx1f(i) * dx2f(j) * dx3f(k) * mean_w2(x1f(i+1), x1f(i));
  }
  return;
}

//----------------------------------------------------------------------------------------
// GetCellVolume: returns cell volume at (i,j,k)

Real WarpedCoordinate::GetCellVolume(const int k, const int j, const int i) {
  return dx1f(i) * dx2f(j) * dx3f(k) * mean_w2(x1f(i+1), x1f(i));
}


//----------------------------------------------------------------------------------------
//! Coordinate (Geometric) source term function

void WarpedCoordinate::AddCoordTermsDivergence(
    const Real dt, const AthenaArray<Real> *flux,
    const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &u) {

  for (int k=pmy_block->ks; k<=pmy_block->ke; ++k) {
    for (int j=pmy_block->js; j<=pmy_block->je; ++j) {
#pragma omp simd
      for (int i=pmy_block->is; i<=pmy_block->ie; ++i) {
        Real m_pp = prim(IDN,k,j,i)*prim(IM2,k,j,i)*prim(IM2,k,j,i);
        m_pp += prim(IEN,k,j,i);

        u(IM1,k,j,i) += dt * (dw2(x1f(i + 1), x1f(i))
          / mean_w2(x1f(i + 1), x1f(i))
        ) * m_pp;

      }
    }
  }
  return;
}
