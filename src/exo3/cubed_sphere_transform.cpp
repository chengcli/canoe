// C/C++
#include <algorithm>

// athena
#include <athena/athena.hpp>

namespace CubedSphereUtility {

#define DBL_EPSILON 1.0e-10

void VecTransABPFromRLL(Real X, Real Y, int blockID, Real U, Real V, Real *V2,
                        Real *V3) {
  Real C = sqrt(1 + X * X);
  Real D = sqrt(1 + Y * Y);
  Real E = sqrt(X * X + Y * Y);
  Real delta = sqrt(1 + X * X + Y * Y);
  switch (blockID) {
    case 2:
    case 3:
    case 4:
    case 6:
      *V2 = X * Y / delta * U - V;
      *V3 = C * D / delta * U;
      break;
    case 1:
      *V2 = -D * Y / (delta * E) * U - D * X / E * V;
      *V3 = C * X / (delta * E) * U - C * Y / E * V;
      break;
    case 5:
      *V2 = D * Y / (delta * E) * U + D * X / E * V;
      *V3 = -C * X / (delta * E) * U + C * Y / E * V;
      break;
  }
};

void VecTransRLLFromABP(Real X, Real Y, int blockID, Real V2, Real V3, Real *U,
                        Real *V) {
  // Calculate the needed parameters
  Real C = sqrt(1 + X * X);
  Real D = sqrt(1 + Y * Y);
  Real E = sqrt(X * X + Y * Y);
  Real delta = sqrt(1 + X * X + Y * Y);

  switch (blockID) {
    case 2:
    case 3:
    case 4:
    case 6:
      *U = delta / (C * D) * V3;
      *V = X * Y / (C * D) * V3 - V2;
      break;
    case 1:
      *U = -delta * Y / (D * E) * V2 + delta * X / (C * E) * V3;
      *V = -X / (D * E) * V2 - Y / (C * E) * V3;
      break;
    case 5:
      *U = delta * Y / (D * E) * V2 - delta * X / (C * E) * V3;
      *V = X / (D * E) * V2 + Y / (C * E) * V3;
      break;
  }
};

void RLLFromXYP(Real dX, Real dY, int nP, Real &lon, Real &lat) {
  switch (nP) {
    // Equatorial panel 2
    case 1:
      lon = atan(dX);
      lat = atan(dY / sqrt(1.0 + dX * dX));
      break;

    // Equatorial panel 4
    case 3:
      lon = atan(dX) + 0.5 * M_PI;
      lat = atan(dY / sqrt(1.0 + dX * dX));
      break;

    // Equatorial panel 6
    case 5:
      lon = atan(dX) + M_PI;
      lat = atan(dY / sqrt(1.0 + dX * dX));
      break;

    // Equatorial panel 3
    case 2:
      lon = atan(dX) + 1.5 * M_PI;
      lat = atan(dY / sqrt(1.0 + dX * dX));
      break;

    // North polar panel
    case 0:
      if (fabs(dX) > DBL_EPSILON) {
        lon = atan2(dX, -dY);
      } else if (dY <= 0.0) {
        lon = 0.0;
      } else {
        lon = M_PI;
      }
      lat = 0.5 * M_PI - atan(sqrt(dX * dX + dY * dY));
      break;

    // South polar panel
    case 4:
      if (fabs(dX) > DBL_EPSILON) {
        lon = atan2(dX, dY);
      } else if (dY > 0.0) {
        lon = 0.0;
      } else {
        lon = M_PI;
      }
      lat = -0.5 * M_PI + atan(sqrt(dX * dX + dY * dY));
      break;
  }

  // Map to the interval [0, 2 pi]
  if (lon < 0.0) {
    lon += 2.0 * M_PI;
  }
}

////////////////////////////////////////////////////////////////////////////////
// May be used later... Note that mo modification of this has been done yet.
void XYPFromRLL(Real lon, Real lat, Real &dX, Real &dY, int &nP) {
  // Default panel to unattainable value
  nP = 6;

  // Translate from RLL coordinates to XYZ space
  Real xx, yy, zz, pm;

  xx = cos(lon) * cos(lat);
  yy = sin(lon) * cos(lat);
  zz = sin(lat);

  pm = std::max(fabs(xx), std::max(fabs(yy), fabs(zz)));

  // Check maxmality of the x coordinate
  if (pm == fabs(xx)) {
    if (xx > 0) {
      nP = 0;
    } else {
      nP = 2;
    }
  }

  // Check maximality of the y coordinate
  if (pm == fabs(yy)) {
    if (yy > 0) {
      nP = 1;
    } else {
      nP = 3;
    }
  }

  // Check maximality of the z coordinate
  if (pm == fabs(zz)) {
    if (zz > 0) {
      nP = 4;
    } else {
      nP = 5;
    }
  }

  // Panel assignments
  Real sx, sy, sz;
  if (nP == 0) {
    sx = yy;
    sy = zz;
    sz = xx;

  } else if (nP == 1) {
    sx = -xx;
    sy = zz;
    sz = yy;

  } else if (nP == 2) {
    sx = -yy;
    sy = zz;
    sz = -xx;

  } else if (nP == 3) {
    sx = xx;
    sy = zz;
    sz = -yy;

  } else if (nP == 4) {
    sx = yy;
    sy = -xx;
    sz = zz;

  } else if (nP == 5) {
    sx = yy;
    sy = xx;
    sz = -zz;
  }

  // Convert to gnomonic coordinates
  dX = sx / sz;
  dY = sy / sz;
}

}  // namespace CubedSphereUtility
