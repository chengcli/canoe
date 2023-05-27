#include <stdlib.h>

#include "interpolation.h"

/*! Multidimensional linear interpolation
 * val[0..nval-1]   : output values
 * coor[0..ndim-1]  : coordinate of the interpolation point
 * data[...]        : points to the start position of a multidimensional data
 * table. len[0..ndim-1]   : length of each dimension axis[...]        :
 * coordinates of each dimesnion is placed sequentially in axis
 */
void interpn(double *val, double const *coor, double const *data,
             double const *axis, size_t const *len, int ndim, int nval) {
  int i1, i2;
  i1 = locate(axis, *coor, *len);

  // if the interpolation value is out of bound
  // use the closest value
  if (i1 == -1) {
    i1 = 0;
    i2 = 0;
  } else if (i1 == *len - 1) {
    i1 = *len - 1;
    i2 = *len - 1;
  } else
    i2 = i1 + 1;

  double *v1 = (double *)malloc(nval * sizeof(double));
  double *v2 = (double *)malloc(nval * sizeof(double));

  double x1 = axis[i1];
  double x2 = axis[i2];

  if (ndim == 1) {
    for (int j = 0; j < nval; ++j) {
      v1[j] = data[i1 * nval + j];
      v2[j] = data[i2 * nval + j];
    }
  } else {
    int s = nval;
    for (int j = 1; j < ndim; ++j) s *= len[j];
    interpn(v1, coor + 1, data + i1 * s, axis + *len, len + 1, ndim - 1, nval);
    interpn(v2, coor + 1, data + i2 * s, axis + *len, len + 1, ndim - 1, nval);
  }

  if (x2 != x1)
    for (int j = 0; j < nval; ++j)
      val[j] = ((*coor - x1) * v2[j] + (x2 - *coor) * v1[j]) / (x2 - x1);
  else
    for (int j = 0; j < nval; ++j) val[j] = (v1[j] + v2[j]) / 2.;

  free(v1);
  free(v2);
}

/*! A handy function for one dimensional interpolation
 * x              : interpolation point
 * data[0..len-1] : data array
 * axis[0..len-1] : coordinates
 */
double interp1(double x, double const *data, double const *axis, size_t len) {
  double value;
  interpn(&value, &x, data, axis, &len, 1, 1);
  return value;
}
