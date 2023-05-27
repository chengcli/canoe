#include <stdio.h>
#include <stdlib.h>

#include "interpolation.h"

/*! Multidimensional linear interpolation
 * val[0..nval-1]   : output values
 * coor[0..ndim-1]  : coordinate of the interpolation point
 * data[...]        : points to the start position of a multidimensional data
 * table. len[0..ndim-1]   : length of each dimension axis[...]        :
 * coordinates of each dimesnion is placed sequentially in axis
 */
void interpnf(double *val, double const *coor, double const *data,
              double const *axis, size_t const *len, int ndim) {
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

  double x1 = axis[i1];
  double x2 = axis[i2];
  double v1, v2;

  if (ndim == 1) {
    v1 = data[i1];
    v2 = data[i2];
  } else {
    int s = 1;
    for (int j = 1; j < ndim; ++j) s *= len[j];
    interpnf(&v1, coor + 1, data + i1 * s, axis + *len, len + 1, ndim - 1);
    interpnf(&v2, coor + 1, data + i2 * s, axis + *len, len + 1, ndim - 1);
  }

  if (x2 != x1)
    *val = ((*coor - x1) * v2 + (x2 - *coor) * v1) / (x2 - x1);
  else
    *val = (v1 + v2) / 2.;
}
