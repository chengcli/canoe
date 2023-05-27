#ifndef INTERPOLATION_H_
#define INTERPOLATION_H_

#ifdef __cplusplus
extern "C" {
#endif
#include <stddef.h>

struct float_triplet;

int locate(double const *xx, double x, int n);
void interpn(double *val, double const *coor, double const *data,
             double const *axis, size_t const *len, int ndim, int nval);
double interp1(double x, double const *data, double const *axis, size_t len);

void spline(int n, struct float_triplet *table, double y1_bot, double y1_top);
int find_place_in_table(int n, struct float_triplet *table, double x,
                        double *dx, int il);
double splint(double xx, struct float_triplet *table, double dx);

void interpnf(double *val, double const *coor, double const *data,
              double const *axis, size_t const *len, int ndim);

#endif

#ifdef __cplusplus
} /* extern "C" */
#endif
