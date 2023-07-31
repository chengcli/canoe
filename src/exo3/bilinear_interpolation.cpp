// C/C++
#include <cmath>

double bilinear_interpolation(double x, double y, double x1, double y1,
                              double x2, double y2, double f_x1_y1,
                              double f_x2_y1, double f_x1_y2, double f_x2_y2) {
  // Calculate the intermediate values
  double f_x_y1 =
      ((x2 - x) / (x2 - x1)) * f_x1_y1 + ((x - x1) / (x2 - x1)) * f_x2_y1;
  double f_x_y2 =
      ((x2 - x) / (x2 - x1)) * f_x1_y2 + ((x - x1) / (x2 - x1)) * f_x2_y2;

  // Calculate the interpolated value at (x, y)
  double interpolatedValue =
      ((y2 - y) / (y2 - y1)) * f_x_y1 + ((y - y1) / (y2 - y1)) * f_x_y2;

  return interpolatedValue;
}
