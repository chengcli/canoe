#include <float.h>
#include <math.h>

/*
 * Derived from fcmp(), version 1.2.2, 
 * Copyright (c) 1998-2000 Theodore C. Belding
 * University of Michigan Center for the Study of Complex Systems
 * <mailto:Ted.Belding@umich.edu>
 * <http://fcmp.sourceforge.net>
 *
 * The major modification we have made is to remove the "epsilon" argument
 * and set epsilon inside the fcmp() function.
 *
 * Description:
 *   It is generally not wise to compare two floating-point values for
 *   exact equality, for example using the C == operator.  The function
 *   fcmp() implements Knuth's suggestions for safer floating-point
 *   comparison operators, from:
 *   Knuth, D. E. (1998). The Art of Computer Programming.
 *   Volume 2: Seminumerical Algorithms. 3rd ed. Addison-Wesley.
 *   Section 4.2.2, p. 233. ISBN 0-201-89684-2.
 *
 * Input parameters:
 *   x1, x2: numbers to be compared
 *
 * This routine may be used for both single and double precision.
 *
 * Returns:
 *   -1 if x1 < x2
 *    0 if x1 == x2
 *    1 if x1 > x2		
 */

int fcmp(double x1, double x2) {
  int 
    exponent;
  double
    delta,
    difference;
  const double
    epsilon = DBL_EPSILON;
  
  /* 
   * Get exponent(max(fabs(x1),fabs(x2))) and store it in exponent. 
   *
   * If neither x1 nor x2 is 0,
   * this is equivalent to max(exponent(x1),exponent(x2)).
   *
   * If either x1 or x2 is 0, its exponent returned by frexp would be 0,
   * which is much larger than the exponents of numbers close to 0 in
   * magnitude. But the exponent of 0 should be less than any number
   * whose magnitude is greater than 0.
   *
   * So we only want to set exponent to 0 if both x1 and x2 are 0. 
   * Hence, the following works for all x1 and x2. 
   */
  frexp(fabs(x1) > fabs(x2) ? x1 : x2,&exponent);

  /* 
   * Do the comparison.
   *
   * delta = epsilon*pow(2,exponent)
   *
   * Form a neighborhood around x2 of size delta in either direction.
   * If x1 is within this delta neighborhood of x2, x1 == x2.
   * Otherwise x1 > x2 or x1 < x2, depending on which side of
   * the neighborhood x1 is on.
   */
  delta      = ldexp(epsilon,exponent); 
  difference = x1-x2;

  if (difference > delta) {
    /* x1 > x2 */
    return 1;
  }
  else if (difference < -delta) {
    /* x1 < x2 */
    return -1;
  }
  else  {
    /* -delta <= difference <= delta */
    return 0;  /* x1 == x2 */
  }
}
