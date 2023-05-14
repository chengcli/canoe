#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

//Returns the incomplete gamma function P (a, x) evaluated by its series representation as gamser .
//Also returns ln Γ(a) as gln .
void gser(double *gamser, double a, double x, double *gln)
{
  int n;
  double sum,del,ap;
  *gln=lgamma(a);
  if (x <= 0.0) {
    if (x < 0.0) fprintf(stderr, "x less than 0 in routine gser");
    *gamser=0.0;
    return;
  } else {
    ap=a;
    del=sum=1.0/a;
    for (n=1;n<=ITMAX;n++) {
      ++ap;
      del *= x/ap;
      sum += del;
      if (fabs(del) < fabs(sum)*EPS) {
        *gamser=sum*exp(-x+a*log(x)-(*gln));
      return;
      }
    }
    fprintf(stderr, "a too large, ITMAX too small in routine gser");
    return;
  }
}

// Returns the incomplete gamma function Q(a, x) evaluated by its continued fraction represen-
// tation as gammcf . Also returns ln Γ(a) as gln .
void gcf(double *gammcf, double a, double x, double *gln) {
  int i;
  double an,b,c,d,del,h;
  *gln=lgamma(a);
  b=x+1.0-a;
  //Set up for evaluating continued fraction
  //by modified Lentz’s method (§5.2)
  c=1.0/FPMIN;
  //with b 0 = 0.
  d=1.0/b;
  h=d;
  for (i=1;i<=ITMAX;i++) {
  //Iterate to convergence.
    an = -i*(i-a);
    b += 2.0;
    d=an*d+b;
    if (fabs(d) < FPMIN) d=FPMIN;
    c=b+an/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (fabs(del-1.0) < EPS) break;
  }
  if (i > ITMAX) fprintf(stderr, "a too large, ITMAX too small in gcf");
  *gammcf=exp(-x+a*log(x)-(*gln))*h;
}

//Returns the incomplete gamma function P (a, x).
double gammp(double a, double x) {
  double gamser,gammcf,gln;
  if (x < 0.0 || a <= 0.0) fprintf(stderr, "Invalid arguments in routine gammp");
  //Use the series representation.
  if (x < (a+1.0)) {
    gser(&gamser,a,x,&gln);
    return gamser;
  } else { // Use the continued fraction representation
    gcf(&gammcf,a,x,&gln);
  return 1.0-gammcf; //and take its complement.
  }
}

// Returns the incomplete gamma function Q(a, x) ≡ 1 − P (a, x).
double gammq(double a, double x) {
  double gamser,gammcf,gln;
  if (x < 0.0 || a <= 0.0) fprintf(stderr, "Invalid arguments in routine gammq");
  // Use the series representation
  if (x < (a+1.0)) {
    gser(&gamser,a,x,&gln);
    return 1.0-gamser;  //and take its complement.
  } else { //Use the continued fraction representation.
    gcf(&gammcf,a,x,&gln);
    return gammcf;
  }
}

#undef ITMAX
#undef EPS
#undef FPMIN
