
// Calculates saturation adjustment of
// NH3 + H2S -> NH4SH
//
void Thermodynamics::condensation_NH3_H2S_NH4SH(Variable *qfrac) const {
  int j1 = m_id[0], j2 = m_id[1], j3 = m_id[2];

  double a = 1, b = 1;
  double a = m_stoi[0] / m_stoi[2], b = m_stoi[1] / m_stoi[2];

  double ptol = air.get_pa(), temp = air.get_temp(),
         gmol = get_gas_mols(air) - air.mgas(j1) - air.mgas(j2);

  double x0 = a * air.mgas(j3) + air.mgas(j1),
         y0 = b * air.mgas(j3) + air.mgas(j2), z0 = 0.,
         pp1 = x0 / (x0 + y0 + gmol) * ptol, pp2 = y0 / (x0 + y0 + gmol) * ptol;

  double mol0 = air.mgas(j1), mol1 = air.mgas(j2);

  if (pp1 * pp2 > sat_vapor_p(temp)) {
    double zeta = sat_vapor_p(temp) / (ptol * ptol), rh;

    int error;

    if (a * y0 > b * x0) {
      error = _root(0., 1., 1.E-8, &rh, [a, b, x0, y0, gmol, zeta](double r) {
        double x = r * x0, y = y0 - b / a * (x0 - x);
        return (x * y) / _sqr(x + y + gmol) - zeta;
      });
      y0 -= b / a * (1. - rh) * x0;
      z0 = (1. - rh) * x0 / a;
      x0 *= rh;
    } else {
      error = _root(0., 1., 1.E-8, &rh, [a, b, x0, y0, gmol, zeta](double r) {
        double y = r * y0, x = x0 - a / b * (y0 - y);
        return (x * y) / _sqr(x + y + gmol) - zeta;
      });
      x0 -= a / b * (1. - rh) * y0;
      z0 = (1. - rh) * y0 / b;
      y0 *= rh;
    }

    /* correction for y0 << x0 and x0 >> y0 */
    if (y0 < 1.E-8 * x0) {
      y0 = (x0 + gmol) * (x0 + gmol) / x0 * zeta;
      z0 = (b * air.mgas(j3) + air.mgas(j2) - y0) / b;
      x0 = a * air.mgas(j3) + air.mgas(j1) - a * z0;
    } else if (x0 < 1.E-8 * y0) {
      x0 = (y0 + gmol) * (y0 + gmol) / y0 * zeta;
      z0 = (a * air.mgas(j3) + air.mgas(j1) - x0) / a;
      y0 = b * air.mgas(j3) + air.mgas(j2) - b * z0;
    }

    if (error) {
      std::cerr << "** AerosolCondensation failed. **" << std::endl;
      std::cerr << "** Before equilibration: **" << std::endl;
      std::cerr << air;
      std::cerr << "** Reaction: **" << std::endl;
      std::cerr << *this;
      assert(0);
    }
  }

  air.mgas(j1) = x0;
  air.mgas(j2) = y0;
  air.mgas(j3) = z0;

  return (fabs(x0 - mol0) < PRECISION) && (fabs(y0 - mol1) < PRECISION);
}
