
void Thermodynamics::enrollSystemJupiterJuno() {
  if (iv == AMMONIA_VAPOR_ID)
    s = sat_vapor_p_NH3_BriggsS(q[IDN]);
  else if (iv == WATER_VAPOR_ID)
    s = sat_vapor_p_H2O_BriggsS(q[IDN]);
  else
    s = SatVaporPresIdeal(t, p3, beta, delta);
  s /= q[IPR];
}
