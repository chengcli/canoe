/**************************************************************************

                    =====  JUNO codes  =====

   Copyright (C) 2019 California Institute of Technology (Caltech)
   All rights reserved.

   Written by Fabiano Oyafuso (fabiano@jpl.nasa.gov)
   Adapted by Cheng Li (cli@gps.caltech.edu) to SNAP structure

**************************************************************************/

// C/C++ headers
#include <cmath>

// h2s data headers
#include "mwr_absorber_h2s.hpp"

// cli, 190801
inline double SQUARE(double x) { return x * x; }
static double const pi = 3.14159265358979;
// end

inline double Lineshape_BR(double obsfreq, double linefreq, double gamma,
                           double zeta, double delta) {
  // double const pi = 2.0*acos(0.0);
  double const c = 2.99792458e10;

  double obsfreqghz = obsfreq / 1000.0;
  double linefreqghz = linefreq / 1000.0;

  double u1 = 2.0 / pi * SQUARE(obsfreqghz / linefreqghz);
  double u2 = SQUARE(linefreqghz + delta) + SQUARE(gamma) - SQUARE(zeta);

  u1 *= ((gamma - zeta) * SQUARE(obsfreqghz) + (gamma + zeta) * u2);
  u2 = SQUARE(obsfreqghz) - SQUARE(linefreqghz + delta) - SQUARE(gamma) +
       SQUARE(zeta);
  u2 = SQUARE(u2) + 4.0 * SQUARE(obsfreqghz) * SQUARE(gamma);

  // The factor c/1.0e9 converts 1/GHz to 1/(cm**-1):
  return (u1 / u2) * c / 1.0e9;
}

double attenuation_H2S_Hofstadter(double freq, double P, double T, double XH2,
                                  double XHe, double XH2S) {
  /*
    code adpated from Mark Hofstadter.
    Here are some of his comments:

  * This subroutine calculates H2S opacity, as described in DeBoer
  * and Steffes 1994 (Icarus 109, 352-366).  Inputs are local atmospheric
  * conditions and observing frequency, output is the opacity in cm-1.
  *
  *
  * COMMENTS:
  * H2S line information is passed in via the parameter list, presumably
  * the calling program has also called subroutine ReadLineInfo.
  *
  * Line widths have this functional form (Eq. 20 of DeBoer and Steffes
  * 1994):
  *     delnu = SUM_i [ delnu_o_i * P_i * (T_o/T)**Xi_i ]
  * where the sum is over all gasses present, the subscript "o" refers
  * to a reference state, the subscript "i" to the broadening gas in
  * question, and T is temperature.  DeBoer and Steffes, however, only
  * specify parameters for H2, He, and H2S, so I am forced to neglect
  * broadening by CH4, and they suggest using Xi_i = 0.7 for all species.
  *
  * The absorption at the line center, A, is calculated with Equation 19
  * from DeBoer and Steffes, having units of cm**-1.  Note that in
  * their equation, P (which is the H2S partial pressure, not the total
  * atmosphere) is in bar (the leading factor of 1E6 will convert
  * bar to dyne/cm**2) and delnu must be converted from GHz
  * to cm**-1.  The line intensity in their equation is in units
  * of cm**-1 / (molecules/cm**2), requiring a conversion from the
  * Poynter and Picket units of nm**2 MHz.  The conversion factor is
  * to multiply the PandP value by 1E-8 (converting nm**2 MHz to cm**2 Hz)
  * and divide by the speed of light (converting Hz to cm**-1).  The
  * "per molecule" is just assumed in some way.  This intensity unit
  * conversion is set to the value PPunitsToDS.
  *
  * This software was tested against Karuthika Devaraj's code during a
  * visit to the VLA in June 2008 (while she was a summer intern there).
  * We were satisfied that our results matched to the precision of our
  * machines.  Our results did differ from those published by DeBoer and
  * Steffes, but we believe this is due to that paper using
  * and older version of the JPL ("Poynter and Pickett" for us
  * old-timers) line catalog.  In my notes, I state that the differences
  * between Karuthika and I were ~1% near line centers, 10x smaller
  * elsewhere (though it's possible that I wrote that before we
  * discovered that we were comparing calculations at slightly difference
  * frequencies, which made us think differences near line centers were
  * greater than they really were).
  *
  *
  * KNOWN BUGS AND LIMITATIONS:
  * 1) Does not account for line broadening due to collisions with
  *    CH4.
  *
  * 2) I can leave F and delnu in GHz in the final expression rather
  *    than converting both to cm**-1 cuz the units will cancel.
  *
  *
  * SUBROUTINES CALLED:
  *   None.
  *
  * FUNCTIONS CALLED:
  *   Lineshape_BR (real):  Returns the Ben Reuven line-shape function
  *                         as defined in DeBoer and Steffes Eq. 23, in
  *                         units of 1/(cm**-1) (note the double
  *                         inverse!).  See the Function for required
  *                         inputs.
  *
  *
  * INTERNALLY SET PARAMETERS
  *   MAXNUMLINE (integer): Max number of spectral lines in input file.
  *
  *   NUMPARTFUNC (integer): The number of temperatures for which
  *                            information on the partition function
  *                            is provided (to dimension qlog and tqlog).
  *
  *   NUMQUANTNUM (integer):  The number of quantum numbers used to
  *                           specify molecular states in the input file.
  *
  *
  * INPUTS:
  *    EXTMAXLINE (integer): Externally used value to dimension arrays.
  *                          Is tested to verify it matches MAXNUMLINE.
  *
  *    EXTPARTFUNC (integer): Externally used value to dimension arrays.
  *                          Is tested to verify it matches NUMPARTFUNC.
  *
  *    EXTQUANTNUM (integer): Externally used value to dimension arrays.
  *                          Is tested to verify it matches NUMQUANTNUM.
  *
  *    nline (integer):  Number of H2S lines to model.
  *
  *    qlog (real):  NUMPARTFUNC-element array of the base 10 log of
  *                  the partition function for the temperatures
  *                  specifed in tqlog.
  *
  *    tqlog (real): NUMPARTFUNC-element array of the
  *                  temperatures (K) at which the partition functions
  *                  in qlog are calculated.
  *
  *    freq (real):  Array of line frequencies, in MHz.
  *
  *    freqerr (real):  Array of line frequency uncertainties, in MHz.
  *
  *    lgint (real):  Array of base 10 log of integrated intensity in
  *                   units of nm^2 MHz at 300 K.
  *
  *    dr (integer):  Array of degrees of freedom in the rotational
  *                   partition function.
  *
  *    elo (real):  Array of lower energy state (cm-1) relative to
  *                 ground state.
  *
  *    gup (integer): Array of upper state degeneracy.
  *
  *    qnfmt (integer): Array specifying the format for QNp and QNpp.
  *
  *    qnp (integer):  Array (NUMQUANTNUMxMAXNUMLINE) of upper-state
  *                    quantum numbers.
  *
  *    qnpp (integer):  Array (NUMQUANTNUMxMAXNUMLINE) of lower-state
  *                     quantum numbers.
  *
  *    ObsFreq (real):  The observing frequency, in MHz.
  *
  *    T (real):  Atmospheric temperature, in Kelvin.
  *
  *    P (real):  Total atmospheric pressure, in mbar.
  *
  *    Ph2/he/ch4/h2o/nh3/h2s (real):  Partial pressures (mbar) of
  *                                    H2, He, CH4, H2O, NH3, and H2S.
  *
  *
  * OUTPUTS:
  *
  *  kh2s (real):  The total H2S opacity, in 1/cm.
  *
  *
  * LOCAL VARIABLES:
  *  pi (real):  3.1415926.....
  *
  *  kb (real):  Boltzmann's constant in cgs units (erg/K).  Taken from
  *              Cohen and Taylor 1995.
  *
  *  c (real):  Speed of light, in cm/s.  From Cohen and Taylor 1995.
  *
  *  h (real):  Planck constant, in erg sec.  From Cohen and Taylor 1995.
  *
  *  I (integer):  Loop counter.
  *
  *  Xi (real):  Line broadening temperature dependence.  As per DeBoer
  *              and Steffes, it is assumed to be 0.7 for all species.
  *              (I see, however, that the Joiner et al. 1992 reference
  *              DS quote really used 0.67.  Using it makes a
  *              negligable difference.)
  *
  *  To_linewidth (real):  The reference temperature for calculating
  *                        delnu (K).  This is never explicitly given
  *                        in the DS paper but the Joiner et al. 1992
  *                        reference (IEEE Trans. Micro. Theory and
  *                        Tech.) indicates they used 296 K.  I use 296,
  *                        but I tried both it and 300 K and it made
  *                        only a trivial difference.
  *
  *  delnu (real):  Total line broadening parameter (GHz) of H2S
  *                  due to H2, He, and H2S.
  *
  *  delnuh2_o (real):  Line broadening parameter of H2S by H2, from
  *                      DeBoer and Steffes (GHz/bar).
  *
  *  delnuhe_o (real):  Line broadening parameter of H2S by He, from
  *                      DeBoer and Steffes (GHz/bar).
  *
  *  delnuh2s_o (real):  Line broadening parameter of H2S by H2S, from
  *                       DeBoer and Steffes (GHz/bar) EXCEPT for 4
  *                       millimeter lines, where other values are used.
  *
  *  delnuh2s_168ghz (real):  Line broadening parameter of H2S by H2S,
  *                       (GHz/bar) from DeBoer and Steffes for the
  *                       168 GHz line.  (They reference Helminger and
  *                       DeLucia 1977, JQSRT 17, 751-754.)
  *
  *  delnuh2s_216ghz (real):  Line broadening parameter of H2S by H2S,
  *                       (GHz/bar) from DeBoer and Steffes for the
  *                       216 GHz line.  (They reference Helminger and
  *                       DeLucia 1977, JQSRT 17, 751-754.)
  *
  *  delnuh2s_300ghz (real):  Line broadening parameter of H2S by H2S,
  *                       (GHz/bar) from DeBoer and Steffes for the
  *                       300 GHz line.  (They reference Helminger and
  *                       DeLucia 1977, JQSRT 17, 751-754.)
  *
  *  delnuh2s_393ghz (real):  Line broadening parameter of H2S by H2S,
  *                       (GHz/bar) from DeBoer and Steffes for the
  *                       393 GHz line.  (They reference Helminger and
  *                       DeLucia 1977, JQSRT 17, 751-754.)
  *
  *  delnuh2s (real):  Line broadening parameter of H2S by H2S to use
  *                    for a particular line, set to one of the above
  *                    H2S values.
  *
  *  eta (real):  Constant reflecting the temperature dependence of
  *               the partition function.  It is 1 for linear molecules,
  *               and ~1.5 for non-linear and symmetric top molecules.
  *               Here it is set to 1.5.  I believe H2S is an
  *               asymmetric top molecule (LeSueur et al., 1992,
  *               Molec. Phys. 76, 1147-1156) so I was not sure whether
  *               1 or 1.5 is a better value to pick.  Since the PandP
  *               catalog indicates there are 3 degrees of freedom
  *               for H2S, it makes sense to use 1.5.
  *
  *  To_intensity (real):  Reference temperature for line intensity
  *                        value in the line catalog (300 K).
  *
  *  A (real):  Line intensity, in cm**-1.
  *
  *  PPunitsToDSunits (real):  Conversion factor for taking the line
  *                            intensity units of the PandP catlog:
  *                                 nm**2 MHz
  *                            to what DS want:
  *                                 cm**-1/(molec/cm**2).
  *                            It is 1 / (c*1.0E8), where the 1E-8
  *                            converts PandP units to cgs, and the
  *                            factor of c converts Hz to cm**-1.
  *
  *  zeta (real):  The coupling parameter in the Ben-Reuven lineshape.
  *                As per DandS, it is set equal to delnu (GHz).
  *
  *  delta (real):  The pressure-shift term in the Ben-Reuven lineshape.
  *                 As per DandS, it is set equal to
  *                    1.28 * Ph2s/1000
  *                 which I think is GHz for Ph2s in bar.
  *
  *  F (real):  The line shape function in cm (think of it as 1/cm**-1).

   */

  double T_in = T;
  // pressures in mbar
  double Ph2s = 1e3 * P * XH2S;
  double Phe = 1e3 * P * XHe;
  double Ph2 = 1e3 * P * XH2;

  int const Nlines(sizeof(h2s_lines) / sizeof(h2s_lines[0]));

  // double const pi = 2.0*acos(0.0);

  // Physical constants in cgs from Cohen and Taylor 1986/1987.
  double const kb = 1.38066e-16;
  double const c = 2.99792458e10;
  double const h = 6.626075e-27;

  // Other constants:
  double const PPunitsToDSunits = 1.0 / (c * 1.0e8);

  // Parameters from DeBoer and Steffes or references therein.
  // (The calculation of delta has a divisor of 1000 to convert
  //  pressure from mbar to bar.)
  double const Xi = 0.7;
  double const To_linewidth = 296.0;
  double const delnuh2_o = 1.96;
  double const delnuhe_o = 1.20;
  double const delnuh2s_o = 5.78;
  double const delnuh2s_168ghz = 5.38;
  double const delnuh2s_216ghz = 6.82;
  double const delnuh2s_300ghz = 5.82;
  double const delnuh2s_393ghz = 5.08;
  double const eta = 1.5;
  double const To_intensity = 300.0;
  double const delta = 1.28 * Ph2s / 1000.0;

  double kh2s = 0.0;

  //  Calculate opacity from each H2S line.
  for (auto I = 0; I < Nlines; I++) {
    //  Calculate linewidth (delnu) as per DeBoer and Steffes Eq. 20.
    //  Factors of 1000 convert mbar to bar.
    double delnu = delnuh2_o * (Ph2 / 1000.0) * pow(To_linewidth / T_in, Xi) +
                   delnuhe_o * (Phe / 1000.0) * pow(To_linewidth / T_in, Xi);

    double delnuh2s;
    if (int(h2s_lines[I].FREQ) == 168762)
      delnuh2s = delnuh2s_168ghz;
    else if (int(h2s_lines[I].FREQ) == 216710)
      delnuh2s = delnuh2s_216ghz;
    else if (int(h2s_lines[I].FREQ) == 300505)
      delnuh2s = delnuh2s_300ghz;
    else if (int(h2s_lines[I].FREQ) == 393450)
      delnuh2s = delnuh2s_393ghz;
    else
      delnuh2s = delnuh2s_o;

    delnu += delnuh2s * (Ph2s / 1000.0) * pow(To_linewidth / T_in, Xi);

    //    Calculate absorption at line center (A), as per DeBoer and
    //    Steffes Eq. 19.  The leading factor of 1E3 converts P in mbar
    //    to dyne/cm**2 (DS had a factor of 1E6 because they give P
    //    in bar).  The factor of 1E9/c converts
    //    delnu from GHz to cm**-1, as specified by DS.
    double A = 1.0e3 * Ph2s / (pi * kb * T_in * delnu * 1.0e9 / c);
    A *= pow(To_intensity / T_in, eta + 1.0);
    A *= pow(10.0, h2s_lines[I].LGINT) * PPunitsToDSunits;
    A *= exp(-(h * c * h2s_lines[I].ELO / kb) *
             (1.0 / T_in - 1.0 / To_intensity));

    //    Calculate lineshape (F)

    //    As per DandS, set zeta to delnu
    double zeta = delnu;
    double F = Lineshape_BR(freq * 1e3, h2s_lines[I].FREQ, delnu, zeta, delta);

    //    Calculate this line's contribution to the final opacity
    //    (the factor of 1E9/c converts delnu in GHz to cm**-1):
    kh2s += A * pi * delnu * (1.0e9 / c) * F;
  }

  return kh2s;  // inv(cm)
}
