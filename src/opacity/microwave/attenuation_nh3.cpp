/**************************************************************************

                    =====  JUNO codes  =====

   Copyright (C) 2019 California Institute of Technology (Caltech)
   All rights reserved.

   Written by Fabiano Oyafuso (fabiano@jpl.nasa.gov)
   Adapted by Cheng Li (cli@gps.caltech.edu) to SNAP structure

**************************************************************************/

// C/C++ headers
#include <cmath>

// nh3 data headers
#include "mwr_absorber_nh3.hpp"

// cli, 190801
inline double SQUARE(double x) { return x*x; }
inline double atm_from_bar(double P) { return 0.9869232667160128 * P; }
static double const pi = 3.14159265358979;
// end

double attenuation_NH3_Bellotti_switch(double freq, double P, double P_idl,
       double T, double XH2, double XHe, double XNH3, double XH2O)
//  sum the contributions of all lines, first the Inversion transitions,
//  then the Rotational transitions

{
   #if defined(FAKE_CUTOFF_TO_AVOID_SINGULARITY_IN_BELLOTTI_ABS)
   // Bellotti's parameterization suffers from two problems at high pressure:
   // (1) for P~7000 bar, absorption becomes negative!!!
   // (2) slope of opacity wrt pressure changes sign around P~1000 bar.
   // We patch this by readjusting opacity for T>T_cutoff.  (We choose a temperature
   // cutoff instead of a pressure cutoff, since the opacity is more strongly dependent
   // on temperature.)
   // Update: this anomalous behavior has been corrected by zeroing zhe
   constexpr double T_Bellotti_cutoff = 1250;
   if (T>T_Bellotti_cutoff) {
      T = T_Bellotti_cutoff;
   }
   #endif

   // TEMPORARY KLUDGE to compare to Hanley, fig19!!!
   //T=295;
   //XNH3=0.004366812227074;
   bool use_mm_parameters(freq>= 30 ? true : false);
   // miscellaneous constants:
   double const ppscale = 2.99792458e18; // converts to cm^2*cm-1 from nm^2*MHz
   double const gcon = 25.923;
   double const cmghz = 29.9792458;    // converts to GHz from cm-1
   double const hckt0 = 0.0047959232;  // hc/kT0, where T0=300K
   double const cdens = 7.242971565e21;  // 1/kB * (1e5 Pa/bar) * (1e-6 m^-3/cm^-3)

   // partial pressures in bars
   double const PNH3 = P_idl* XNH3;
   double const PH2  = P_idl* XH2;
   double const PHe  = P_idl* XHe;
   double const PH2O = P_idl* XH2O;

   double const dens = cdens*(P_idl*XNH3)/T;
   //double const dens = cdens*PNH3/T;
   double const beta300 = 300.0/T;
   double const beta295 = 295.0/T;
   double const beta300_25 = pow(beta300,2.5);
   double gh2, zh2, ghe, zhe, gnh3, znh3, gh2o, zh2o, dd, Dinv;

   if(freq> 30){
      gh2  = 1.7465 * PH2  * pow(beta300,0.8202);
      ghe  = 0.9779 * PHe  * pow(beta300,1.0);
      gnh3 = 0.7298 * PNH3 * pow(beta295,1.0);
      gh2o = 4.8901 * PH2O * pow(beta300,0.5);

      zh2  = 1.2163 * PH2  * pow(beta300,0.8873);
      // 161014: Bellotti says zhe should be zero to prevent anomolous behavior at high temp
      //zhe  = 0.0291 * PHe  * pow(beta300,0.8994);
      zhe  = 0.0;
      znh3 = 0.5152 * PNH3 * pow(beta295,0.6667);
      zh2o = 2.731  * PH2O * pow(beta300,1.0);

      dd = -0.0627;      // pressure shift factor from linewidth
      Dinv = 0.9862;      // unitless emperical scale factor
   }else{
      // coefficients from Bellotti (private correspondence 160407)

      gh2  = 1.69369709 * PH2  * pow(beta300,0.80847492);
      ghe  = 0.69968615 * PHe  * pow(beta300,1.);
      gnh3 = 0.75232086 * PNH3 * pow(beta295,1.0);
      //gh2o = 4.8901 * PH2O * pow(beta300,0.5);
      gh2o = 5.2333 * PH2O * pow(beta300,0.6224);

      zh2  = 1.32634109 * PH2  * pow(beta300,0.81994711);
      // 161014: Bellotti says zhe should be zero to prevent anomolous behavior at high temp
      //zhe  = 0.16068662 * PHe  * pow(beta300,-0.72691007);
      zhe  = 0.0;
      znh3 = 0.61622276 * PNH3 * pow(beta295,1.38318269);
      //zh2o = 2.731  * PH2O * pow(beta300,1.0);
      zh2o = 5.2333  * PH2O * pow(beta300,2.1248);

      dd = -0.01386851;      // pressure shift factor from linewidth
      Dinv = 0.96190579;      // unitless emperical scale factor
   }



   //  Inversion:
   double alpha_inv(0.0);
   int const Ninv = sizeof(NH3_inv_JPL5)/sizeof(NH3_inv_JPL5[0]);
   for(int i=0; i<Ninv; i++) {
      LineParameters_NH3_inversion_JPL5 const & x = NH3_inv_JPL5[i];
      double f = x.f0;  // convert to GHz
      double fx = freq/f;
      double sjkt = x.I0 * beta300_25 * exp(hckt0*(x.E0)*(1.0-beta300));

      // pressure-broadening widths:
      double gam = gh2 + ghe + gnh3 * x.gamma_self;
      double zet = zh2 + zhe + znh3 * x.gamma_self;

      if (PH2O != 0.0) { // Broadening line width's due to water Vapor:
         gam += gh2o;
         zet += zh2o;
      }

      // Ben-Reuven line profile (eqn 2.6):
      double num = (gam-zet)*SQUARE(freq) + (gam+zet)*( SQUARE(f+dd*gam) + SQUARE(gam) - SQUARE(zet) );
      double den = SQUARE(SQUARE(freq) - SQUARE(f+dd*gam) - SQUARE(gam) + SQUARE(zet)) +
         SQUARE(2.0*freq*gam);
      double fline = cmghz*2.0*fx*num/(pi*den);

      alpha_inv += Dinv*dens*sjkt*fx*fline;
   }

   double alpha_final;

   if(freq> 30){
      //  Rotational:
      double alpha_rot(0.0);
      double const gH2_rot  = 0.2984 * PH2  * pow(beta300,0.8730);
      double const gHe_rot  = 0.75   * PHe  * pow(beta300,0.6667);
      double const gNH3_rot = 3.1789 * PNH3 * pow(beta300,1.0);
      double const Drot = 2.4268;

      int const Nrot = sizeof(NH3_rot_JPL5)/sizeof(NH3_rot_JPL5[0]);
      for(int i=0; i<Nrot; i++) {
         LineParameters_NH3_rotational_JPL5 const & x = NH3_rot_JPL5[i];
         double f = x.f0;  // convert to GHz
         double fx = freq/f;
         // pressure-broadening width:
         double gam = gH2_rot*(x.gamma_H2) + gHe_rot*(x.gamma_He) + gNH3_rot*(x.gamma_self);

         if (PH2O != 0.0) {
            gam += gh2o;
         }
         // Gross line profile (from Joiner at al.):
         double fline = 4.0/pi*cmghz*fx*gam/(SQUARE(f-fx*freq) + SQUARE(2.0*fx*gam));
         // the rest is as for Inversion case:
         alpha_rot += Drot*dens*(x.I0)*fx*beta300_25 * exp(hckt0*(x.E0)*(1.0-beta300))*fline;
      }

      //  Rotovibrational:
      double alpha_rv(0.0);
      double const gH2_nu2  = 1.4  * PH2  * pow(beta300,0.73);
      double const gHe_nu2  = 0.68 * PHe  * pow(beta300,0.5716);
      double const gNH3_nu2 = 9.5  * PNH3 * pow(beta300,1.0);
      double const Drv = 1.1206;

      int const Nrv = sizeof(NH3_rv_JPL5)/sizeof(NH3_rv_JPL5[0]);
      for(int i=0; i<Nrv; i++) {
         LineParameters_NH3_rotovibrational_JPL5 const & x = NH3_rv_JPL5[i];
         double f = x.f0;  // convert to GHz
         double fx = freq/f;
         // pressure-broadening width:
         double gam = gH2_nu2 + gHe_nu2 + gNH3_nu2;

         if (PH2O != 0.0) {
            gam += gh2o;
         }
         // Gross line profile (from Joiner at al.):
         double fline = 4.0/pi*cmghz*fx*gam/(SQUARE(f-fx*freq) + SQUARE(2.0*fx*gam));
         // the rest is as for Inversion case:
         alpha_rv += Drv*dens*(x.I0)*fx*beta300_25 * exp(hckt0*(x.E0)*(1.0-beta300))*fline;
      }
      alpha_final = alpha_inv + alpha_rot + alpha_rv; // in cm^-1

   }else{
      alpha_final = alpha_inv;
   }

   return alpha_final;
}

double attenuation_NH3_Bellotti(double freq, double P, double P_idl, double T,
       double XH2, double XHe, double XNH3, double XH2O)
//  sum the contributions of all lines, first the Inversion transitions,
//  then the Rotational transitions

{
   #if defined(FAKE_CUTOFF_TO_AVOID_SINGULARITY_IN_BELLOTTI_ABS)
   // Bellotti's parameterization suffers from two problems at high pressure:
   // (1) for P~7000 bar, absorption becomes negative!!!
   // (2) slope of opacity wrt pressure changes sign around P~1000 bar.
   // We patch this by readjusting opacity for T>T_cutoff.  (We choose a temperature
   // cutoff instead of a pressure cutoff, since the opacity is more strongly dependent
   // on temperature.)
   // Update: this anomalous behavior has been corrected by zeroing zhe
   constexpr double T_Bellotti_cutoff = 1250;
   if (T>T_Bellotti_cutoff) {
      T = T_Bellotti_cutoff;
   }
   #endif

   // TEMPORARY KLUDGE to compare to Hanley, fig19!!!
   //T=295;
   //XNH3=0.004366812227074;

   // miscellaneous constants:
   double const ppscale = 2.99792458e18; // converts to cm^2*cm-1 from nm^2*MHz
   double const gcon = 25.923;
   double const cmghz = 29.9792458;    // converts to GHz from cm-1
   double const hckt0 = 0.0047959232;  // hc/kT0, where T0=300K
   double const cdens = 7.242971565e21;  // 1/kB * (1e5 Pa/bar) * (1e-6 m^-3/cm^-3)

   // partial pressures in bars
   double const PNH3 = P_idl* XNH3;
   double const PH2  = P_idl* XH2;
   double const PHe  = P_idl* XHe;
   double const PH2O = P_idl* XH2O;

   double const dens = cdens*(P_idl*XNH3)/T;
   //double const dens = cdens*PNH3/T;
   double const beta300 = 300.0/T;
   double const beta295 = 295.0/T;
   double const beta300_25 = pow(beta300,2.5);
   double gh2, zh2, ghe, zhe, gnh3, znh3, gh2o, zh2o, dd, Dinv;

   // coefficients from Bellotti (private correspondence 160407)
   gh2  = 1.6948 * PH2  * pow(beta300,1.0);
   ghe  = 0.6385 * PHe  * pow(beta300,0.7995);
   gnh3 = 0.7 * PNH3 * pow(beta295,1.0);
   //gh2o = 4.8901 * PH2O * pow(beta300,0.5);
   gh2o = 5.2333 * PH2O * pow(beta300,0.6224);

   zh2  = 1.2718 * PH2  * pow(beta300,0.9988);
   zhe  = 0.1316 * PHe  * pow(beta300,-0.9129);
   znh3 = 0.4879 * PNH3 * pow(beta295,0.7885);
   //zh2o = 2.731  * PH2O * pow(beta300,1.0);
   zh2o = 5.2333  * PH2O * pow(beta300,2.1248);

   dd = 0.;      // pressure shift factor from linewidth
   Dinv = 1.0323;      // unitless emperical scale factor

   //  Inversion:
   double alpha_inv(0.0);
   int const Ninv = sizeof(NH3_inv_JPL5)/sizeof(NH3_inv_JPL5[0]);
   for(int i=0; i<Ninv; i++) {
      LineParameters_NH3_inversion_JPL5 const & x = NH3_inv_JPL5[i];
      double f = x.f0;  // convert to GHz
      double fx = freq/f;
      double sjkt = x.I0 * beta300_25 * exp(hckt0*(x.E0)*(1.0-beta300));

      // pressure-broadening widths:
      double gam = gh2 + ghe + gnh3 * x.gamma_self;
      double zet = zh2 + zhe + znh3 * x.gamma_self;

      if (PH2O != 0.0) { // Broadening line width's due to water Vapor:
         gam += gh2o;
         zet += zh2o;
      }

      // Ben-Reuven line profile (eqn 2.6):
      double num = (gam-zet)*SQUARE(freq) + (gam+zet)*( SQUARE(f+dd*gam) + SQUARE(gam) - SQUARE(zet) );
      double den = SQUARE(SQUARE(freq) - SQUARE(f+dd*gam) - SQUARE(gam) + SQUARE(zet)) +
         SQUARE(2.0*freq*gam);
      double fline = cmghz*2.0*fx*num/(pi*den);

      alpha_inv += Dinv*dens*sjkt*fx*fline;
   }
   //  Rotational:
   double alpha_rot(0.0);

   double const gH2_rot  = 1.7761 * PH2  * pow(beta300,0.5);
   double const gHe_rot  = 0.6175   * PHe  * pow(beta300,0.5663);
   double const gNH3_rot = 3.1518 * PNH3 * pow(beta300,1.0);
   double const Drot = 2.7252;

   int const Nrot = sizeof(NH3_rot_JPL5)/sizeof(NH3_rot_JPL5[0]);
   for(int i=0; i<Nrot; i++) {
      LineParameters_NH3_rotational_JPL5 const & x = NH3_rot_JPL5[i];
      double f = x.f0;  // convert to GHz
      double fx = freq/f;
      // pressure-broadening width:
      double gam = gH2_rot*(x.gamma_H2) + gHe_rot*(x.gamma_He) + gNH3_rot*(x.gamma_self);
      if (PH2O != 0.0) {
         gam += gh2o;
      }
      // Gross line profile (from Joiner at al.):
      double fline = 4.0/pi*cmghz*fx*gam/(SQUARE(f-fx*freq) + SQUARE(2.0*fx*gam));
      // the rest is as for Inversion case:
      alpha_rot += Drot*dens*(x.I0)*fx*beta300_25 * exp(hckt0*(x.E0)*(1.0-beta300))*fline;
   }

   //  Rotovibrational:
   double alpha_rv(0.0);

   double const gH2_nu2  = 0.5982 * PH2  * pow(beta300,0.5);
   double const gHe_nu2  = 0.6175 * PHe  * pow(beta300,0.5505);
   double const gNH3_nu2 = 5.0894  * PNH3 * pow(beta300,0.9996);
   double const Drv = 0.7286;

   int const Nrv = sizeof(NH3_rv_JPL5)/sizeof(NH3_rv_JPL5[0]);
   for(int i=0; i<Nrv; i++) {
      LineParameters_NH3_rotovibrational_JPL5 const & x = NH3_rv_JPL5[i];
      double f = x.f0;  // convert to GHz
      double fx = freq/f;
      // pressure-broadening width:
      double gam = gH2_nu2 + gHe_nu2 + gNH3_nu2;
      if (PH2O != 0.0) {
         gam += gh2o;
      }
      // Gross line profile (from Joiner at al.):
      double fline = 4.0/pi*cmghz*fx*gam/(SQUARE(f-fx*freq) + SQUARE(2.0*fx*gam));
      // the rest is as for Inversion case:
      alpha_rv += Drv*dens*(x.I0)*fx*beta300_25 * exp(hckt0*(x.E0)*(1.0-beta300))*fline;
   }

   return alpha_inv + alpha_rot + alpha_rv; // in cm^-1
}


double attenuation_NH3_Devaraj(double freq, double P, double P_idl, double T,
       double XH2, double XHe, double XNH3, double XH2O, int version)
   //  sum the contributions of all lines, first the Inversion transitions,
   //  then the Rotational transitions

{
   // TEMPORARY KLUDGE to compare to Hanley, fig19!!!
   //T=295;
   //XNH3=0.004366812227074;

   // versioning
   bool use_v2011(false);
   bool use_high_pressure_parameters(P_idl>= 15 ? true : false);
   bool use_low_pressure_parameters(!use_high_pressure_parameters);

   switch(version) {
       case 0: // use defaults; do not override
          break;
       case -1: // use version 2011
          use_v2011 = true;
          use_high_pressure_parameters = false;
          use_low_pressure_parameters = false;
          break;
       case 1: // use high pressure parameters
          use_high_pressure_parameters = true;
          use_low_pressure_parameters = false;
          break;
       case 2: // use low pressure parameters
          use_high_pressure_parameters = false;
          use_low_pressure_parameters = true;
          break;
   }

   // miscellaneous constants:
   double const ppscale = 2.99792458e18; // converts to cm^2*cm-1 from nm^2*MHz
   double const gcon = 25.923;
   double const cmghz = 29.9792458;    // converts to GHz from cm-1
   double const hckt0 = 0.0047959232;  // hc/kT0, where T0=300K
   double const cdens = 7.242971565e21;  // 1/kB * (1e5 Pa/bar) * (1e-6 m^-3/cm^-3)

   // partial pressures in bars
   double const PNH3 = P_idl* XNH3;
   double const PH2  = P_idl* XH2;
   double const PHe  = P_idl* XHe;
   double const PH2O = P_idl* XH2O;

   double const dens = cdens*(P_idl*XNH3)/T;
   //double const dens = cdens*PNH3/T;
   double const beta300 = 300.0/T;
   double const beta295 = 295.0/T;
   double const beta300_25 = pow(beta300,2.5);

   //  Inversion:
   double alpha_inv(0.0);
   double gh2, zh2, ghe, zhe, gnh3, znh3, gh2o, zh2o, dd, Dinv;
   if (use_v2011) {
      gh2  = 1.7947 * PH2  * pow(beta300,0.8357);
      ghe  = 0.75   * PHe  * pow(beta300,0.6667);
      gnh3 = 0.7719 * PNH3 * pow(beta295,1.0);
      gh2o = 4.8901 * PH2O * pow(beta300,0.5);

      zh2  = 1.2031 * PH2  * pow(beta300,0.8610);
      zhe  = 0.3    * PHe  * pow(beta300,0.6667);
      znh3 = 0.5620 * PNH3 * pow(beta295,0.6206);
      zh2o = 2.731  * PH2O * pow(beta300,1.0);

      dd = -0.0404;      // pressure shift factor from linewidth
      Dinv = 0.9903;      // unitless emperical scale factor
   }
   else {
      if (use_low_pressure_parameters) {
         gh2  = 1.7465 * PH2  * pow(beta300,0.8202);
         ghe  = 0.9779 * PHe  * pow(beta300,1.0);
         gnh3 = 0.7298 * PNH3 * pow(beta295,1.0);
         gh2o = 4.8901 * PH2O * pow(beta300,0.5);

         zh2  = 1.2163 * PH2  * pow(beta300,0.8873);
         zhe  = 0.0291 * PHe  * pow(beta300,0.8994);
         znh3 = 0.5152 * PNH3 * pow(beta295,0.6667);
         zh2o = 2.731  * PH2O * pow(beta300,1.0);

         dd = -0.0627;      // pressure shift factor from linewidth
         Dinv = 0.9862;      // unitless emperical scale factor
      }
      else if (use_high_pressure_parameters) {
         gh2  = 1.6361 * PH2  * pow(beta300,0.8);
         ghe  = 0.4555 * PHe  * pow(beta300,0.5);
         gnh3 = 0.7298 * PNH3 * pow(beta295,1.0);
         gh2o = 4.8901 * PH2O * pow(beta300,0.5);

         zh2  = 1.1313 * PH2  * pow(beta300,0.6234);
         zhe  = 0.1    * PHe  * pow(beta300,0.5);
         znh3 = 0.5152 * PNH3 * pow(beta295,0.6667);
         zh2o = 2.731  * PH2O * pow(beta300,1.0);

         dd = 0.2;      // pressure shift factor from linewidth
         Dinv = 1.3746;      // unitless emperical scale factor
      }
   }

   int const Ninv = sizeof(NH3_inv_JPL5)/sizeof(NH3_inv_JPL5[0]);
   for(int i=0; i<Ninv; i++) {
      LineParameters_NH3_inversion_JPL5 const & x = NH3_inv_JPL5[i];
      double f = x.f0;  // convert to GHz
      double fx = freq/f;
      double sjkt = x.I0 * beta300_25 * exp(hckt0*(x.E0)*(1.0-beta300));

      // pressure-broadening widths:
      double gam = gh2 + ghe + gnh3 * x.gamma_self;
      double zet = zh2 + zhe + znh3 * x.gamma_self;
      if (PH2O != 0.0) { // Broadening line width's due to water Vapor:
         gam += gh2o;
         zet += zh2o;
      }

      // Ben-Reuven line profile (eqn 2.6):
      double num = (gam-zet)*SQUARE(freq) + (gam+zet)*( SQUARE(f+dd*gam) + SQUARE(gam) - SQUARE(zet) );
      double den = SQUARE(SQUARE(freq) - SQUARE(f+dd*gam) - SQUARE(gam) + SQUARE(zet)) +
         SQUARE(2.0*freq*gam);
      double fline = cmghz*2.0*fx*num/(pi*den);

      alpha_inv += Dinv*dens*sjkt*fx*fline;
   }

   //  Rotational:
   double alpha_rot(0.0);
   double const gH2_rot  = 0.2984 * PH2  * pow(beta300,0.8730);
   double const gHe_rot  = 0.75   * PHe  * pow(beta300,0.6667);
   double const gNH3_rot = 3.1789 * PNH3 * pow(beta300,1.0);
   double const Drot = 2.4268;

   int const Nrot = sizeof(NH3_rot_JPL5)/sizeof(NH3_rot_JPL5[0]);
   for(int i=0; i<Nrot; i++) {
      LineParameters_NH3_rotational_JPL5 const & x = NH3_rot_JPL5[i];
      double f = x.f0;  // convert to GHz
      double fx = freq/f;
      // pressure-broadening width:
      double gam = gH2_rot*(x.gamma_H2) + gHe_rot*(x.gamma_He) + gNH3_rot*(x.gamma_self);
      if (PH2O != 0.0) {
         gam += gh2o;
      }
      // Gross line profile (from Joiner at al.):
      double fline = 4.0/pi*cmghz*fx*gam/(SQUARE(f-fx*freq) + SQUARE(2.0*fx*gam));
      // the rest is as for Inversion case:
      alpha_rot += Drot*dens*(x.I0)*fx*beta300_25 * exp(hckt0*(x.E0)*(1.0-beta300))*fline;
   }

   //  Rotovibrational:
   double alpha_rv(0.0);
   double const gH2_nu2  = 1.4  * PH2  * pow(beta300,0.73);
   double const gHe_nu2  = 0.68 * PHe  * pow(beta300,0.5716);
   double const gNH3_nu2 = 9.5  * PNH3 * pow(beta300,1.0);
   double const Drv = 1.1206;

   int const Nrv = sizeof(NH3_rv_JPL5)/sizeof(NH3_rv_JPL5[0]);
   for(int i=0; i<Nrv; i++) {
      LineParameters_NH3_rotovibrational_JPL5 const & x = NH3_rv_JPL5[i];
      double f = x.f0;  // convert to GHz
      double fx = freq/f;
      // pressure-broadening width:
      double gam = gH2_nu2 + gHe_nu2 + gNH3_nu2;
      if (PH2O != 0.0) {
         gam += gh2o;
      }
      // Gross line profile (from Joiner at al.):
      double fline = 4.0/pi*cmghz*fx*gam/(SQUARE(f-fx*freq) + SQUARE(2.0*fx*gam));
      // the rest is as for Inversion case:
      alpha_rv += Drv*dens*(x.I0)*fx*beta300_25 * exp(hckt0*(x.E0)*(1.0-beta300))*fline;
   }

   return alpha_inv + alpha_rot + alpha_rv; // in cm^-1
}


double attenuation_NH3_Hanley(double freq, double P, double P_idl, double T,
       double XH2, double XHe, double XNH3, double XH2O, double power)
{
   // KLUDGE to compare to Hanley, fig19!!!
   //P=P_idl;
   //T=295;
   //XNH3=0.004366812227074;

   static double hif0[NINV], hie0[NINV], hiint[NINV],
      hrf0[NROT], hre0[NROT], hrint[NROT];

   int hij[NINV], hik[NINV], hrj[NROT], hrk[NROT];
   static bool first=true;

   // for line broadening (Table 4.7):
   double const gh2 =1.64;
   double const cgh2 = 0.7756;
   double const zh2 = 1.262;
   double const czh2 = 0.7964;

   double const ghe = 0.75;
   double const cghe = 0.6667;
   double const zhe = 0.3;
   double const czhe = 0.6667;

   double const gnh3 = 0.852;
   double const cgnh3 = 1.0;
   double const znh3 = 0.5296;
   double const cznh3 = 1.554;

   // miscellaneous constants:
   double const ppscale = 3.00e18; // converts to cm^2*cm-1 from nm^2*MHz
   double const gcon = 25.923;
   double const cmghz = 29.979;    // converts to GHz from cm-1

   double const dd = -0.0498;      // pressure shift factor from linewidth
   double const ddd = 0.9301;      // unitless emperical scale factor
   double const hckt0 = 0.004796;  // hc/kT0, where T0=300K
   double const cdens = 7.244e21;  // 1/kB * (1e5 Pa/bar) * (1e-6 m^-3/cm^-3)
   //double const cdens = 7.242971565e21;  // 1/kB * (1e5 Pa/bar) * (1e-6 m^-3/cm^-3)

   // on first call, load the arrays from the tables:
   if (first || 1) {
      for(int i=0; i<NINV; i++) {
	 hif0[i] = 1.0e-3*hireals[i][0];	// convert from MHz to GHz
	 // intensities from Poynter Pickett catalog:
	 hiint[i] = pow(10.0,hireals[i][2])/ppscale;
	 hie0[i] = hireals[i][3];
	 hij[i] = hiints[i][3];
	 hik[i] = hiints[i][4];
      }
      for(int i=0; i<NROT; i++) {
	 hrf0[i] = 1.0e-3*hrreals[i][0];	// convert from MHz to GHz
	 // intensities from Poynter Pickett catalog:
	 hrint[i] = pow(10.0,hrreals[i][2])/ppscale;
	 hre0[i] = hrreals[i][3];
	 hrj[i] = hrints[i][3];
	 hrk[i] = hrints[i][4];
      }
      first = false;
   }

   // partial pressures in bars
   double const PNH3 = P_idl * XNH3;
   double const PH2  = P_idl * XH2;
   double const PHe  = P_idl * XHe;
   double const PH2O = P_idl * XH2O;

   //double const dens = cdens*PNH3/T;
   double const dens = cdens*(P_idl*XNH3)/T;
   double const beta300 = 300.0/T;
   double const beta295 = 295.0/T;

   double const beta300_25 = pow(beta300,2.5);
   double const beta300_cgh2 = pow(beta300,cgh2);
   double const beta300_cghe = pow(beta300,cghe);
   double const beta300_czh2 = pow(beta300,czh2);
   double const beta300_czhe = pow(beta300,czhe);
   double const beta295_cgnh3 = pow(beta295,cgnh3);
   double const beta295_cznh3 = pow(beta295,cznh3);

   // water broadening parameters
   double const gh2o = 4.8901;
   double const cgh2o = 0.5;
   double const zh2o = 2.731;
   double const czh2o = 1;
   double const beta300_cgh2o = pow(beta300,cgh2o);
   double const beta300_czh2o = pow(beta300,czh2o);

   //  sum the contributions of all lines, first innversion then rotational
   double dtau(0.0);
   double dtau_inv(0.0);
   double dtau_rot(0.0);

   //  Inversion:
   for(int i=0; i<NINV; ++i) {
      double fx = freq/hif0[i];
      double sjkt = hiint[i] * beta300_25 * exp(hckt0*hie0[i]*(1.0-beta300));
      // pressure-broadening widths:
      double gam = gh2*PH2 * beta300_cgh2 + ghe * PHe * beta300_cghe;
      double zet = zh2*PH2 * beta300_czh2 + zhe * PHe * beta300_czhe;
      if (PH2O!=0.0) {
         gam += gh2o * PH2O * beta300_cgh2o;
         zet += zh2o * PH2O * beta300_czh2o;
      }

      // self-broadening line width:
      double gam0 = gcon*double(hik[i])/sqrt(double(hij[i]*(hij[i]+1)));
      gam += gnh3 * gam0 * PNH3 * beta295_cgnh3;
      zet += znh3 * gam0 * PNH3 * beta295_cznh3;

      // Ben-Reuven line profile (eqn 2.6):
      double num = (gam-zet)*SQUARE(freq) +
         (gam+zet)*( SQUARE(hif0[i]+dd*gam) + SQUARE(gam) - SQUARE(zet) );
      double den = SQUARE(SQUARE(freq) - SQUARE(hif0[i]+dd*gam) - SQUARE(gam) + SQUARE(zet)) +
	 SQUARE(2.0*freq*gam);
      double fline = cmghz*2.0*fx*num/(pi*den);
      double alpha = ddd*dens*sjkt*fx*fline;
      dtau_inv += alpha;
      dtau += alpha;
   }

   double const beta296 = 296/T;
   double const gH2_rot  = PH2  * pow(beta296,cgh2);
   double const gHe_rot  = PHe  * pow(beta296,cghe);
   double const gNH3_rot = PNH3 * pow(beta296,cgnh3);

   //  Rotational:
   for(int i=0; i<NROT; i++) {
      LineParameters_NH3_rotational_JPL5 const & x = NH3_rot_Hanley[i];
      double fx = freq/hrf0[i];

      double gam = atm_from_bar(gH2_rot*(x.gamma_H2) + gHe_rot*(x.gamma_He) + gNH3_rot*(x.gamma_self));
      if (PH2O!=0.0) {
         gam += gh2o * PH2O * beta300_cgh2o;
      }

      // Gross line profile (from Joiner at al.):
      double fline = 4.0/pi*cmghz*fx*gam/(SQUARE(hrf0[i]-fx*freq) + SQUARE(2.0*fx*gam));
      // the rest is as for Inversion case:
      double dtaua = ddd*dens*hrint[i]*fx*beta300_25 * exp(hckt0*hre0[i]*(1.0-beta300))*fline;
      dtau_rot += dtaua;
      dtau += dtaua;
   }

   #ifdef ABSORPTION_KLUDGE
   constexpr double P_top = 150;
   constexpr double P_bot = 4000;
   constexpr double scale_bot = 2.0;
   if (P>P_bot)
      dtau *= scale_bot;
   else if (P>P_top)
      dtau *= ((P-P_top)*scale_bot + (P_bot-P)) / (P_bot-P_top);
   #endif

   if (T > 750.) dtau *= pow(750./T, power);

   return dtau; // in cm^-1
}

// helper function for attenuation_NH3_radtran() -- array ordering
inline int Index_NH3(int I, int J)
{
   int INDX = J;

   if (I > 1) {
      for (int L=1; L <= I-1; L++) {
         INDX += L;
      }
   }
   return INDX;
}

// helper function for attenuation_NH3_radtran() -- computes Ben-Reuven lineshape
inline double shape(double NU0, double GAMMA, double NU, double DELTA, double ZETA)
{
   double const NUSQ = SQUARE(NU);
   double const NU0SQ = SQUARE(NU0);
   double const GSQ = SQUARE(GAMMA);
   double const ZSQ = SQUARE(ZETA);
   double const T = (NU0 + DELTA)*(NU0 + DELTA);
   return 2.0*GAMMA*(NUSQ/NU0SQ)*
      (((GAMMA-ZETA)*NUSQ + (GAMMA+ZETA)*
        (T+GSQ-ZSQ))/((NUSQ-T-GSQ+ZSQ)*(NUSQ-T-GSQ+ZSQ) + 4.0*NUSQ*GSQ));
}

double attenuation_NH3_radtran(double freq, double P, double T, double XH2,
       double XHe, double XNH3)
{
   /*
      SUBROUTINE NH3ABS CALCULATES THE ABSORPTION DUE TO AMMONIA
      INCLUDING THE PRESSURE BROADENING DUE TO HYDROGEN.

       PNH3 .....      PARTIAL PRESSURE OF NH3 IN ATM
       PH2 ......      PARTIAL PRESSURE OF H2 IN ATM
       T ........      TEMPERATURE IN DEGREES KELVIN
       FQ .......      FREQUENCY IN GHZ AT WHICH ABSORPTION
                       IS TO BE CALCULATED
       DTAU .....      ABSORPTION IN 1/CM AT A PARTICULAR FREQ DUE TO ALL
                       NH3 LINES.
   */

   double DTAU=0.0,DTAUA=0.0;

   double const PNH3 = atm_from_bar(P * XNH3);
   double const PH2  = atm_from_bar(P * XH2);
   double const PHe  = atm_from_bar(P * XHe);

   double const TESTPR=1.0e-4;
   double const EXPON=4.0e0;
   double const GAMFAC=1.0e0;
   double const BOLTZ=2.084e4;
   double const XMU=1.476;
   double const DROT=-1.09107e5;
   double const B=2.98107e5;

   // COLLSW IS A FUDGE FACTOR WHICH INCREASES THE SELF-BROADENING
   // OF NH3 INVERSION LINES ...       ...     FACNH3 = COLLSW*FACNH3
   double const collsw = 1.0;
   double const ZFAC = 1.0;
   double const DFAC = 1.0;

   // COMMON/ATMCA/WT,TAU,DTAU,DTB,freq,ANH3SV,AH2OSV,acldsv

   if (T <= 0.0)  throw "Non-positive temperature";

   double const PCROSS = pow(((PNH3 + PH2/8.0)/TESTPR),EXPON);
   double const SWITCH = ( (PCROSS <= 9.0) ? 1.0 - exp(-PCROSS) : 1.0 );
   double const FACT = pow(T,3.5);
   double const TBOLTZ = T*BOLTZ;
   double const FACH2 = 2.318 * PH2 * pow((300.0/T),(2./3.));
   double const FACHE = (35.403*PHe) / pow(T,(2./3.));
   double const FNH3 = collsw * 224.2 * PNH3 / T;
   double const DELTA = -SWITCH*0.45*PNH3 * DFAC;

   // FACMP IS A FUDGE FACTOR DETERMINED BY COMPARISON OF THIS ROUTINE
   // WITH THE RESULTS OF MORRIS AND PARSONS.
   double const FACMP = 1.00745 + 3.08301e-2*(PH2/T) +
      5.52505e-2*pow((PH2/T),2.0);

   for (int J=1; J<=16; J++) {
      double ENJ = B*J*(J+1);
      double FJMUSQ = (2*J+1)*XMU*XMU/(J*(J+1));

      for (int K=1; K<=J; K++) {
         int I = Index_NH3(J,K) - 1;
         if (FI[I] <= 1.0 || GAMIN[I] <= 1.0) continue;

         double EN = -(ENJ + DROT*K*K)/TBOLTZ;
         double FXMUSQ = FJMUSQ*K*K;
         double FACNH3 = FNH3*GAMIN[I];
         double GAMMA = GAMFAC*(FACNH3 + FACH2 + FACHE);
         double ZETA = SWITCH*(0.655*FACNH3 + 0.83*FACH2 + 0.3797*FACHE) * ZFAC;
         double GK = (K%3==0 ? 3.0 : 1.5);

         double TEXP = 1.132e3*GK*FXMUSQ*FI[I]*FI[I]*PNH3*exp(EN)/(GAMMA*FACT);
         double R = shape(FI[I],GAMMA,freq,DELTA,ZETA);

         // MAKE FURTHER ADJUSTMENT BECAUSE ONLY 1/2 OF MOLECULES ABSORB
         DTAUA = 0.5*FACMP*TEXP*R;
         DTAU += DTAUA;
      }
   }

   // disable this check, as it excludes these transitions for our freqs:
   // if (freq <= FI[6])  return DTAU;

   // INCLUDE THE ABSORPTION DUE TO ROTATIONAL TRANSITIONS
   double ROTFAC = 1.0;
   double FACNH3 = 20.0*FNH3;
   double GAMMA = GAMFAC*(FACNH3 + FACH2 + FACHE);
   double ZETA = SWITCH*(0.655*FACNH3 + 0.83*FACH2 + 0.3797*FACHE) * ZFAC;

   for (int L=0; L<49; L++) {
      int const J = JROT[L];
      int const K = KROT[L];
      double const F = FROT[L];

      double R = shape(F,GAMMA,freq,DELTA,ZETA);
      double FXMUSQ = XMU*XMU*((J+1)*(J+1)-K*K)/(J+1);
      double EN = -(B*J*(J+1) + DROT*K*K)/TBOLTZ;

      double GK = (K%3==0 ? 3.0 : 1.5);
      if (K==0)  GK *= 0.5;

      double TEXP = 1.132e3*GK*FXMUSQ*F*F*PNH3*exp(EN)/(GAMMA*FACT);
      DTAUA = 0.5*ROTFAC*FACMP*TEXP*R;
      DTAU += DTAUA;
   }

   return DTAU;
}
