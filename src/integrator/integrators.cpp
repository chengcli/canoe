//! \file  integrators.cpp
//! \brief Contains the implementation of the integrators.

void MultiStageIntegrator::SetIntegrator(std::string integrator) {
  if (integrator == "vl2") {
    //! \note `integrator == "vl2"`
    //! - VL: second-order van Leer integrator (Stone & Gardiner, NewA 14, 139
    //! 2009)
    //! - Simple predictor-corrector scheme similar to MUSCL-Hancock
    //! - Expressed in 2S or 3S* algorithm form

    // set number of stages and time coeff.
    nstages_ = 2;
  } else if (integrator == "rk3") {
    //! \note `integrator == "rk3"`
    //! - SSPRK (3,3): Gottlieb (2009) equation 3.2
    //! - Optimal (in error bounds) explicit three-stage, third-order SSPRK

    // set number of stages and time coeff.
    nstages_ = 3;
    stage_wghts_[0].sbeta = 0.0;
    stage_wghts_[0].ebeta = 1.0;
    stage_wghts_[0].delta = 1.0;
    stage_wghts_[0].gamma_1 = 0.0;
    stage_wghts_[0].gamma_2 = 1.0;
    stage_wghts_[0].gamma_3 = 0.0;

    stage_wghts_[1].sbeta = 1.0;
    stage_wghts_[1].ebeta = 0.5;
    stage_wghts_[1].delta = 0.0;
    stage_wghts_[1].gamma_1 = 0.25;
    stage_wghts_[1].gamma_2 = 0.75;
    stage_wghts_[1].gamma_3 = 0.0;

    stage_wghts_[2].sbeta = 0.5;
    stage_wghts_[2].ebeta = 1.0;
    stage_wghts_[2].delta = 0.0;
    stage_wghts_[2].gamma_1 = 2. / 3.;
    stage_wghts_[2].gamma_2 = 1. / 3.;
    stage_wghts_[2].gamma_3 = 0.0;

    stage_wghts_[0].beta = 1.0;
    stage_wghts_[1].beta = 0.25;
    stage_wghts_[2].beta = 2. / 3.;
  };
