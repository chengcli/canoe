#pragma once

// C/C++
#include <future>

// torch
#include <torch/nn/cloneable.h>
#include <torch/nn/module.h>
#include <torch/nn/modules/common.h>

// share
// clang-format off
#include <configure.hpp>
#include <add_arg.h>
// clang-format on

using SharedData = std::shared_ptr<
    std::unordered_map<std::string, std::shared_future<torch::Tensor>>>;

struct SedimentationOptions {
  //! radius and density of particles
  ADD_ARG(std::vector<double>, radius) = {10.0e-6};
  ADD_ARG(std::vector<double>, density) = {1.0e3};

  //! additional constant sedimentation velocity
  ADD_ARG(std::vector<double>, const_vsed) = {};

  ADD_ARG(double, gravity) = -10.0;

  //! default H2-atmosphere properties
  //! diameter of molecule [m]
  ADD_ARG(double, a_diameter) = 2.827e-10;

  //! Lennard-Jones potential in J [J]
  ADD_ARG(double, a_epsilon_LJ) = 59.7e-7;

  //! molecular mass of H2 [kg]
  ADD_ARG(double, a_mass) = 3.34e-27;

  //! minimum radius of particles subject to sedimentation [m]
  ADD_ARG(double, min_radius) = 1.e-6;

  //! upper limit of sedimentation velocity [m/s]
  ADD_ARG(double, upper_limit) = 5.e3;
};

class SedimentationImpl : public torch::nn::Cloneable<SedimentationImpl> {
 public:
  //! particle radius and density
  //! 1D tensor of number of particles
  //! radius and density must have the same size
  torch::Tensor radius, density;

  //! options with which this `Sedimentation` was constructed
  SedimentationOptions options;

  //! Constructor to initialize the layers
  explicit SedimentationImpl(SedimentationOptions const& options_)
      : options(options_) {
    reset();
  }
  void reset() override;

  //! Calculate sedimentation velocites
  /*!
   * Calling this function requires the shared future diagnostic data
   * `temperature` to be set. Otherwise, the function will throw an error.
   *
   * \param hydro_w 4D tensor of hydro variables
   * \return 4D tensor of sedimentation velocities. The first dimension is the
   *         number of particles.
   */
  torch::Tensor forward(torch::Tensor hydro_w);

  //! Set shared future diagnostic data
  /*!
   * \param data shared future diagnostic data
   */
  void set_shared_data(SharedData data) { shared_ = data; }

 private:
  SharedData shared_;
};
TORCH_MODULE(Sedimentation);
