// external
#include <gtest/gtest.h>

// opacity
#include <opacity/microwave/mwr_absorbers.hpp>

class MicrowaveOpacityCIA : public ::testing::Test {
 protected:
  MwrAbsorberCIA* mwr_cia_;
  /*MwrAbsorberNH3 mwr_nh3_;
  MwrAbsorberH2O mwr_h2o_;
  MwrAbsorberPH3 mwr_ph3_;
  MwrAbsorberH2S mwr_h2s_;
  MwrAbsorberElectron mwr_o2_;*/

  virtual void SetUp() {
    // Set the temperature and pressure
    mwr_cia_.set_temperature(300.0);
    mwr_cia_.set_pressure(1.0);

    // Set the frequency grid
    mwr_cia_.set_frequency_grid(1.0, 1000.0, 1000);
  }
};

int main(int argc, char* argv[]) {
  testing::InitGoogleTest(&argc, argv);

  // Create the microwave opacity calculator
  MicrowaveOpacityCalculator opacity_calculator;

  // Set the gas composition
  opacity_calculator.set_gas_composition(
      "H2O:0.1, CO2:0.1, CO:0.1, CH4:0.1, N2O:0.1, O3:0.1, O2:0.1, NO:0.1, "
      "NO2:0.1, NH3:0.1");

  // Set the temperature and pressure
  opacity_calculator.set_temperature(300.0);
  opacity_calculator.set_pressure(1.0);

  // Set the frequency grid
  opacity_calculator.set_frequency_grid(1.0, 1000.0, 1000);

  // Calculate the opacity
  opacity_calculator.calculate_opacity();

  // Write the opacity to a file
  opacity_calculator.write_opacity("test_microwave_opacity.dat");

  return RUN_ALL_TESTS();
}
