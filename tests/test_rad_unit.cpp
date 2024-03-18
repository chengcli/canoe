#include <iostream>
// harp
#include <fstream>
#include <harp/radiation.hpp>
#include <index_map.hpp>
#include <opacity/hitran_absorber.hpp>
#include <vector>
// #include <virtual_groups.hpp>
// #include "rt_solvers.hpp"

void calc_and_save_fluxes(RadiationBand& newRadBand, std::string name_band,
                          IndexMap const* index_map);

int main() {
  YAML::Node yfile = YAML::LoadFile("amars.yaml");

  RadiationBand b1("b1", yfile);
  RadiationBand b2("b2", yfile);
  RadiationBand b3("b3", yfile);
  RadiationBand b4("b4", yfile);
  RadiationBand b5("b5", yfile);
  RadiationBand b6("b6", yfile);
  RadiationBand b7("b7", yfile);
  RadiationBand b8("b8", yfile);

  std::map<std::string, std::vector<std::string>> species_map = {
      {"vapor", {"H2S", "SO2", "H2O", "CO2"}}};
  IndexMap const* the_index_map = IndexMap::InitFromSpeciesMap(species_map);

  calc_and_save_fluxes(b1, "b1", the_index_map);
  calc_and_save_fluxes(b2, "b2", the_index_map);
  calc_and_save_fluxes(b3, "b3", the_index_map);
  calc_and_save_fluxes(b4, "b4", the_index_map);
  calc_and_save_fluxes(b5, "b5", the_index_map);
  calc_and_save_fluxes(b6, "b6", the_index_map);
  calc_and_save_fluxes(b7, "b7", the_index_map);
  calc_and_save_fluxes(b8, "b8", the_index_map);

  IndexMap::Destroy();

  return 0;
}

void calc_and_save_fluxes(RadiationBand& newRadBand, std::string name_band,
                          IndexMap const* index_map) {
  auto co2_abs = newRadBand.GetAbsorberByName("CO2");
  co2_abs->LoadCoefficient("kcoeff-B1.nc", 0);

  auto h2s_abs = newRadBand.GetAbsorberByName("H2S");
  h2s_abs->LoadCoefficient("kcoeff-B1.nc", 0);

  auto so2_abs = newRadBand.GetAbsorberByName("SO2");
  so2_abs->LoadCoefficient("kcoeff-B1.nc", 0);

  auto h2o_abs = newRadBand.GetAbsorberByName("H2O");
  h2o_abs->LoadCoefficient("kcoeff-B1.nc", 0);

  AirParcel newParcel;

  // co2_abs -> SetSpeciesIndex(species_map["vapor"]);
  co2_abs->SetSpeciesIndex({"vapor.CO2"});
  h2s_abs->SetSpeciesIndex({"vapor.H2S"});
  so2_abs->SetSpeciesIndex({"vapor.SO2"});
  h2o_abs->SetSpeciesIndex({"vapor.H2O"});

  int n_layers = 40;

  std::string directory =
      "/home/cometz/Desktop/rce/ancientMars_so2Cond_andwater/profiles/";
  std::vector<std::string> filenames = {
      directory + "CO2ppmVals.txt", directory + "H2OppmVals.txt",
      directory + "H2SppmVals.txt", directory + "hVals.txt",
      directory + "pVals.txt",      directory + "SO2ppmVals.txt",
      directory + "TVals.txt"};

  std::vector<double> CO2ppmVals(n_layers), H2OppmVals(n_layers),
      H2SppmVals(n_layers), hVals(n_layers), pVals(n_layers),
      SO2ppmVals(n_layers), TVals(n_layers);

  std::vector<std::vector<double>*> data = {
      &CO2ppmVals, &H2OppmVals, &H2SppmVals, &hVals,
      &pVals,      &SO2ppmVals, &TVals};

  for (int i = 0; i < filenames.size(); ++i) {
    std::ifstream file(filenames[i]);
    for (int j = 0; j < n_layers; ++j) {
      file >> (*data[i])[j];
    }
  }

  AirColumn ac(n_layers);
  for (int i = 0; i < n_layers; ++i) {
    newParcel.w[IDN] = TVals[i];        // kelvin
    newParcel.w[IPR] = pVals[i] * 100;  // Pa
    newParcel.w[index_map->GetVaporId("H2S")] =
        H2SppmVals[i] / 1e6;  // mixing ratio
    newParcel.w[index_map->GetVaporId("SO2")] = SO2ppmVals[i] / 1e6;
    newParcel.w[index_map->GetVaporId("H2O")] = H2OppmVals[i] / 1e6;
    newParcel.w[index_map->GetVaporId("CO2")] = CO2ppmVals[i] / 1e6;
    ac[i] = newParcel;

    hVals[i] = hVals[i] * 1000;  // km -> m
  }

  ac.pop_back();  // remove top layer CURRENTLY, ALL VALUES ARE GIVEN AT EACH
                  // HVAL LINE IN THE ATMOSPHERE THIS LINE OF CODE PROJECTS THE
                  // VALUE AT A GIVEN LINE TO THE LAYER ABOVE IT, UP UNTIL THE
                  // NEXT LINE REALLY, EACH AIRPARCEL OBJECT IS A LAYER

  // for (int i = 0; i < n_layers; ++i)
  //{
  //     std::cout << ac[i] << std::endl;
  // }

  Real* x1f = &hVals[0];
  newRadBand.Resize(ac.size(), 1, 1, 4);
  newRadBand.SetSpectralProperties(ac, x1f, 0, 0, 0);

  newRadBand.CalBandFlux(nullptr, 0, 0, 0, ac.size());

  // Open a file output stream
  std::ofstream outFile("output_" + name_band + ".txt");

  // Write the bfluxup and bflxdn arrays to the file in two columns
  for (int i = 0; i < newRadBand.bflxup.GetDim1(); ++i) {
    outFile << newRadBand.bflxup(i) << " " << newRadBand.bflxdn(i) << " "
            << newRadBand.btoa(i) << " " << newRadBand.btau(i) << '\n';
  }

  // Close the file
  outFile.close();

  // auto parcel = ac[0];

  // Open a file output stream
  // std::ofstream outFile2("attenuation_coefficients.txt");

  // double attenuationco2;
  // double attenuationso2;
  // double attenuationh2o;
  // double attenuationh2s;

  // Loop over wav1 from 1 to 3000
  // for (int wav1 = 1; wav1 <= 3000; ++wav1)
  //{
  // Get the attenuation value
  // attenuationco2 = co2_abs->GetAttenuation(wav1, 5, parcel);
  // attenuationso2 = so2_abs->GetAttenuation(wav1, 5, parcel);
  // attenuationh2o = h2o_abs->GetAttenuation(wav1, 5, parcel);
  // attenuationh2s = h2s_abs->GetAttenuation(wav1, 5, parcel);

  // outFile2 << attenuationco2 << " " << attenuationso2 << " " <<
  // attenuationh2o << " " << attenuationh2s << '\n';
  //}
  // outFile2.close();

  // std::cout << ac[0].w[mySpeciesId(0)] << std::endl;

  // newParcel.w[index_map -> GetVaporId("H2S")] = 0.1;
  // newParcel.w[index_map -> GetVaporId("SO2")] = 0.1;
  // newParcel.w[index_map -> GetVaporId("H2O")] = 0.1;
  // newParcel.w[index_map -> GetVaporId("CO2")] = 0.7;

  // std::cout << index_map -> GetVaporId("H2S") << std::endl;
  // std::cout << index_map -> GetVaporId("SO2") << std::endl;
  // std::cout << index_map -> GetVaporId("H2O") << std::endl;
  // std::cout << index_map -> GetVaporId("CO2") << std::endl;

  // std::cout << newParcel << std::endl;
  // std::cout << "IPR: " << IPR << std::endl;
  // std::cout << "IDN: " << IDN << std::endl;

  // IndexMap::Destroy();
}
