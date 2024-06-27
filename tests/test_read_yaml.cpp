// C/C++
#include <fstream>
#include <vector>

// external
#include <gtest/gtest.h>
#include <yaml-cpp/yaml.h>

class YamlReadTests : public ::testing::Test {
 protected:
  std::string filename = "example_infrared.yaml";  // The name of your YAML file
  YAML::Node node;

  virtual void SetUp() {
    std::ifstream stream(filename);
    ASSERT_TRUE(stream.good()) << "Failed to open file: " << filename;
    ASSERT_NO_THROW(node = YAML::Load(stream)) << "Failed to parse YAML file";
  }

  virtual void TearDown() {
    // You can perform teardown actions here
  }
};

TEST_F(YamlReadTests, OpacitySource) {
  // Checking 'opacity-source'
  ASSERT_TRUE(node["opacity-sources"])
      << "'opacity-sources' not found in YAML file";

  auto sources = node["opacity-sources"].as<std::vector<YAML::Node>>();
  ASSERT_EQ(6, sources.size())
      << "Unexpected number of sources in 'opacity-sources'";

  // Check first source
  EXPECT_EQ("H2-H2-CIA", sources[0]["name"].as<std::string>());
  EXPECT_EQ("Hydrogen-Hydrogen collisional absorption",
            sources[0]["long-name"].as<std::string>());

  // Check third source
  EXPECT_EQ("CH4", sources[2]["name"].as<std::string>());
  EXPECT_EQ("Methane line absorption",
            sources[2]["long-name"].as<std::string>());
  EXPECT_EQ("Voigt", sources[2]["model"].as<std::string>());
}

TEST_F(YamlReadTests, BandsList) {
  // Checking 'bands'
  ASSERT_TRUE(node["bands"]) << "'bands' not found in YAML file";
  std::vector<std::string> bands = node["bands"].as<std::vector<std::string>>();
  ASSERT_EQ(2, bands.size()) << "Unexpected number of bands";
}

TEST_F(YamlReadTests, RadiationBand) {
  // Checking 'ir'
  ASSERT_TRUE(node["ir"]) << "'ir' not found in YAML file";
  EXPECT_EQ((std::vector<double>{10., 200.}),
            node["ir"]["wavenumber-range"].as<std::vector<double>>());
  EXPECT_EQ(0.01, node["ir"]["resolution"].as<double>());
  EXPECT_EQ((std::vector<std::string>{"H2-H2-CIA", "H2-He-CIA", "CH4", "C2H2",
                                      "C2H4", "C2H6"}),
            node["ir"]["opacity"].as<std::vector<std::string>>());
  EXPECT_EQ(true, node["ir"]["heating-flux"].as<bool>());
  EXPECT_EQ(false, node["ir"]["spectrum"].as<bool>());
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
