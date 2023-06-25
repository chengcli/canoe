// external
#include <gtest/gtest.h>
#include <yaml-cpp/yaml.h>

// harp
#include <harp/absorber.hpp>

TEST(TestAbsorber, Construct) {
  Absorber ab("dummy");

  EXPECT_EQ(ab.GetName(), "dummy");
};

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  auto app = Application::GetInstance();

  // app->InstallMonitor("harp", "harp.out", "harp.err");

  int result = RUN_ALL_TESTS();

  Application::Destroy();

  return result;
}
