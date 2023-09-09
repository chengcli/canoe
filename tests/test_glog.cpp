// C/C++
#include <algorithm>
#include <cmath>

// external
#include <gtest/gtest.h>

// athena
#include <athena/reconstruct/interpolation.hpp>

// canoe
#include <configure.hpp>

#ifdef ENABLE_GLOG
#include <glog/logging.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

TEST(glob, test_log_info) {
  LOG(INFO) << "test INFO";
  LOG(WARNING) << "test WARNING";
  LOG(ERROR) << "test ERROR";
}

#endif

/*TEST(interp_weno5, test_case2) {
  double phim2 = 1.0;
  double phim1 = 2.0;
  double phi = 3.0;
  double phip1 = 4.0;
  double phip2 = 5.0;
  double result = interp_weno5(phim2, phim1, phi, phip1, phip2);
  double expected_result = 2.5000000000000004;
  EXPECT_NEAR(result, expected_result, 1.E-10);
}*/

int main(int argc, char **argv) {
#ifdef ENABLE_GLOG
  std::string prog_name = argv[0];
  // Extract the base name
  size_t last_slash_pos = prog_name.find_last_of("/");
  std::string base_name = (last_slash_pos == std::string::npos)
                              ? prog_name
                              : prog_name.substr(last_slash_pos + 1);
  std::string log_dir_name = base_name + ".glog";

  struct stat st = {0};
  if (stat(log_dir_name.c_str(), &st) ==
      -1) {  // Check if the directory doesn't exist.
    if (mkdir(log_dir_name.c_str(), 0755) ==
        -1) {  // Mode 0755 gives rwx permissions for everyone.
      perror("Error creating directory");
    } else {
      std::cout << "Directory " << log_dir_name << " created." << std::endl;
    }
  } else {
    std::cout << "Directory " << log_dir_name << " already exists."
              << std::endl;
  }

  google::InitGoogleLogging(argv[0]);
  FLAGS_log_dir = log_dir_name;
  int my_rank = 0;
  char log_name[FILENAME_MAX];
  snprintf(log_name, FILENAME_MAX, "%s.%d", base_name.c_str(), my_rank);
  google::SetLogSymlink(google::GLOG_INFO, log_name);
  google::SetLogSymlink(google::GLOG_WARNING, log_name);
  google::SetLogSymlink(google::GLOG_ERROR, log_name);
#endif

  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
