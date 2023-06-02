// C/C++ headers
#include <glob.h>
#include <stdio.h>

#include <cstring>
#include <iostream>
#include <sstream>  // stringstream
#include <stdexcept>

// athena
#include <athena/globals.hpp>

// canoe
#include <configure.hpp>

// outputs
#include "user_outputs.hpp"

int mppnccombine(int argc, char *argv[]);

void NetcdfOutput::CombineBlocks() {
// Only proceed if NETCDF output enabled
#ifdef NETCDFOUTPUT

  std::stringstream msg;
#ifdef MPI_PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  if (Globals::my_rank == 0) {
    char number[64];
    snprintf(number, sizeof(number), "%05d", output_params.file_number - 1);

    std::string infile;
    infile.assign(output_params.file_basename);
    infile.append(".block*.");
    infile.append(output_params.file_id);
    infile.append(".");
    infile.append(number);
    infile.append(".nc");

    glob_t glob_result;
    int err = glob(infile.c_str(), GLOB_TILDE, NULL, &glob_result);
    if (err != 0) {
      globfree(&glob_result);
      msg << "### FATAL ERROR in function [NetcdfOutput::CombineBlocks]"
          << std::endl
          << "glob() failed with error " << err << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }

    std::string outfile;
    outfile.assign(output_params.file_basename);
    outfile.append(".");
    outfile.append(output_params.file_id);
    outfile.append(".");
    outfile.append(number);
    outfile.append(".nc");

    int argc = 3 + glob_result.gl_pathc;
    // char argv[][2048] = {"CombineBlocks", "-r", outfile.c_str(),
    // infile.c_str()};
    char **argv = new char *[argc];
    for (int i = 0; i < argc; ++i) argv[i] = new char[2048];
    snprintf(argv[0], sizeof(argv[0]), "%s", "CombineBlocks");
    snprintf(argv[1], sizeof(argv[1]), "%s", "-r");
    snprintf(argv[2], sizeof(argv[2]), "%s", outfile.c_str());
    for (int i = 3; i < argc; ++i)
      snprintf(argv[i], sizeof(argv[i]), "%s", glob_result.gl_pathv[i - 3]);

    remove(outfile.c_str());
    err = mppnccombine(argc, argv);
    if (err) {
      std::cerr << "### WARNING in function [NetcdfOutput::CombineBlocks]"
                << std::endl
                << "mppnccombine returns none zero.";
      // throw std::runtime_error(msg.str().c_str());
    }

    globfree(&glob_result);
    for (int i = 0; i < argc; ++i) delete[] argv[i];
    delete[] argv;
  }

#endif  // NETCDFOUTPUT
}
