// C/C++ headers
#include <iostream>

// Athena++ headers
#include <globals.hpp>

// debugger
#include <debugger.hpp>

// cliutils headers
#include <configure.hpp>
#include "command_line.hpp"

#ifdef MPI_PARALLEL
  #include <mpi.h>
#endif

char* CommandLine::input_filename = nullptr;
char* CommandLine::restart_filename = nullptr;
char* CommandLine::prundir = nullptr;
int CommandLine::res_flag = 0;
int CommandLine::narg_flag = 0;
int CommandLine::iarg_flag = 0;
int CommandLine::mesh_flag = 0;
int CommandLine::wtlim = 0;
int CommandLine::argc = 0;
char** CommandLine::argv = nullptr;

CommandLine::CommandLine(int argc, char **argv)
{
  CommandLine::argc = argc;
  CommandLine::argv = argv;

  for (int i=1; i<argc; i++) {
    // If argv[i] is a 2 character string of the form "-?" then:
    if (*argv[i] == '-'  && *(argv[i]+1) != '\0' && *(argv[i]+2) == '\0') {
      // check validity of command line options + arguments:
      char opt_letter = *(argv[i]+1);
      switch(opt_letter) {
        // options that do not take arguments:
        case 'n':
        case 'c':
        case 'h':
          break;
          // options that require arguments:
        default:
          if ((i+1 >= argc) // flag is at the end of the command line options
              || (*argv[i+1] == '-') ) { // flag is followed by another flag
            if (Globals::my_rank == 0) {
              Debugger::Fatal("main",
                "-" + std::to_string(opt_letter) 
                + " must be followed by a valid argument");
            }
          }
      }
      switch(*(argv[i]+1)) {
        case 'i':                      // -i <input_filename>
          input_filename = argv[++i];
          iarg_flag = 1;
          break;
        case 'r':                      // -r <restart_file>
          res_flag = 1;
          restart_filename = argv[++i];
          break;
        case 'd':                      // -d <run_directory>
          prundir = argv[++i];
          break;
        case 'n':
          narg_flag = 1;
          break;
        case 'm':                      // -m <nproc>
          mesh_flag = static_cast<int>(std::strtol(argv[++i], nullptr, 10));
          break;
        case 't':                      // -t <hh:mm:ss>
          int wth, wtm, wts;
          std::sscanf(argv[++i], "%d:%d:%d", &wth, &wtm, &wts);
          wtlim = wth*3600 + wtm*60 + wts;
          break;
        case 'c':
          //if (Globals::my_rank == 0) ShowConfig();
#ifdef MPI_PARALLEL
          MPI_Finalize();
#endif
          return;
          break;
        case 'h':
        default:
          if (Globals::my_rank == 0) {
            std::cout << "HARP 2.0.1" << std::endl;
            std::cout << "Usage: " << argv[0] << " [options] [block/par=value ...]\n";
            std::cout << "Options:" << std::endl;
            std::cout << "  -i <file>       specify input file [athinput]\n";
            std::cout << "  -r <file>       restart with this file\n";
            std::cout << "  -d <directory>  specify run dir [current dir]\n";
            std::cout << "  -n              parse input file and quit\n";
            std::cout << "  -c              show configuration and quit\n";
            std::cout << "  -m <nproc>      output mesh structure and quit\n";
            std::cout << "  -t hh:mm:ss     wall time limit for final output\n";
            std::cout << "  -h              this help\n";
            //ShowConfig();
          }
#ifdef MPI_PARALLEL
          MPI_Finalize();
#endif
          return;
          break;
      }
    } // else if argv[i] not of form "-?" ignore it here (tested in ModifyFromCmdline)
  }

  return;
}
