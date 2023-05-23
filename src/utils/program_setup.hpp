#ifndef PROGRAM_SETUP_HPP
#define PROGRAM_SETUP_HPP

// C/C++ headers
#include <cstdint>

// forward declaration
class Mesh;
class ParameterInput;
struct CommandLine;

namespace Globals {
extern int mpi_tag_ub;
extern clock_t tstart;
extern std::uint64_t mbcnt;
extern CommandLine *cli;
}  // namespace Globals

void program_start(int, char **);

void program_end();

void program_end(Mesh *);

#endif
