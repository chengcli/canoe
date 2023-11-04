// C/C++
#include <functional>
#include <string>

// athena++
#include <athena/bvals/bvals.hpp>

namespace MessageHelper {
int mpi_tag_ub;

int create_mpi_tag(int lid, int tid, std::string name) {
  int tag = BoundaryBase::CreateBvalsMPITag(lid, tid, 0);

  std::string str = name + std::to_string(tag);
  return std::hash<std::string>{}(str) % (mpi_tag_ub);
}
}  // namespace MessageHelper
