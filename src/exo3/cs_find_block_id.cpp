// C/C++
#include <iostream>
#include <sstream>

// athena
#include <athena/athena.hpp>

// exo3
#include "cubed_sphere.hpp"

int CubedSphere::FindBlockID(LogicalLocation const& loc) {

#ifdef NBLOCKS
    // Updated method, need to manually setup in configure.hpp, allow 6*n^2 blocks
    int bound_lim = (int)(sqrt(NBLOCKS / 6) + 0.5);
    // Sanity check of NBLOCKS
    static_assert(NBLOCKS == 6 * bound_lim * bound_lim,
                  "NBLOCKS must be 6*(2n)^2 for the cubed sphere.");
#else
    // Old method, suitable for 6*4^n blocks
    int bound_lim = (1 << (loc.level - 2));
#endif
  // Infer the location in the 2*3 configuration
  int lv2_lx2 = loc.lx2 / bound_lim;
  int lv2_lx3 = loc.lx3 / bound_lim;
  // std::cout << "loc.level=" << loc.level << ", loc.lx2=" << loc.lx2 << ", loc.lx3=" << loc.lx3 << std::endl;

  // Determine the block number
  int block_id;
  switch (lv2_lx3) {
    case 0:
      switch (lv2_lx2) {
        case 0:
          block_id = 1;
          break;
        case 1:
          block_id = 2;
          break;
        default:
          std::stringstream msg;
          msg << "Error: something wrong, check the geometry setup of the "
                 "cubed sphere. \n";
          msg << "----------------------------------" << std::endl;
          ATHENA_ERROR(msg);
      }
      break;
    case 1:
      switch (lv2_lx2) {
        case 0:
          block_id = 3;
          break;
        case 1:
          block_id = 4;
          break;
        default:
          std::stringstream msg;
          msg << "Error: something wrong, check the geometry setup of the "
                 "cubed sphere. \n";
          msg << "----------------------------------" << std::endl;
          ATHENA_ERROR(msg);
      }
      break;
    case 2:
      switch (lv2_lx2) {
        case 0:
          block_id = 5;
          break;
        case 1:
          block_id = 6;
          break;
        default:
          std::stringstream msg;
          msg << "Error: something wrong, check the geometry setup of the "
                 "cubed sphere. \n";
          msg << "----------------------------------" << std::endl;
          ATHENA_ERROR(msg);
      }
      break;
    default:
      std::stringstream msg;
      msg << "Error: something wrong, check the geometry setup of the cubed "
             "sphere. \n";
      msg << "----------------------------------" << std::endl;
      ATHENA_ERROR(msg);
  }

  return block_id;
}
