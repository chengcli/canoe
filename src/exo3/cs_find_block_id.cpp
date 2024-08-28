// C/C++
#include <iostream>
#include <sstream>

// athena
#include <athena/athena.hpp>

// exo3
#include "cubed_sphere.hpp"

int CubedSphere::FindBlockID(LogicalLocation const& loc) {

  int lv2_lx2, lv2_lx3, local_lx2, local_lx3, bound_lim;
  GetLocalIndex(&lv2_lx2, &lv2_lx3, &local_lx2, &local_lx3, &bound_lim, loc);

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
          msg << "lv2_lx2: " << lv2_lx2 << " lv2_lx3: " << lv2_lx3 << "bound_lim: " << bound_lim << "loc.lx2: " << loc.lx2 << "loc.lx3: " << loc.lx3 << std::endl;
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

void CubedSphere::GetLocalIndex(int *lv2_lx2, int *lv2_lx3, int *local_lx2, int *local_lx3, int *bound_lim, LogicalLocation const& loc) {
#ifdef USE_NBLOCKS
  // Updated method, need to manually setup in configure.hpp, allow 6*n^2 blocks
  bound_lim = (int)(sqrt(NBLOCKS / 6) - 0.5);
  // Find relative location within block
  lv2_lx2 = loc.lx2 / (bound_lim + 1);
  lv2_lx3 = loc.lx3 / (bound_lim + 1);
  local_lx2 = loc.lx2 - (lv2_lx2 * (bound_lim + 1));
  local_lx3 = loc.lx3 - (lv2_lx3 * (bound_lim + 1));
#else
  // Old method, suitable for 6*4^n blocks
  bound_lim = (1 << (loc.level - 2)) - 1;
    // Find relative location within block
  lv2_lx2 = loc.lx2 >> (loc.level - 2);
  lv2_lx3 = loc.lx3 >> (loc.level - 2);
  local_lx2 = loc.lx2 - (lv2_lx2 << (loc.level - 2));
  local_lx3 = loc.lx3 - (lv2_lx3 << (loc.level - 2));
#endif
}