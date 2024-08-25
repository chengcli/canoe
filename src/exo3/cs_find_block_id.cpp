// C/C++
#include <iostream>
#include <sstream>

// athena
#include <athena/athena.hpp>

// exo3
#include "cubed_sphere.hpp"

int CubedSphere::FindBlockID(LogicalLocation const& loc) {
  int lv2_lx2 = loc.lx2 >> (loc.level - 2);
  int lv2_lx3 = loc.lx3 >> (loc.level - 2);
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
