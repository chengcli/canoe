// C/C++
#include <iostream>
#include <sstream>

// athena
#include <athena/athena.hpp>

// exo3
#include "cubed_sphere.hpp"

void CubedSphere::TransformOX(int *ox2, int *ox3, int *tox2, int *tox3,
                              LogicalLocation const &loc) {
  // Find the block ID
  int block_id = FindPanelID(loc);

  // Find relative location within block
  int local_lx2, local_lx3, bound_lim;
  FindPanelPosition(&local_lx2, &local_lx3, &bound_lim, loc);

  // Hard code the cases...
  // No need to consider the corner cases, abandon in reading buffers.
  int target_block = -1;           // Block id of target
  int target_loc_2, target_loc_3;  // local x2 and x3 in target block

  switch (block_id) {
    case 1:                                    // Special: left, top, right
      if ((local_lx3 == 0) && (*ox3 == -1)) {  // Left Boundary
        target_block = 3;
        target_loc_2 = 0;
        target_loc_3 = local_lx2;
        *tox2 = -1;
        *tox3 = 0;
      }
      if ((local_lx2 == 0) && (*ox2 == -1)) {  // Top Boundary
        target_block = 6;
        target_loc_2 = 0;
        target_loc_3 = bound_lim - local_lx3;
        *tox2 = -1;
        *tox3 = 0;
      }
      if ((local_lx3 == bound_lim) && (*ox3 == 1)) {  // Right Boundary
        target_block = 4;
        target_loc_2 = 0;
        target_loc_3 = bound_lim - local_lx2;
        *tox2 = -1;
        *tox3 = 0;
      }
      break;
    case 2:
      if ((local_lx3 == 0) && (*ox3 == -1)) {  // Left Boundary
        target_block = 3;
        target_loc_2 = local_lx2;
        target_loc_3 = bound_lim;
        *tox2 = 0;
        *tox3 = 1;
      }
      if ((local_lx2 == bound_lim) && (*ox2 == 1)) {  // Bottom Boundary
        target_block = 5;
        target_loc_2 = 0;
        target_loc_3 = local_lx3;
        *tox2 = -1;
        *tox3 = 0;
      }
      break;
    case 3:
      if ((local_lx3 == 0) && (*ox3 == -1)) {  // Left Boundary
        target_block = 6;
        target_loc_2 = local_lx2;
        target_loc_3 = bound_lim;
        *tox2 = 0;
        *tox3 = 1;
      }
      if ((local_lx2 == 0) && (*ox2 == -1)) {  // Top Boundary
        target_block = 1;
        target_loc_2 = local_lx3;
        target_loc_3 = 0;
        *tox2 = 0;
        *tox3 = -1;
      }
      if ((local_lx2 == bound_lim) && (*ox2 == 1)) {  // Bottom Boundary
        target_block = 5;
        target_loc_2 = bound_lim - local_lx3;
        target_loc_3 = 0;
        *tox2 = 0;
        *tox3 = -1;
      }
      if ((local_lx3 == bound_lim) && (*ox3 == 1)) {  // Right Boundary
        target_block = 2;
        target_loc_2 = local_lx2;
        target_loc_3 = 0;
        *tox2 = 0;
        *tox3 = -1;
      }
      break;
    case 4:
      if ((local_lx2 == 0) && (*ox2 == -1)) {  // Top Boundary
        target_block = 1;
        target_loc_2 = bound_lim - local_lx3;
        target_loc_3 = bound_lim;
        *tox2 = 0;
        *tox3 = 1;
      }
      if ((local_lx2 == bound_lim) && (*ox2 == 1)) {  // Bottom Boundary
        target_block = 5;
        target_loc_2 = local_lx3;
        target_loc_3 = bound_lim;
        *tox2 = 0;
        *tox3 = 1;
      }
      break;
    case 5:
      if ((local_lx3 == 0) && (*ox3 == -1)) {  // Left Boundary
        target_block = 3;
        target_loc_2 = bound_lim;
        target_loc_3 = bound_lim - local_lx2;
        *tox2 = 1;
        *tox3 = 0;
      }
      if ((local_lx3 == bound_lim) && (*ox3 == 1)) {  // Right Boundary
        target_block = 4;
        target_loc_2 = bound_lim;
        target_loc_3 = local_lx2;
        *tox2 = 1;
        *tox3 = 0;
      }
      if ((local_lx2 == bound_lim) && (*ox2 == 1)) {  // Bottom Boundary
        target_block = 6;
        target_loc_2 = bound_lim;
        target_loc_3 = bound_lim - local_lx3;
        *tox2 = 1;
        *tox3 = 0;
      }
      if ((local_lx2 == 0) && (*ox2 == -1)) {  // Top Boundary
        target_block = 2;
        target_loc_2 = bound_lim;
        target_loc_3 = local_lx3;
        *tox2 = 1;
        *tox3 = 0;
      }
      break;
    case 6:
      if ((local_lx2 == 0) && (*ox2 == -1)) {  // Top Boundary
        target_block = 1;
        target_loc_2 = 0;
        target_loc_3 = bound_lim - local_lx3;
        *tox2 = -1;
        *tox3 = 0;
      }
      if ((local_lx2 == bound_lim) && (*ox2 == 1)) {  // Bottom Boundary
        target_block = 5;
        target_loc_2 = bound_lim;
        target_loc_3 = bound_lim - local_lx3;
        *tox2 = 1;
        *tox3 = 0;
      }
      if ((local_lx3 == bound_lim) && (*ox3 == 1)) {  // Right Boundary
        target_block = 3;
        target_loc_2 = local_lx2;
        target_loc_3 = 0;
        *tox2 = 0;
        *tox3 = -1;
      }
      break;
    default:
      std::stringstream msg;
      msg << "Error: something wrong, check the geometry setup of the cubed "
             "sphere. \n";
      msg << "----------------------------------" << std::endl;
      ATHENA_ERROR(msg);
  }

  // Calculate ox1 and ox2
  if (target_block >
      0) {  // Need to change only when a special boundary is crossed
    // Calculate the lx1 and lx2 positions of the neighbor block
    // First calculate the top left corner position of the block
    int lx3_0, lx2_0;
    switch (target_block) {
      case 1:
        lx3_0 = 0;
        lx2_0 = 0;
        break;
      case 2:
        lx3_0 = 0;
        lx2_0 = 1;
        break;
      case 3:
        lx3_0 = 1;
        lx2_0 = 0;
        break;
      case 4:
        lx3_0 = 1;
        lx2_0 = 1;
        break;
      case 5:
        lx3_0 = 2;
        lx2_0 = 0;
        break;
      case 6:
        lx3_0 = 2;
        lx2_0 = 1;
        break;
      default:
        std::stringstream msg;
        msg << "Error: something wrong, check the geometry setup of the cubed "
               "sphere. \n";
        msg << "----------------------------------" << std::endl;
        ATHENA_ERROR(msg);
    }
    lx3_0 = (lx3_0 << (loc.level - 2));
    lx2_0 = (lx2_0 << (loc.level - 2));
    // Add up first block and local positions
    int lx3_t = lx3_0 + target_loc_3;
    int lx2_t = lx2_0 + target_loc_2;
    // Calculate and pass the differences
    *ox3 = lx3_t - loc.lx3;
    *ox2 = lx2_t - loc.lx2;
  } else {
    *tox2 = -*ox2;
    *tox3 = -*ox3;
  }
  // std::cout << "|Block ID: " << block_id << "|Target Block: " << target_block
  // << "||ox2: " << *ox2 << "|ox3: " << *ox3 << "||tox2:" << *tox2 << "|tox3:"
  // << *tox3 << std::endl;
  return;
}
