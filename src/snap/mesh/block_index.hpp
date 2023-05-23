#ifndef BLOCK_INDEX_HPP
#define BLOCK_INDEX_HPP

#include <athena/defs.hpp>
#include <configure.hpp>

class MeshBlock;

class BlockIndex {
 public:
  enum Size { nghost = NGHOST };

  // index range
  int is, ie, js, je, ks, ke;

  // block ids
  int gid, lid;

  // array size
  int ncells1, ncells2, ncells3;

  // boundary condition flag
  bool apply_bc[6];

  // constructor
  BlockIndex(int nk, int nj, int ni, int mygid = 0, int mylid = 0);

  BlockIndex(MeshBlock *pmb);
};

#endif
