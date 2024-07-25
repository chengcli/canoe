// C/C++
#include <functional>
#include <iostream>
#include <sstream>
#include <string>

// athena
#include <athena/bvals/bvals.hpp>
#include <athena/globals.hpp>
#include <athena/mesh/mesh.hpp>

// application
#include <application/application.hpp>
#include <application/globals.hpp>

// canoe
#include <constants.hpp>
#include <impl.hpp>

// snap
#include <snap/stride_iterator.hpp>
#include <snap/thermodynamics/thermodynamics.hpp>

#include "decomposition.hpp"

Decomposition::Decomposition(MeshBlock *pmb)
    : pmy_block_(pmb), has_top_neighbor(false), has_bot_neighbor(false) {
  Application::Logger app("snap");
  app->Log("Initialize Decomposition");
  int nc1 = pmb->ncells1, nc2 = pmb->ncells2, nc3 = pmb->ncells3;
  // allocate hydrostatic and nonhydrostatic pressure
  psf_.NewAthenaArray(nc3, nc2, nc1 + 1);
  psv_.NewAthenaArray(nc3, nc2, nc1);
  dsf_.NewAthenaArray(nc3, nc2, nc1 + 1);
  dsv_.NewAthenaArray(nc3, nc2, nc1);
  pres_.NewAthenaArray(nc3, nc2, nc1);
  dens_.NewAthenaArray(nc3, nc2, nc1);

  buffer_ = new Real[3 * (NGHOST + 1) * nc3 * nc2];
  send_buffer_ = new Real[(NGHOST + 1) * nc3 * nc2];
  wsend_top_ = new Real[2 * NGHOST * nc3 * nc2];
  wrecv_top_ = new Real[2 * NGHOST * nc3 * nc2];
  wsend_bot_ = new Real[2 * NGHOST * nc3 * nc2];
  wrecv_bot_ = new Real[2 * NGHOST * nc3 * nc2];
  brank_ = new int[Globals::nranks];
  color_ = new int[Globals::nranks];

  // allocate polytropic index and pseudo entropy
  entropy_.NewAthenaArray(2, nc3, nc2);
  gamma_.NewAthenaArray(nc3, nc2, nc1);
}

Decomposition::~Decomposition() {
  Application::Logger app("snap");
  app->Log("Destroy Decomposition");

  delete[] buffer_;
  delete[] send_buffer_;
  delete[] wsend_top_;
  delete[] wrecv_top_;
  delete[] wsend_bot_;
  delete[] wrecv_bot_;
  delete[] brank_;
  delete[] color_;
}

int Decomposition::CreateMPITag(int recvid, int sendid, int phys) {
  std::string str = std::to_string(recvid);
  str += "x";
  str += std::to_string(sendid);
  str += "x";
  str += std::to_string(phys);
  return std::hash<std::string>{}(str) % Globals::mpi_tag_ub;
  // return (recvid<<11) | (sendid<<5) | phys;
}

void Decomposition::FindNeighbors() {
  MeshBlock *pmb = pmy_block_;
  has_top_neighbor = false;
  has_bot_neighbor = false;
  for (int n = 0; n < pmb->pbval->nneighbor; ++n) {
    NeighborBlock &nb = pmb->pbval->neighbor[n];
    if ((nb.ni.ox1 == -1) && (nb.ni.ox2 == 0) && (nb.ni.ox3 == 0)) {
      bblock = nb;
      has_bot_neighbor = true;
    }
    if ((nb.ni.ox1 == 1) && (nb.ni.ox2 == 0) && (nb.ni.ox3 == 0)) {
      tblock = nb;
      has_top_neighbor = true;
    }
  }

  if (!has_bot_neighbor) {
    bblock.snb.gid = -1;
    bblock.snb.rank = -1;
  }
  if (!has_top_neighbor) {
    tblock.snb.gid = -1;
    tblock.snb.rank = -1;
  }

#ifdef MPI_PARALLEL
  MPI_Allgather(&bblock.snb.rank, 1, MPI_INT, brank_, 1, MPI_INT,
                MPI_COMM_WORLD);
#else
  brank_[0] = -1;
#endif
}

// FIXME: local boundary has not been implemented
// Ordering the meshblocks need to be worked out such that
// the upper boundary executes before the lower boundary
void Decomposition::RecvFromTop(AthenaArray<Real> &psf, int kl, int ku, int jl,
                                int ju) {
  MeshBlock *pmb = pmy_block_;
  int ssize = (ju - jl + 1) * (ku - kl + 1) * (NGHOST + 1);

  std::stringstream msg;
#ifdef MPI_PARALLEL
  MPI_Status status;
#endif

  if (tblock.snb.rank != Globals::my_rank) {  // MPI boundary
#ifdef MPI_PARALLEL
    int tag = CreateMPITag(pmb->gid, tblock.snb.gid, 22);
    MPI_Recv(buffer_, ssize, MPI_ATHENA_REAL, tblock.snb.rank, tag,
             MPI_COMM_WORLD, &status);
#endif
  } else {  // local boundary
    // need to wait for the top boundary to finish
    msg << "### FATAL ERROR in Decompositin::RecvFromTop" << std::endl
        << "Local boundary not yet implemented" << std::endl;
    ATHENA_ERROR(msg);
  }
  int p = 0;
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      for (int i = 0; i <= NGHOST; ++i)
        psf(k, j, pmb->ie + i + 1) = buffer_[p++];
}

void Decomposition::SendToBottom(AthenaArray<Real> const &psf, int kl, int ku,
                                 int jl, int ju) {
  MeshBlock *pmb = pmy_block_;
  int ssize = 0;
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      for (int i = 0; i <= NGHOST; ++i)
        buffer_[ssize++] = psf(k, j, pmb->is + i);

  if (bblock.snb.rank != Globals::my_rank) {  // MPI boundary
#ifdef MPI_PARALLEL
    int tag = CreateMPITag(bblock.snb.gid, pmb->gid, 22);
    MPI_Isend(buffer_, ssize, MPI_ATHENA_REAL, bblock.snb.rank, tag,
              MPI_COMM_WORLD, &req_send_bot_);
#endif
  } else {  // local boundary
    MeshBlock *pbl = pmy_block_->pmy_mesh->FindMeshBlock(bblock.snb.gid);
    std::memcpy(pbl->pimpl->pdec->buffer_, buffer_, ssize * sizeof(Real));
  }
}

void Decomposition::SyncNewVariables(AthenaArray<Real> const &w, int kl, int ku,
                                     int jl, int ju) {
  MeshBlock *pmb = pmy_block_;

  if (has_bot_neighbor) {
    int sbot = 0;
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j)
        for (int i = 0; i < NGHOST; ++i) {
          wsend_bot_[sbot++] = w(IDN, k, j, pmb->is + i);
          wsend_bot_[sbot++] = w(IPR, k, j, pmb->is + i);
        }
    if (bblock.snb.rank != Globals::my_rank) {  // MPI boundary
#ifdef MPI_PARALLEL
      int tag = CreateMPITag(bblock.snb.gid, pmb->gid, 17);
      MPI_Isend(wsend_bot_, sbot, MPI_ATHENA_REAL, bblock.snb.rank, tag,
                MPI_COMM_WORLD, &req_send_sync_bot_);
      tag = CreateMPITag(pmb->gid, bblock.snb.gid, 19);
      MPI_Irecv(wrecv_bot_, sbot, MPI_ATHENA_REAL, bblock.snb.rank, tag,
                MPI_COMM_WORLD, &req_recv_sync_bot_);
#endif
    } else {  // local boundary
      MeshBlock *pbl = pmy_block_->pmy_mesh->FindMeshBlock(bblock.snb.gid);
      std::memcpy(pbl->pimpl->pdec->wrecv_top_, wsend_bot_,
                  sbot * sizeof(Real));
      std::memcpy(wrecv_bot_, pbl->pimpl->pdec->wsend_top_,
                  sbot * sizeof(Real));
    }
  }

  if (has_top_neighbor) {
    int stop = 0;
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j)
        for (int i = 0; i < NGHOST; ++i) {
          wsend_top_[stop++] = w(IDN, k, j, pmb->ie - i);
          wsend_top_[stop++] = w(IPR, k, j, pmb->ie - i);
        }
    if (tblock.snb.rank != Globals::my_rank) {  // MPI boundary
#ifdef MPI_PARALLEL
      int tag = CreateMPITag(tblock.snb.gid, pmb->gid, 19);
      MPI_Isend(wsend_top_, stop, MPI_ATHENA_REAL, tblock.snb.rank, tag,
                MPI_COMM_WORLD, &req_send_sync_top_);
      tag = CreateMPITag(pmb->gid, tblock.snb.gid, 17);
      MPI_Irecv(wrecv_top_, stop, MPI_ATHENA_REAL, tblock.snb.rank, tag,
                MPI_COMM_WORLD, &req_recv_sync_top_);
#endif
    } else {  // local boundary
      MeshBlock *pbl = pmy_block_->pmy_mesh->FindMeshBlock(bblock.snb.gid);
      std::memcpy(pbl->pimpl->pdec->wrecv_bot_, wsend_top_,
                  stop * sizeof(Real));
      std::memcpy(wrecv_top_, pbl->pimpl->pdec->wsend_bot_,
                  stop * sizeof(Real));
    }
  }
}

void Decomposition::WaitToFinishSend() {
#ifdef MPI_PARALLEL
  MPI_Status status;
  if (has_bot_neighbor && (bblock.snb.rank != Globals::my_rank))
    MPI_Wait(&req_send_bot_, &status);

  if (has_top_neighbor && (tblock.snb.rank != Globals::my_rank))
    MPI_Wait(&req_send_top_, &status);
#endif
}

void Decomposition::WaitToFinishSync(AthenaArray<Real> &w, int kl, int ku,
                                     int jl, int ju) {
  MeshBlock *pmb = pmy_block_;
#ifdef MPI_PARALLEL
  MPI_Status status;
  if (has_bot_neighbor && (bblock.snb.rank != Globals::my_rank)) {
    MPI_Wait(&req_send_sync_bot_, &status);
    MPI_Wait(&req_recv_sync_bot_, &status);
  }
#endif

  if (has_bot_neighbor) {
    int p = 0;
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j)
        for (int i = 1; i <= NGHOST; ++i) {
          w(IDN, k, j, pmb->is - i) = wrecv_bot_[p++];
          w(IPR, k, j, pmb->is - i) = wrecv_bot_[p++];
        }
  }

#ifdef MPI_PARALLEL
  if (has_top_neighbor && (tblock.snb.rank != Globals::my_rank)) {
    MPI_Wait(&req_send_sync_top_, &status);
    MPI_Wait(&req_recv_sync_top_, &status);
  }
#endif

  if (has_top_neighbor) {
    int p = 0;
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j)
        for (int i = 1; i <= NGHOST; ++i) {
          w(IDN, k, j, pmb->ie + i) = wrecv_top_[p++];
          w(IPR, k, j, pmb->ie + i) = wrecv_top_[p++];
        }
  }
}

// FIXME: local boundary has not been implemented
// Ordering the meshblocks need to be worked out such that
// the upper boundary executes before the lower boundary
void Decomposition::RecvBuffer(AthenaArray<Real> &psf, int kl, int ku, int jl,
                               int ju, int il, int iu, NeighborBlock nb) {
  MeshBlock *pmb = pmy_block_;
  int ssize = (iu - il + 1) * (ju - jl + 1) * (ku - kl + 1);

  std::stringstream msg;
#ifdef MPI_PARALLEL
  MPI_Status status;
#endif

  if (nb.snb.rank != Globals::my_rank) {  // MPI boundary
#ifdef MPI_PARALLEL
    int tag = CreateMPITag(pmb->gid, nb.snb.gid, 22);
    MPI_Recv(buffer_, ssize, MPI_ATHENA_REAL, nb.snb.rank, tag, MPI_COMM_WORLD,
             &status);
#endif
  } else {  // local boundary
    // need to wait for the top boundary to finish
    msg << "### FATAL ERROR in Decompositin::RecvBuffer" << std::endl
        << "Local boundary not yet implemented" << std::endl;
    ATHENA_ERROR(msg);
  }
  int p = 0;
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      for (int i = il; i <= iu; ++i) psf(k, j, i) = buffer_[p++];
}

void Decomposition::SendBuffer(AthenaArray<Real> const &psf, int kl, int ku,
                               int jl, int ju) {
  MeshBlock *pmb = pmy_block_;
  int s1 = 0;
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      for (int i = pmb->is; i <= pmb->is + NGHOST; ++i)
        send_buffer_[s1++] = psf(k, j, i);

  int s2 = s1;
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      for (int i = pmb->ie - NGHOST + 1; i <= pmb->ie + 1; ++i)
        send_buffer_[s2++] = psf(k, j, i);

  if (has_bot_neighbor) {
    if (bblock.snb.rank != Globals::my_rank) {  // MPI boundary
      int tag = CreateMPITag(bblock.snb.gid, pmb->gid, 22);
#ifdef MPI_PARALLEL
      MPI_Isend(send_buffer_, s1, MPI_ATHENA_REAL, bblock.snb.rank, tag,
                MPI_COMM_WORLD, &req_send_bot_);
#endif
    } else {  // local boundary, place holder, may not work
      MeshBlock *pbl = pmy_block_->pmy_mesh->FindMeshBlock(bblock.snb.gid);
      std::memcpy(pbl->pimpl->pdec->buffer_, send_buffer_, s1 * sizeof(Real));
    }
  }

  if (has_top_neighbor) {
    if (tblock.snb.rank != Globals::my_rank) {
      int tag = CreateMPITag(tblock.snb.gid, pmb->gid, 22);
#ifdef MPI_PARALLEL
      MPI_Isend(send_buffer_ + s1, s1, MPI_ATHENA_REAL, tblock.snb.rank, tag,
                MPI_COMM_WORLD, &req_send_top_);
#endif
    } else {  // local boundary, place holder, may not work
      MeshBlock *pbl = pmy_block_->pmy_mesh->FindMeshBlock(tblock.snb.gid);
      std::memcpy(pbl->pimpl->pdec->buffer_, send_buffer_ + s1,
                  s1 * sizeof(Real));
    }
  }
}

void Decomposition::PopulateBotEntropy(AthenaArray<Real> const &w, int kl,
                                       int ku, int jl, int ju) {
  MeshBlock *pmb = pmy_block_;
  int c = 0;
  for (int i = 0; i < Globals::nranks; ++i) {
    if (brank_[i] == -1)
      color_[i] = c++;
    else
      color_[i] = color_[brank_[i]];
  }

  int ssize = 2 * (ku - kl + 1) * (ju - jl + 1), p = 0;
  auto pthermo = Thermodynamics::GetInstance();
  if (!has_bot_neighbor) {
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j) {
        Real gamma = pthermo->GetGamma(w.at(k, j, pmb->is));
        buffer_[p++] = gamma;
        buffer_[p++] =
            log(w(IPR, k, j, pmb->is)) - gamma * log(w(IDN, k, j, pmb->is));
      }
  }

  /*if (Globals::my_rank == 0) {
    for (int i = 0; i < Globals::nranks; ++i)
      std::cout << brank_[i] << " ";
    std::cout << std::endl;
    for (int i = 0; i < Globals::nranks; ++i)
      std::cout << color_[i] << " ";
    std::cout << std::endl;
  }*/

#ifdef MPI_PARALLEL
  MPI_Comm comm_col;
  MPI_Comm_split(MPI_COMM_WORLD, color_[Globals::my_rank], Globals::my_rank,
                 &comm_col);
  // assuming correct ordering
  MPI_Bcast(buffer_, ssize, MPI_ATHENA_REAL, 0, comm_col);
  MPI_Comm_free(&comm_col);
#endif

  p = 0;
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j) {
      entropy_(0, k, j) = buffer_[p++];
      entropy_(1, k, j) = buffer_[p++];
    }
}
