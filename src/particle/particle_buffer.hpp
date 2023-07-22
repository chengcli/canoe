/** @file particle_buffer.hpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Sunday May 30, 2021 15:05:42 PDT
 * @bug No known bugs.
 */

#ifndef PARTICLE_BUFFER_HPP
#define PARTICLE_BUFFER_HPP

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// C/C++ header
#include <vector>

// Athena++ header

class MaterialPoint;
class Particles;

//! Defines the class for managing buffers for transporting particles.
class ParticleBuffer {
 public:
  // data
  Particles *pmy_particle;

  // functions
  ParticleBuffer(Particles *ppart);
  ~ParticleBuffer();
  int CreateMPITag(int lid, int tid);
  void DetachParticle(std::vector<MaterialPoint> &mp);
  void SendParticle();
  void RecvParticle();
  bool AttachParticle(std::vector<MaterialPoint> &mp);
  void ClearBoundary();

 protected:
  enum BoundaryStatus particle_flag_[56];

  // resizeable send/recv buffer
  std::vector<MaterialPoint> particle_send_[56];
  std::vector<MaterialPoint> particle_recv_[56];

#ifdef MPI_PARALLEL
  MPI_Request req_particle_send_[56];
  MPI_Request req_particle_recv_[56];
#endif
};

#endif /* end of include guard PARTICLE_BUFFER_HPP */
