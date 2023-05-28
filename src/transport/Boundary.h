#ifndef BOUNDARY_H_
#define BOUNDARY_H_

enum BoundaryType { Dirichlet, Neumann, Periodic };

template <typename Scalar>
struct BoundaryInfo {
  BoundaryType m_type;
  Scalar m_value;

  BoundaryInfo(BoundaryType type = Dirichlet, Scalar value = 0)
      : m_type(type), m_value(value) {}
};

#endif
