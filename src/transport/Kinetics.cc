#include "Kinetics.h"

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_gmres.h>

Kinetics::Kinetics(int nrows, int ncols, int nhalo)
    : m_grid(nrows, ncols, nhalo),
      m_psi(m_grid, "psi"),
      m_eddy(m_grid, "eddy"),
      m_q(m_grid, "q"),
      m_advection(m_grid, m_psi),
      m_diffusion(m_grid, m_eddy),
      m_pattern(m_grid.rank(), m_grid.rank(), 9),
      m_src(m_grid.rank()),
      m_rhs(m_grid.rank()),
      m_bnd(m_grid.rank()),
      m_buffer(m_grid.pattern()) {
  m_pattern.compress();
  m_adj.reinit(m_pattern);
  m_buffer = 0.;

  initialize();
  m_advection.update();
  m_diffusion.update();

  m_psi.finish();
  m_eddy.finish();
  m_q.finish();
}

void Kinetics::initialize() {
  m_psi.setZero();
  m_q.setZero();
  m_eddy.setOnes();
  m_eddy *= 0.1;
  for (int i = 0; i < m_grid.rows(); i++)
    for (int j = 0; j < m_grid.cols(); j++) {
      m_psi.val(i, j) = -(j + 1) * m_grid.dd();
      // m_psi.val(i, j) = 0.;
      if (i > m_grid.rows() / 4. && i < 3. * m_grid.rows() / 4. &&
          j > m_grid.cols() / 4. && j < 3. * m_grid.cols() / 4.)
        m_q.val(i, j) = 1.;
    }
  m_q.set_boundary(3, Dirichlet, 1.);
  for (int i = 0; i < m_grid.rank(); i++)
    m_src(i) = m_q.val(i / m_grid.cols(), i % m_grid.cols());
}

void Kinetics::assemble(double theta, double dt) {
  m_buffer.add(1., m_diffusion());
  m_buffer.add(-1., m_advection());
  m_buffer.mmult(m_adj, m_q.neumann());
  m_mass.reinit(m_pattern);
  m_force.reinit(m_pattern);

  m_mass = dealii::IdentityMatrix(m_grid.rank());
  m_mass.add(-theta * dt, m_adj);

  m_force = dealii::IdentityMatrix(m_grid.rank());
  m_force.add((1. - theta) * dt, m_adj);

  m_buffer.vmult(m_bnd, m_q.dirichlet());
  m_bnd *= dt;
}

void Kinetics::checkout() {
  std::cout << "== begin ==" << std::endl;
  std::cout << m_psi << std::endl;
  std::cout << m_eddy << std::endl;
  std::cout << m_q << std::endl;
  std::cout << "== source ==" << std::endl;
  std::cout << m_src << std::endl;
  /*
  std::cout << "== advection ==" << std::endl;
  m_advection().print_formatted(std::cout, 2, false);
  std::cout << "== diffusion ==" << std::endl;
  m_diffusion().print_formatted(std::cout, 2, false);
  std::cout << "== mass matrix ==" << std::endl;
  m_mass.print_formatted(std::cout, 2, false);
  std::cout << "== force matrix ==" << std::endl;
  m_force.print_formatted(std::cout, 2, false);
  */
  std::cout << "== end ==" << std::endl;
}

void Kinetics::run(int nt) {
  dealii::SolverControl control(1000, 1e-12);
  dealii::SolverGMRES<> solver(control);

  for (int t = 0; t < nt; t++) {
    m_force.vmult(m_rhs, m_src);
    m_rhs += m_bnd;
    solver.solve(m_mass, m_src, m_rhs, dealii::PreconditionIdentity());

    if (t % 10 == 0) {
      std::cout << m_src << std::endl;
      // std::cout << "Iteration converges after " << control.last_step()
      //     << " steps." << std::endl;
    }
  }
}
