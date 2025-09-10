#pragma once

#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "adcpp.hpp"

// class definition

namespace IceShell {
  template<class Real>
  struct PoissonValueBoundary {
    int cell_i;
    int boundary_i;
    Real center_dist;
    Real area;
  };

  template<class Real>
  struct PoissonZeroValueBoundary {
    int i;
    Real center_dist;
    Real area;
  };

  template<class Real>
  struct PoissonLink {
    int i;
    int j;
    Real center_dist;
    Real area;
  };

  template<class Real>
  class PoissonEquation{
    public:

      const int nx;
      const int nz;
      const int nc;
      const int nb;

      PoissonEquation(
          const std::vector<Real> & dx, const std::vector<Real> & dz);
      auto get_boundary_flux(int j) const;
      auto get_greens_func();
      void solve(const std::vector<Real> &boundary_x);

    private:

      typedef Eigen::Vector<Real, Eigen::Dynamic> EigenVector;
      typedef Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> EigenMatrix;
      typedef Eigen::Triplet<Real> EigenTriplet;
      typedef Eigen::SparseMatrix<Real> EigenSparseMatrix;
      typedef Eigen::SimplicialLDLT<EigenSparseMatrix> EigenSparseMatrixSolver;

      typedef PoissonValueBoundary<Real> ValueBoundary;
      typedef PoissonZeroValueBoundary<Real> ZeroValueBoundary;
      typedef PoissonLink<Real> Link;

      std::vector<Real> dx;
      std::vector<Real> dz;

      std::vector<ValueBoundary> value_boundaries;
      std::vector<ZeroValueBoundary> zero_value_boundaries;
      std::vector<Link> links;

      EigenVector boundary_x;
      EigenVector interior_x;
      EigenVector interior_b;
      EigenSparseMatrix laplacian;
      EigenSparseMatrixSolver inverse_laplacian;

      inline int cell_index(int i, int k) {
        return i * nz + k;
      }

      void clear_buffer();
      void init_grid();
      void init_laplacian();
      void init_inverse_laplacian();
  };

  template <class Real>
  class TemperatureDependentDiffusivity {
    public:
      const Real temp_ref; // Kelvin
      const Real k0; // Watts per meter

      TemperatureDependentDiffusivity(Real temp_ref, Real k0):
        temp_ref(temp_ref), k0(k0) {};

      template<class Number>
      inline auto tau(Number temp) const {
        return k0 * log(temp / temp_ref);
      }

      template<class Number>
      inline auto temp(Number tau) const {
        return temp_ref * exp(tau / k0);
      }
  };

  template <class Real>
  class IceDiffusion {
    public:

      typedef Eigen::Vector<Real, Eigen::Dynamic> EigenVector;
      typedef Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> EigenMatrix;

      IceDiffusion(
        const TemperatureDependentDiffusivity<Real> tdd,
        const std::vector<Real> & dx, const std::vector<Real> & dz
      );

      void compute_flux(EigenVector & f, EigenMatrix & df_dt,
          const EigenVector & t);

    private:

      const EigenMatrix greens_func;
      const TemperatureDependentDiffusivity<Real> tdd;
      const int nb;
      EigenVector tau;
      EigenVector dtau_dt;

      static auto get_Poisson_equation_greens_function(
          const std::vector<Real> & dx, const std::vector<Real> & dz);
      static auto get_tdd(Real temp_ref, Real ice_d);
  };

  template<class Real>
  class VaporCondensation {
    public:
      Real p3;
      Real temp3;
      Real gas_constant;
      Real gamma;
      Real beta;
      Real delta;
      VaporCondensation(Real p3, Real temp3,
        Real gas_constant, Real gamma, Real beta, Real delta):
          p3(p3), temp3(temp3), gas_constant(gas_constant),
          gamma(gamma), beta(beta), delta(delta) {}

      template<class R>
      inline auto p_sat(const R &temp) const {
        auto t3 = temp / temp3;
        return p3 * exp(beta * (1. - 1./t3) - delta * log(t3));
      }

      template<class R1, class R2>
      inline auto specific_enthalpy_diff(
          const R1 &ice_temp, const R2 &air_temp) const {
        return gas_constant * (
            gamma / (gamma - 1.) * (air_temp - ice_temp)
            + beta * temp3 - delta * ice_temp
        );
      }

      template<class R>
      inline auto one_side_vapor_flux(const R &temp) const {
        return p_sat(temp) / sqrt(2 * M_PI * gas_constant * temp);
      }

      template<class R>
      inline auto one_side_vapor_flux(const R &temp, const R &pres) const {
        return pres / sqrt(2 * M_PI * gas_constant * temp);
      }

      template<class R1, class R2, class R3>
      inline auto net_vapor_flux(
          const R1 &ice_temp, const R2 &air_temp, const R3 &vapor_p) const {
        return (
          one_side_vapor_flux(air_temp, vapor_p)
          - one_side_vapor_flux(ice_temp)
        );
      }

      template<class R1, class R2, class R3>
      inline auto energy_flux(
          const R1 &ice_temp, const R2 &air_temp, const R3 &vapor_p) const {
        return (
            net_vapor_flux(ice_temp, air_temp, vapor_p)
            * specific_enthalpy_diff(ice_temp, air_temp)
        );
      }
  };

  template<class Real>
  class OuterSurfaceRadiation {
    public:
      Real effective_temp;
      Real stefan_boltzmann_const;

      OuterSurfaceRadiation(Real effective_temp, Real stefan_boltzmann_const):
        effective_temp(effective_temp),
        stefan_boltzmann_const(stefan_boltzmann_const) {}

      template<class R>
      inline auto pow4(const R &x) const {
        auto x2 = x * x;
        return x2 * x2;
      }

      template<class R>
      inline auto energy_flux(const R &ice_temp) const {
        return stefan_boltzmann_const * (pow4(effective_temp) - pow4(ice_temp));
      }
  };

  template<class Real>
  class BoundaryValue {
    public:

      typedef Eigen::Vector<Real, Eigen::Dynamic> EigenVector;

      const int nx;
      const int nz;
      const int nb;
      EigenVector data;

      BoundaryValue(int nx, int nz);

      inline int index_side(int k) const {
        return k;
      }

      inline int index_top(int i) const {
        return i + nz;
      }

      inline bool is_top(int i) const {
        return i >= nz;
      }

      inline void set_side(int k, Real t) {
        data(index_side(k)) = t;
      }

      inline void set_top(int i, Real t) {
        data(index_top(i)) = t;
      }

      inline Real get_side(int k) {
        return data(index_side(k));
      }

      inline Real get_top(int i) {
        return data(index_top(i));
      }
  };

  template<class Real>
  class IceAirBoundary {
    public:
      typedef Eigen::Vector<Real, Eigen::Dynamic> EigenVector;
      typedef Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> EigenMatrix;
      typedef BoundaryValue<Real> BV;

      IceAirBoundary(int nx, int nz,
          VaporCondensation<Real> cond, OuterSurfaceRadiation<Real> rad);

      void compute_flux(EigenVector & f, EigenVector & df_dt,
            const EigenVector & ice_temp,
            const BV & air_temp,
            const BV & vapor_p) const;

      VaporCondensation<Real> cond;
      OuterSurfaceRadiation<Real> rad;
    private:
  };

  template<class Real>
  class IceBoundaryModel {
    public:

      typedef IceAirBoundary<Real> AirB;
      typedef IceDiffusion<Real> IceDiff;

      typedef BoundaryValue<Real> BV;
      typedef VaporCondensation<Real> Cond;
      typedef OuterSurfaceRadiation<Real> Rad;
      typedef TemperatureDependentDiffusivity<Real> TDD;

      typedef Eigen::Vector<Real, Eigen::Dynamic> EigenVector;
      typedef Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> EigenMatrix;

      const int nx;
      const int nz;
      const int nb;

      AirB ice_air_boundary;
      IceDiff ice_diffusion;

      IceBoundaryModel(
          const std::vector<Real> & dx, const std::vector<Real> & dz,
          Cond cond, Rad rad, TDD tdd);

      void solve(BV & ice_t, const BV & air_t, const BV &vapor_p,
          const Real abs_tol=1e-5,  const int max_iter = 20);

    private:
      EigenVector t;
      EigenVector fa;
      EigenVector fi;
      EigenVector r;
      EigenVector dfa;
      EigenMatrix dfi;
      EigenMatrix dr;

      void init_guess(const BV & air_t);
      void nr_iterate(const BV & air_t, const BV & vapor_p);
  };
}

// class definition

namespace IceShell {
  template<class Real>
  PoissonEquation<Real>::PoissonEquation(
      const std::vector<Real> & dx, const std::vector<Real> & dz):
    dx(dx), dz(dz), nx(dx.size()), nz(dz.size()), nc(nx * nz), nb (nx + nz),
    boundary_x(nb), interior_x(nc), interior_b(nc), laplacian(nc, nc) {

    clear_buffer();
    init_grid();
    init_laplacian();
    init_inverse_laplacian();
  }

  template<class Real>
  auto PoissonEquation<Real>::get_boundary_flux(int j) const {
    auto i = value_boundaries[j];
    return (
        interior_x(i.cell_i) - boundary_x(i.boundary_i)
      ) / i.center_dist;
  }

  template<class Real>
  auto PoissonEquation<Real>::get_greens_func() {
    EigenMatrix greens_func(nb, nb);
    std::vector<Real> x(nb, 0.);
    for (int i = 0; i < nb; ++i) {
      std::fill(x.begin(), x.end(), 0.);
      x[i] = 1.;
      solve(x);
      for (int j = 0; j < nb; ++j) {
        greens_func(j, i) = get_boundary_flux(j);
      }
    }
    clear_buffer();
    return greens_func;
  }

  template<class Real>
  void PoissonEquation<Real>::solve(const std::vector<Real> &boundary_x) {
    for (int i = 0; i < nb; ++i) {
      this->boundary_x(i) = boundary_x[i];
    }

    interior_b.setZero();
    for (auto & i : value_boundaries) {
      interior_b(i.cell_i) += (
          boundary_x[i.boundary_i] * i.area / i.center_dist
      );
    }

    interior_x = inverse_laplacian.solve(interior_b);
  }

  template<class Real>
  void PoissonEquation<Real>::clear_buffer() {
    boundary_x.setZero();
    interior_x.setZero();
    interior_b.setZero();
  }

  template<class Real>
  void PoissonEquation<Real>::init_grid() {
    for (int i = 0; i < nx; ++i) {
      ZeroValueBoundary bc = {
        cell_index(i, 0), dz[0], dx[i]
      };
      zero_value_boundaries.push_back(bc);
    }

    int l = 0;

    for (int k = 0; k < nz; ++k) {
      ValueBoundary bc = {
        cell_index(0, k), l++, dx[0], dz[k]
      };
      value_boundaries.push_back(bc);
    }

    for (int i = 0; i < nx; ++i) {
      ValueBoundary bc = {
        cell_index(i, nz - 1), l++, dz[nz - 1], dx[i]
      };
      value_boundaries.push_back(bc);
    }

    assert(l == nb);

    for (int i = 0; i < nx - 1; ++i) {
      for (int k = 0; k < nz; ++k) {
        int i_n = i + 1;
        Link link = {
            cell_index(i, k), cell_index(i_n, k),
            (dx[i] + dx[i_n]) / 2, dz[k]
        };
        links.push_back(link);
      }
    }

    for (int i = 0; i < nx; ++i) {
      for (int k = 0; k < nz - 1; ++k) {
        int k_n = k + 1;
        Link link = {
            cell_index(i, k), cell_index(i, k_n),
            (dz[k] + dz[k_n]) / 2, dx[i]
        };
        links.push_back(link);
      }
    }
  }

  template<class Real>
  void PoissonEquation<Real>::init_laplacian() {
    std::vector<EigenTriplet> coefs;

    for (auto & i : value_boundaries) {
      coefs.push_back(EigenTriplet(i.cell_i, i.cell_i, i.area / i.center_dist));
    }

    for (auto & i : zero_value_boundaries) {
      coefs.push_back(EigenTriplet(i.i, i.i, i.area / i.center_dist));
    }

    for (auto & link : links) {
      Real f = link.area / link.center_dist;
      coefs.push_back(EigenTriplet(link.i, link.j, -f));
      coefs.push_back(EigenTriplet(link.j, link.i, -f));
      coefs.push_back(EigenTriplet(link.i, link.i, f));
      coefs.push_back(EigenTriplet(link.j, link.j, f));
    }

    laplacian.setFromTriplets(coefs.begin(), coefs.end());

    laplacian.makeCompressed();
  }

  template<class Real>
  void PoissonEquation<Real>::init_inverse_laplacian() {
    inverse_laplacian.compute(laplacian);
  }

  template <class Real>
  IceDiffusion<Real>::IceDiffusion(
      const TemperatureDependentDiffusivity<Real> tdd,
      const std::vector<Real> & dx, const std::vector<Real> & dz):
    tdd(tdd),
    greens_func(get_Poisson_equation_greens_function(dx, dz)),
    nb(dx.size() + dz.size()) {
      tau = EigenVector::Zero(nb);
      dtau_dt = EigenVector::Zero(nb);
  }

  template <class Real>
  auto IceDiffusion<Real>::get_tdd(Real temp_ref, Real ice_d) {
    return TemperatureDependentDiffusivity<Real>(temp_ref, ice_d);
  }

  template <class Real>
  auto IceDiffusion<Real>::get_Poisson_equation_greens_function(
          const std::vector<Real> & dx, const std::vector<Real> & dz) {
    PoissonEquation<Real> p (dx, dz);
    return p.get_greens_func();
  }

  template <class Real>
  void IceDiffusion<Real>::compute_flux(EigenVector & f, EigenMatrix & df_dt,
      const EigenVector & t) {
    typedef adcpp::fwd::Number<Real> Dual;
    for (int i = 0; i < nb; ++i) {
        Dual t_ad(t(i), 1.);
        auto tau_ad = tdd.tau(t_ad);
        tau(i) = tau_ad.value();
        dtau_dt(i) = tau_ad.derivative();
    }
    f = greens_func * tau;
    df_dt = greens_func * dtau_dt.asDiagonal();
  }

  template<class Real>
  BoundaryValue<Real>::BoundaryValue(int nx, int nz):
    nx(nx), nz(nz), nb(nx + nz),
    data(nb) {
      data.setZero();
  }

  template<class Real>
  IceAirBoundary<Real>::IceAirBoundary(int nx, int nz,
      VaporCondensation<Real> cond, OuterSurfaceRadiation<Real> rad):
    cond(cond), rad(rad) {
  }

  template<class Real>
  void IceAirBoundary<Real>::compute_flux(EigenVector & f, EigenVector & df_dt,
        const EigenVector & ice_temp,
        const BV & air_temp, const BV & vapor_p) const {
    typedef adcpp::fwd::Number<Real> Dual;
    for (int i = 0; i < air_temp.nb; ++i) {

      Dual temp_ad(ice_temp(i), 1.);

      Dual f_ad = cond.energy_flux(temp_ad, air_temp.data(i), vapor_p.data(i));
      if (air_temp.is_top(i)) {
        f_ad += rad.energy_flux(temp_ad);
      }

      f(i) = f_ad.value();
      df_dt(i) = f_ad.derivative();
    }
  }

  template<class Real>
  IceBoundaryModel<Real>::IceBoundaryModel(
      const std::vector<Real> & dx, const std::vector<Real> & dz,
      Cond cond, Rad rad, TDD tdd):
      nx(dx.size()), nz(dz.size()), nb(nx + nz),
      ice_air_boundary(AirB(nx, nz, cond, rad)),
      ice_diffusion(IceDiff(tdd, dx, dz)),
      t(nb), fa(nb), fi(nb), r(nb), dfa(nb), dfi(nb, nb), dr(nb, nb) {
    t.setZero();
    fa.setZero();
    fi.setZero();
    r.setZero();
    dfa.setZero();
    dfi.setZero();
    dr.setZero();
  }

  template<class Real>
  void IceBoundaryModel<Real>::solve(BV & ice_t,
          const BV & air_t, const BV & vapor_p,
          const Real abs_tol, const int max_iter) {

    init_guess(air_t);

    bool solved = false;

    for (int i = 0; i < max_iter; ++i) {
      nr_iterate(air_t, vapor_p);
      Real abs_error = r.template lpNorm<Eigen::Infinity>();
      if (abs_error < abs_tol) {
        solved = true;
        break;
      }
    }

    if (!solved) {
      Real abs_error = r.template lpNorm<Eigen::Infinity>();
      std::cout << "IceBoundaryModel: abs_error = " << abs_error << std::endl;
    }

    ice_t.data = t;
  }

  template<class Real>
  void IceBoundaryModel<Real>::init_guess(const BV & air_t) {
    t = air_t.data;
  }

  template<class Real>
  void IceBoundaryModel<Real>::nr_iterate(const BV &air_t, const BV &vapor_p) {
    const Real t_max = 300.;
    const Real t_min = 60.;

    ice_air_boundary.compute_flux(fa, dfa, t, air_t, vapor_p);
    ice_diffusion.compute_flux(fi, dfi, t);
    r = fa + fi;
    dr = dfi;
    dr += dfa.asDiagonal();
    t -= dr.completeOrthogonalDecomposition().solve(r);
    t = t.cwiseMin(t_max).cwiseMax(t_min);
  }
}
