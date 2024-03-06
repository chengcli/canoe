#ifndef SRC_AIR_PARCEL_HPP_
#define SRC_AIR_PARCEL_HPP_

// C/C++
#include <array>
#include <iostream>
#include <vector>

// athena
#include <athena/athena.hpp>  // Real

// canoe
#include <configure.hpp>

//! \class AirParcel
//  \brief a collection of all physical data in a computational cell
class AirParcel {
 public:
  enum { Size = NHYDRO + NCLOUD + NCHEMISTRY + NTRACER + NTURBULENCE };

  enum class Type { MassFrac = 0, MassConc = 1, MoleFrac = 2, MoleConc = 3 };

  friend std::ostream &operator<<(std::ostream &os, Type const &type);
  friend std::ostream &operator<<(std::ostream &os, AirParcel const &var);

 protected:
  // data holder
  std::array<Real, Size> data_;

  // type
  Type mytype_;

 public:
  //! data pointers
  //! hydro data
  Real *const w;

  //! cloud data
  Real *const c;

  //! chemistry data
  Real *const q;

  //! tracer data
  Real *const x;

  //! turbulence data
  Real *const t;

  //! particle data
  Real const *d;

  // constructor
  explicit AirParcel(Type type = Type::MoleFrac)
      : mytype_(type),
        w(data_.data()),
        c(w + NHYDRO),
        q(c + NCLOUD),
        x(q + NCHEMISTRY),
        t(x + NTRACER),
        d(t + NTURBULENCE) {
    std::fill(data_.begin(), data_.end(), 0.0);
  }

  // copy constructor
  AirParcel(AirParcel const &other)
      : mytype_(other.mytype_),
        data_(other.data_),
        w(data_.data()),
        c(w + NHYDRO),
        q(c + NCLOUD),
        x(q + NCHEMISTRY),
        t(x + NTRACER),
        d(t + NTURBULENCE) {}

  // Assignment operator
  AirParcel &operator=(const AirParcel &other) {
    // Check for self-assignment
    if (this == &other) {
      return *this;
    }

    data_ = other.data_;
    mytype_ = other.mytype_;
    return *this;
  }

  void SetType(Type type) { mytype_ = type; }

  Type GetType() const { return mytype_; }

  void SetZero() { std::fill(data_.begin(), data_.end(), 0.0); }

  AirParcel &ConvertTo(AirParcel::Type type);

  AirParcel &ToMassFraction();
  AirParcel &ToMassConcentration();
  AirParcel &ToMoleFraction();
  AirParcel &ToMoleConcentration();

  // thermodynamic functions
  Real gammad() const;
  Real chi() const;
  Real theta(Real p0) const;

 protected:
  void massFractionToMassConcentration();
  void massConcentrationToMassFraction();

  void massFractionToMoleFraction();
  void moleFractionToMassFraction();

  void massConcentrationToMoleFraction();
  void moleFractionToMassConcentration();

  void moleFractionToMoleConcentration();
  void moleConcentrationToMoleFraction();

  void massConcentrationToMoleConcentration();
  void moleConcentrationToMassConcentration();

  void massFractionToMoleConcentration();
  void moleConcentrationToMassFraction();
};

using AirColumn = std::vector<AirParcel>;

namespace AirParcelHelper {
AirParcel gather_from_primitive(MeshBlock const *pmb, int k, int j, int i);
AirColumn gather_from_primitive(MeshBlock const *pmb, int k, int j);

inline AirColumn gather_from_primitive(MeshBlock const *pmb, int k, int j,
                                       int il, int iu) {
  AirColumn ac(iu - il + 1);
  for (int i = il; i <= iu; ++i) {
    ac[i - il] = gather_from_primitive(pmb, k, j, i);
  }

  return ac;
}

AirParcel gather_from_conserved(MeshBlock const *pmb, int k, int j, int i);
AirColumn gather_from_conserved(MeshBlock const *pmb, int k, int j);

inline AirColumn gather_from_conserved(MeshBlock const *pmb, int k, int j,
                                       int il, int iu) {
  AirColumn ac(iu - il + 1);
  for (int i = il; i <= iu; ++i) {
    ac[i - il] = gather_from_conserved(pmb, k, j, i);
  }

  return ac;
}

void distribute_to_primitive(MeshBlock *pmb, int k, int j, int i,
                             AirParcel const &air_in);
void distribute_to_primitive(MeshBlock *pmb, int k, int j, AirColumn const &ac);

inline void distribute_to_primitive(MeshBlock *pmb, int k, int j, int il,
                                    int iu, AirColumn const &ac) {
  for (int i = il; i <= iu; ++i) {
    distribute_to_primitive(pmb, k, j, i, ac[i - il]);
  }
}

void distribute_to_conserved(MeshBlock *pmb, int k, int j, int i,
                             AirParcel const &air_in);
void distribute_to_conserved(MeshBlock *pmb, int k, int j, AirColumn const &ac);

inline void distribute_to_conserved(MeshBlock *pmb, int k, int j, int il,
                                    int iu, AirColumn const &ac) {
  for (int i = il; i <= iu; ++i) {
    distribute_to_conserved(pmb, k, j, i, ac[i - il]);
  }
}

}  // namespace AirParcelHelper

#endif  // SRC_AIR_PARCEL_HPP_
