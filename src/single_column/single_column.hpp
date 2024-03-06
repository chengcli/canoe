// C/C++
#include <array>
#include <vector>

// athena
#include <athena/athena.hpp>

// canoe
#include <air_parcel.hpp>
#include <virtual_groups.hpp>

class MeshBlock;

class SingleColumn : public ParameterGroup {
 public:
  friend void find_tp_bottom(int n, double *x, double *f, void *arg);

  SingleColumn(MeshBlock *pmb, ParameterInput *pin);
  virtual ~SingleColumn();

  //! \param ac vector of air parcels [Mass Fraction]
  //! \param il lower index of the range
  //! \param iu upper index of the range
  void FindUnstableRange(AirColumn const &ac, int il, int iu,
                         std::vector<std::array<int, 2>> &ranges);

  std::array<Real, 2> FindMassEnthalpy(AirColumn const &ac, int k, int j,
                                       int il, int iu);

  //! \brief perform convective adjustment over a column range
  //! \param ac vector of air parcels [Mass Fraction]
  //! \param il lower index of the range
  //! \param iu upper index of the range
  void ConvectiveAdjustment(AirColumn &ac, int k, int j, int il, int iu);

 protected:  // convective adjustment functions
  //! Find the mass and enthalpy of an adiabatic air column
  //! The air parcel at the starting level is given by ac
  std::array<Real, 2> findAdiabaticMassEnthalpy(AirParcel air, int il, int iu);

 protected:
  //! ptr to meshblock
  MeshBlock *pmy_block_;

  //! scrach arrays for volumn
  AthenaArray<Real> vol_;
};
