/* -------------------------------------------------------------------------------------
 * SNAP Example Program
 *
 * Contributer: Cheng Li, University of Michigan, 2021
 * Contact: chengcli@umich.edu
 * Reference: Straka et al., 1993
 * -------------------------------------------------------------------------------------
 */

// @sect3{Include files}

// First we include some Athena++ headers files. They contain declarations
// of various components of the solver system.
// Athena++ is able to run both single precision (float) and double precision (double)
// applications.
// A macro <code>Real</code> is defined to indicate the precision and it is define in
// the following header file.
#include <athena/athena.hpp>
// Then the basic multi-dimension array that stores fluid dynamic data is declared in
// the AthenaArray class.
#include <athena/athena_arrays.hpp>
#include <athena/stride_iterator.hpp>
// The model executes according to the parameters specified in an input file.
// The input file and the parameters within are managed by the ParameterInput class.
#include <athena/parameter_input.hpp>
// Coordinate related information are stored in the Coordinates class.
#include <athena/coordinates/coordinates.hpp>
// Everying regarding the equation of state is treated by the EquationOfState class and
// declared in the following file.
#include <athena/eos/eos.hpp>
// Evolution of hydrodynamic fields like pressure, density and velocities are managed by
// the Hydro class.
#include <athena/hydro/hydro.hpp>
// Athena++ can alo simulate magnetohydrodynamics (MHD).
// However, since this example is a hydro-only application, the MHD file is only
// included as a place holder.
#include <athena/field/field.hpp>
// Dynamic fields are evolving on a mesh and this is the file working with meshes and
// the partition of meshes among processors
#include <athena/mesh/mesh.hpp>

// The following header file contains the definition of the MeshBlock::Impl class
#include <impl.hpp>

// Finally, the Thermodynamics class works with thermodynamic aspects of the problem
// such as the temperature, potential temperature, condensation of vapor, etc.
#include <snap/thermodynamics/thermodynamics.hpp>
#include <snap/thermodynamics/thermodynamics_helper.hpp>

// Functions in the math library are protected by a specific namespace because math
// functions are usually concise in the names, such as <code>min</code>,
// <code>sqr</code>. It is better to protect them in a namespace to avoid conflicts.
#include <climath/core.h>

// Now, we get to the real program.
// First, we define several global variables specific to this appliciation
// For example, here <code>K</code> is the kinematic visocisty, <code>p0</code> is the
// surface pressure, <code>cp</code> is the specific heat capacity, and <code>Rd</code>
// is the ideal gas constant of dry air.
// The <a href="https://en.wikipedia.org/wiki/Prandtl_number">Prandtl number</a> is 1.
Real K, p0, cp, Rd;

// @sect3{User-defined output variables}

// The hydrodynamic solver evolves density, pressure and velocities with time, meaning that the
// (potential) temperature is a diagnostic quantity for output only. This block of code
// allocates the memory for the outputs of temperature and potential temperature.
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  AllocateUserOutputVariables(2);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
}

// Since temperature and potential temperature are only used for output purpose, they do
// not need to be calculated every dynamic time step. The subroutine below loops over
// all grids and calculates the value of temperature and potential temperature before the output time step.
// Particularly, the pointer to the Thermodynamics class
// <code>pthermo</code> calculates the temperature and potential temperature using its
// own member function Thermodynamics::GetTemp and PotentialTemp.
// <code>pthermo</code>, is a member of a
// higher level management class MeshBlock. So, you can use it directly inside a member
// function of class MeshBlock.
// In fact, all physics modules are managed by MeshBlock. As long as you have a pointer to a
// MeshBlock, you can access all physics in the simulation.
void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
{
  auto pthermo = pimpl->pthermo;

  // Loop over the grids. <code>js,je,is,ie</code> are members of the MeshBlock class.
  // They are integer values representing the start index and the end index of the grids
  // in each dimension.
  for (int j = js; j <= je; ++j)
    for (int i = is; i <= ie; ++i) {
      // <code>user_out_var</code> stores the actual data.
      // <code>phydro</code> is a pointer to the Hydro class, which has a member
      // <code>w</code> that stores density, pressure, and velocities at each grid.
      user_out_var(0,j,i) = pthermo->GetTemp(phydro->w.at(j,i));
      user_out_var(1,j,i) = pthermo->PotentialTemp(phydro->w.at(j,i), p0);
    }
}

// @sect3{Forcing function}

// The sinking bubble is forced by the gravity and the viscous and thermal dissipation.
// The gravitational forcing is taken care of by other parts of the program and we only
// need to write a few lines of code to facilitate dissipation.
// The first step is to define a function that takes a pointer to the MeshBlock as an
// argument such that we can access all physics via this pointer.
// The name of this function is not important. It can be anything.
// But the order and types of the arguments must be <code>(MeshBlock *, Real const, Real
// const, AthenaArray<Real> const&, AthenaArray<Real> const&, AthenaArray<Real> &)</code>.
// They are called the signature of the function.
void Diffusion(MeshBlock *pmb, Real const time, Real const dt,
  AthenaArray<Real> const& w, AthenaArray<Real> const& r,
  AthenaArray<Real> const& bcc, AthenaArray<Real> &u, AthenaArray<Real> &s)
{
  // <code>pcoord</code> is a pointer to the Coordinates class and it is a member of
  // the MeshBlock class.
  // We use the pointer to the MeshBlock class, <code>pmb</code> to access
  // <code>pcoord</code> and use its member function to get the spacing of the grid.
  Real dx = pmb->pcoord->dx1f(pmb->is);
  Real dy = pmb->pcoord->dx2f(pmb->js);

  // Similarly, we use <code>pmb</code> to find the pointer to the Thermodynamics class,
  // <code>pthermo</code>.
  auto pthermo = pmb->pimpl->pthermo;

  // Loop over the grids.
  for (int j = pmb->js; j <= pmb->je; ++j)
    for (int i = pmb->is; i <= pmb->ie; ++i) {

      // Similar to what we have done in MeshBlock::UserWorkBeforeOutput, we use the
      // Thermodynamics class to calculate temperature and potential temperature.
      Real temp = pthermo->GetTemp(w.at(j,i));
      Real theta = pthermo->PotentialTemp(w.at(j,i), p0);

      // The thermal diffusion is applied to the potential temperature field, which is
      // not exactly correct. But this is the setting of the test program.
      Real theta_ip1_j = pthermo->PotentialTemp(w.at(j+1,i), p0);
      Real theta_im1_j = pthermo->PotentialTemp(w.at(j-1,i), p0);
      Real theta_i_jp1 = pthermo->PotentialTemp(w.at(j,i+1), p0);
      Real theta_i_jm1 = pthermo->PotentialTemp(w.at(j,i-1), p0);

      // Add viscous dissipation to the velocities. Now you encounter another variable
      // called <code>u</code>. <code>u</code> stands for `conserved variables`, which
      // are density, momentums and total energy (kinetic + internal). In contrast,
      // <code>w</code> stands for `primitive variables`, which are density, velocities,
      // and pressure. Their relation is handled by the EquationOfState class. The
      // `conserved variables` are meant for internal calculation. Solving for
      // `conserved variables` guarantees the conservation properties. However,
      // `conserved variables` are not easy to use for diagnostics. Therefore, another
      // group of variables called the `primitive variables` are introduced for external
      // calculations, like calculating transport fluxes, radiative fluxes and
      // interacting with other physical components of the system. In this case, the
      // diffusion is calculated by using the `primitive variabes` and the result is
      // added to the `conserved variables` to ensure conservation properties.
      u(IM1,j,i) += dt*K*w(IDN,j,i)/(dx*dy)*(
        w(IVX,j,i-1) + w(IVX,j,i+1) + w(IVX,j-1,i) + w(IVX,j+1,i) - 4.*w(IVX,j,i));
      u(IM2,j,i) += dt*K*w(IDN,j,i)/(dx*dy)*(
        w(IVY,j,i-1) + w(IVY,j,i+1) + w(IVY,j-1,i) + w(IVY,j+1,i) - 4.*w(IVY,j,i));

      // Adding thermal dissipation is similar to adding viscous dissipation.
      u(IEN,j,i) += dt*K*w(IDN,j,i)/(dx*dy)*cp*temp/theta*(theta_ip1_j + theta_im1_j +
        theta_i_jp1 + theta_i_jm1 - 4.*theta);
    }
}

// @sect3{Program specific variables}

// This is the place where program specific variables are initialized.
// Not that the function is a member function of the Mesh class rather than the
// MeshBlock class we have been working with. The difference between class Mesh and
// class MeshBlock is that class Mesh is an all-encompassing class that manages multiple
// MeshBlocks while class MeshBlock manages all physics modules. During the
// instantiation of the classes. class Mesh is instantiated first and then it
// instantiates all MeshBlocks inside it. Therefore, this subroutine runs before any
// MeshBlock.
void Mesh::InitUserMeshData(ParameterInput *pin)
{
  // The program specific forcing parameters are set here.
  // They come from the input file, which is parsed by the ParameterInput class
  Real gamma = pin->GetReal("hydro", "gamma");
  K  = pin->GetReal("problem", "K");
  p0 = pin->GetReal("problem", "p0");
  Rd = pin->GetReal("thermodynamics", "Rd");
  cp = gamma/(gamma - 1.)*Rd;

  // This line code enrolls the forcing function we wrote in
  // section <a href="#Forcingfunction">Forcing function</a>
  EnrollUserExplicitSourceFunction(Diffusion);
}

// @sect3{Initial condition}

// This is the final part of the program, in which we set the initial condition.
// It is called the `problem generator` in a general Athena++ application.
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  // Get the value of gravity. Positive is upward and negative is downward.
  // We take the negative to get the absolute value.
  Real grav = -phydro->hsrc.GetG1();
  // These lines read the parameter values in the input file, which
  // are organized into sections. The problem specific parameters are usually placed
  // under section `problem`.
  Real Ts = pin->GetReal("problem", "Ts");
  Real xc = pin->GetReal("problem", "xc");
  Real xr = pin->GetReal("problem", "xr");
  Real zc = pin->GetReal("problem", "zc");
  Real zr = pin->GetReal("problem", "zr");
  Real dT = pin->GetReal("problem", "dT");

  // Loop over the grids. The purpose is to set the temperature, pressure and density
  // fields at each cell-centered grid.
  for (int j = js; j <= je; ++j) {
    for (int i = is; i <= ie; ++i) {
      // Get the Cartesian coordinates in the vertical and horizontal directions.
      // `x1` is usually the vertical direction and `x2` is the horizontal direction.
      // The meaning of `x1` and `x2` may change with the coordinate system.
      // <code>x1v</code> and <code>x2v</code> retrieve the coordinate at cell centers.
      Real x1 = pcoord->x1v(i);
      Real x2 = pcoord->x2v(j);
      // Distance to the center of a cold air bubble at (xc,zc).
      Real L = sqrt(sqr((x2 - xc)/xr) + sqr((x1 - zc)/zr));
      // Adiabatic temperature gradient.
      Real temp = Ts - grav*x1/cp;
      // Once we know the temperature, we can calculate the adiabatic and hydrostatic pressure as
      // $p=p_0(\frac{T}{T_s})^{c_p/R_d}$, where $p_0$ is the surface pressure and $T_s$ is the
      // surface temperature.
      phydro->w(IPR,j,i) = p0*pow(temp/Ts, cp/Rd);
      if (L <= 1.)
        temp += dT*(cos(M_PI*L) + 1.)/2.;
      // Set density using ideal gas law.
      phydro->w(IDN,j,i) = phydro->w(IPR,j,i)/(Rd*temp);
      // Initialize velocities to zero.
      phydro->w(IVX,j,i) = phydro->w(IVY,j,i) = 0.;
    }
  }

  // We have set all `primitive variables`. The last step is to calculate the `conserved
  // variables` based on the `primitive variables`. This is done by calling the member
  // function EquationOfState::PrimitiveToConserved. Note that <code>pfield->bcc</code>
  // is a placeholder because there is not MHD involved in this example.
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);
}
