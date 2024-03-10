// C/C++
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

// canoe
#include <configure.hpp>
#include <impl.hpp>

// athena
#include <athena/athena.hpp>
#include <athena/coordinates/coordinates.hpp>
#include <athena/globals.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/outputs/user_outputs.hpp>

// math
#include <climath/core.h>

// harp
#include <harp/radiation.hpp>
#include <harp/radiation_band.hpp>

// utils
#include <utils/vectorize.hpp>

// outputs
#include "output_utils.hpp"

// Only proceed if PNETCDF output enabled
#ifdef PNETCDFOUTPUT

// External library headers
#include <mpi.h>
#include <pnetcdf.h>

#define ERR                                               \
  {                                                       \
    if (err != NC_NOERR) {                                \
      printf("Error at %s:%d : %s\n", __FILE__, __LINE__, \
             ncmpi_strerror(err));                        \
    }                                                     \
  }

//----------------------------------------------------------------------------------------
// PnetcdfOutput constructor
// destructor - not needed for this derived class

PnetcdfOutput::PnetcdfOutput(OutputParameters oparams) : OutputType(oparams) {}

//----------------------------------------------------------------------------------------
//! \fn void PnetcdfOutput:::WriteOutputFile(Mesh *pm)
//  \brief Cycles over all MeshBlocks and writes OutputData in PNETCDF format
//         One timestep per file

void PnetcdfOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag) {
  // create filename: "file_basename"+"."+"file_id"+"."+XXXXX+".nc",
  // where XXXXX = 5-digit file_number
  auto pmeta = MetadataTable::GetInstance();

  std::string fname;
  char number[6];
  int err;
  snprintf(number, sizeof(number), "%05d", output_params.file_number);

  fname.assign(output_params.file_basename);
  fname.append(".");
  fname.append(output_params.file_id);
  fname.append(".");
  fname.append(number);
  fname.append(".nc");

  // 0. reference radius for spherical polar geometry
  float radius;
  if (strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    try {
      radius = pin->GetReal("problem", "radius");
    } catch (std::runtime_error &e) {
      std::stringstream msg;
      msg << "### FATAL ERROR in PnetcdfOutput::WriteOutputFile" << std::endl
          << "Spherical polar coordinate system must define a reference radius "
             "in section <problem>"
          << std::endl
          << "Use radius = XXX " << std::endl;
      ATHENA_ERROR(msg);
    }
  }

  // 1. open file for output
  int ifile;
  err = ncmpi_create(MPI_COMM_WORLD, fname.c_str(), NC_CLOBBER, MPI_INFO_NULL,
                     &ifile);
  ERR;

  // 2. coordinate structure
  size_t nx1 = pm->mesh_size.nx1;
  size_t nx2 = pm->mesh_size.nx2;
  size_t nx3 = pm->mesh_size.nx3;
  size_t nx1f = nx1;
  if (nx1 > 1) nx1f++;
  size_t nx2f = nx2;
  if (nx2 > 1) nx2f++;
  size_t nx3f = nx3;
  if (nx3 > 1) nx3f++;
  //! \todo This applies to the first block. Does it work for all blocks?

  size_t nrays = pm->my_blocks(0)->pimpl->prad->GetNumOutgoingRays();

  // 2. define coordinate
  int idt, idx1, idx2, idx3, idx1f, idx2f, idx3f, iray;
  int ndims = 0;

  ncmpi_def_dim(ifile, "time", NC_UNLIMITED, &idt);

  ncmpi_def_dim(ifile, "x1", nx1, &idx1);
  ndims++;
  if (nx1 > 1) {
    ncmpi_def_dim(ifile, "x1f", nx1f, &idx1f);
    ndims++;
  }

  ncmpi_def_dim(ifile, "x2", nx2, &idx2);
  ndims++;
  if (nx2 > 1) {
    ncmpi_def_dim(ifile, "x2f", nx2f, &idx2f);
    ndims++;
  }

  ncmpi_def_dim(ifile, "x3", nx3, &idx3);
  ndims++;
  if (nx3 > 1) {
    ncmpi_def_dim(ifile, "x3f", nx3f, &idx3f);
    ndims++;
  }

  if (nrays > 0) {
    ncmpi_def_dim(ifile, "ray", nrays, &iray);
    ndims += 2;
  }

  // 3. define variables
  int ivt, ivx1, ivx2, ivx3, ivx1f, ivx2f, ivx3f, imu, iphi;

  ncmpi_def_var(ifile, "time", NC_FLOAT, 1, &idt, &ivt);
  ncmpi_put_att_text(ifile, ivt, "axis", 1, "T");
  ncmpi_put_att_text(ifile, ivt, "long_name", 4, "time");
  if (strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    ncmpi_put_att_text(ifile, ivt, "units", 25, "seconds since 1-1-1 0:0:0");
  } else {
    ncmpi_put_att_text(ifile, ivt, "units", 7, "seconds");
  }

  ncmpi_def_var(ifile, "x1", NC_FLOAT, 1, &idx1, &ivx1);
  ncmpi_put_att_text(ifile, ivx1, "axis", 1, "Z");
  if (strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    ncmpi_put_att_text(ifile, ivx1, "units", 6, "meters");
    ncmpi_put_att_text(ifile, ivx1, "long_name", 8, "altitude");
  } else {
    ncmpi_put_att_text(ifile, ivx1, "units", 6, "meters");
    ncmpi_put_att_text(ifile, ivx1, "long_name", 27,
                       "Z-coordinate at cell center");
  }
  if (nx1 > 1) {
    ncmpi_def_var(ifile, "x1f", NC_FLOAT, 1, &idx1f, &ivx1f);
    ncmpi_put_att_text(ifile, ivx1f, "axis", 1, "Z");
    if (strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
      ncmpi_put_att_text(ifile, ivx1f, "units", 6, "meters");
      ncmpi_put_att_text(ifile, ivx1f, "long_name", 8, "altitude");
    } else {
      ncmpi_put_att_text(ifile, ivx1f, "units", 6, "meters");
      ncmpi_put_att_text(ifile, ivx1f, "long_name", 25,
                         "Z-coordinate at cell face");
    }
  }

  ncmpi_def_var(ifile, "x2", NC_FLOAT, 1, &idx2, &ivx2);
  ncmpi_put_att_text(ifile, ivx2, "axis", 1, "Y");
  if (strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    ncmpi_put_att_text(ifile, ivx2, "units", 13, "degrees_north");
    ncmpi_put_att_text(ifile, ivx2, "long_name", 8, "latitude");
  } else {
    ncmpi_put_att_text(ifile, ivx2, "units", 6, "meters");
    ncmpi_put_att_text(ifile, ivx2, "long_name", 27,
                       "Y-coordinate at cell center");
  }
  if (nx2 > 1) {
    ncmpi_def_var(ifile, "x2f", NC_FLOAT, 1, &idx2f, &ivx2f);
    ncmpi_put_att_text(ifile, ivx2f, "axis", 1, "Y");
    if (strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
      ncmpi_put_att_text(ifile, ivx2f, "units", 13, "degrees_north");
      ncmpi_put_att_text(ifile, ivx2f, "long_name", 8, "latitude");
    } else {
      ncmpi_put_att_text(ifile, ivx2f, "units", 6, "meters");
      ncmpi_put_att_text(ifile, ivx2f, "long_name", 25,
                         "Y-coordinate at cell face");
    }
  }

  ncmpi_def_var(ifile, "x3", NC_FLOAT, 1, &idx3, &ivx3);
  ncmpi_put_att_text(ifile, ivx3, "axis", 1, "X");
  if (strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    ncmpi_put_att_text(ifile, ivx3, "units", 12, "degrees_east");
    ncmpi_put_att_text(ifile, ivx3, "long_name", 9, "longitude");
  } else {
    ncmpi_put_att_text(ifile, ivx3, "units", 6, "meters");
    ncmpi_put_att_text(ifile, ivx3, "long_name", 27,
                       "X-coordinate at cell center");
  }
  if (nx3 > 1) {
    ncmpi_def_var(ifile, "x3f", NC_FLOAT, 1, &idx3f, &ivx3f);
    ncmpi_put_att_text(ifile, ivx3f, "axis", 1, "X");
    if (strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
      ncmpi_put_att_text(ifile, ivx3f, "units", 12, "degrees_east");
      ncmpi_put_att_text(ifile, ivx3f, "long_name", 9, "longitude");
    } else {
      ncmpi_put_att_text(ifile, ivx3f, "units", 6, "meters");
      ncmpi_put_att_text(ifile, ivx3f, "long_name", 25,
                         "X-coordinate at cell face");
    }
  }

  if (nrays > 0) {
    ncmpi_def_var(ifile, "mu_out", NC_FLOAT, 1, &iray, &imu);
    ncmpi_put_att_text(ifile, imu, "units", 1, "1");
    ncmpi_put_att_text(ifile, imu, "long_name", 18, "cosine polar angle");
    ncmpi_def_var(ifile, "phi_out", NC_FLOAT, 1, &iray, &iphi);
    ncmpi_put_att_text(ifile, iphi, "units", 3, "rad");
    ncmpi_put_att_text(ifile, iphi, "long_name", 15, "azimuthal angle");
  }

  MeshBlock *pmb = pm->my_blocks(0);
  LoadOutputData(pmb);
  OutputData *pdata = pfirst_data_;

  // count total variables (vector variables are expanded into flat scalars)
  int total_vars = 0;
  while (pdata != nullptr) {
    auto names = Vectorize<std::string>(pdata->name.c_str(), ",");
    std::string grid = pmeta->GetGridType(names[0]);
    int nvar = get_num_variables(grid, pdata->data);

    total_vars += nvar;
    pdata = pdata->pnext;
  }

  int iaxis[4] = {idt, idx1, idx2, idx3};
  int iaxis1[4] = {idt, idx1f, idx2, idx3};
  int iaxis2[4] = {idt, idx1, idx2f, idx3};
  int iaxis3[4] = {idt, idx1, idx2, idx3f};
  int iaxisr[4] = {idt, iray, idx2, idx3};
  int *var_ids = new int[total_vars];
  int *ivar = var_ids;

  pdata = pfirst_data_;
  while (pdata != nullptr) {
    auto names = Vectorize<std::string>(pdata->name.c_str(), ",");
    std::string grid = pmeta->GetGridType(names[0]);
    int nvar = get_num_variables(grid, pdata->data);

    std::vector<std::string> varnames;
    if (names.size() >= nvar) {
      for (int n = 0; n < nvar; ++n) {
        varnames.push_back(names[n]);
      }
    } else {
      for (int n = 0; n < nvar; ++n) {
        size_t pos = pdata->name.find('?');
        if (nvar == 1) {                     // SCALARS
          if (pos < pdata->name.length()) {  // find '?'
            varnames.push_back(pdata->name.substr(0, pos) +
                               pdata->name.substr(pos + 1));
          } else {
            varnames.push_back(pdata->name);
          }
        } else {  // VECTORS
          char c[16];
          snprintf(c, sizeof(c), "%d", n + 1);
          if (pos < pdata->name.length()) {  // find '?'
            varnames.push_back(pdata->name.substr(0, pos) + c +
                               pdata->name.substr(pos + 1));
          } else {
            varnames.push_back(pdata->name + c);
          }
        }
      }
    }

    for (int n = 0; n < nvar; ++n) {
      auto name = varnames[n];

      if (grid == "RCC")  // radiation rays
        ncmpi_def_var(ifile, name.c_str(), NC_FLOAT, 4, iaxisr, ivar);
      else if (grid == "CCF")
        ncmpi_def_var(ifile, name.c_str(), NC_FLOAT, 4, iaxis1, ivar);
      else if ((grid == "CFC") && (nx2 > 1))
        ncmpi_def_var(ifile, name.c_str(), NC_FLOAT, 4, iaxis2, ivar);
      else if ((grid == "FCC") && (nx3 > 1))
        ncmpi_def_var(ifile, name.c_str(), NC_FLOAT, 4, iaxis3, ivar);
      else if (grid == "--C")
        ncmpi_def_var(ifile, name.c_str(), NC_FLOAT, 2, iaxis, ivar);
      else if (grid == "--F")
        ncmpi_def_var(ifile, name.c_str(), NC_FLOAT, 2, iaxis1, ivar);
      else if (grid == "---")
        ncmpi_def_var(ifile, name.c_str(), NC_FLOAT, 1, iaxis, ivar);
      else
        ncmpi_def_var(ifile, name.c_str(), NC_FLOAT, 4, iaxis, ivar);

      // set units
      auto attr = pmeta->GetUnits(pdata->name);
      if (attr != "") {
        ncmpi_put_att_text(ifile, *ivar, "units", attr.length(), attr.c_str());
      }

      // set long_name
      attr = pmeta->GetLongName(pdata->name);
      if (attr != "") {
        ncmpi_put_att_text(ifile, *ivar, "long_name", attr.length(),
                           attr.c_str());
      }

      ivar++;
    }
    pdata = pdata->pnext;
  }

  err = ncmpi_put_att_text(ifile, NC_GLOBAL, "Conventions", 6, "COARDS");
  ERR;
  if (strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    err = ncmpi_put_att_float(ifile, NC_GLOBAL, "PlanetRadius", NC_FLOAT, 1,
                              &radius);
    ERR;
  }

  err = ncmpi_enddef(ifile);
  ERR;

  // 4. allocate data buffer, MPI requests and status
  int max_ncells = 0, nmb = 0;
  for (int b = 0; b < pm->nblocal; ++b) {
    MeshBlock *pmb = pm->my_blocks(b);
    int nf1 = pmb->ie - pmb->is + 2 * NGHOST;
    int nf2 = pmb->je - pmb->js + 2 * NGHOST;
    int nf3 = pmb->ke - pmb->ks + 2 * NGHOST;
    max_ncells = std::max(max_ncells, nf1 * nf2 * nf3);
    nmb++;  // number of meshblocks this rank
  }
  int nbufs = nmb * (ndims + total_vars);
  int *reqs = new int[nbufs];
  int *stts = new int[nbufs];
  float **buf = new float *[nbufs];
  buf[0] = new float[nbufs * max_ncells];
  for (int i = 0; i < nbufs; ++i) buf[i] = buf[0] + i * max_ncells;
  int *ir = reqs;
  float **ib = buf;

  // 5. first meshblock writes time
  float time = (float)pm->time;
  MPI_Offset stime = 0, etime = 1;
  err = ncmpi_put_vara_float_all(ifile, ivt, &stime, &etime, &time);
  ERR;

  ClearOutputData();  // required when LoadOutputData() is used.

  // Loop over MeshBlocks
  // do {
  for (int b = 0; b < pm->nblocal; ++b) {
    MeshBlock *pmb = pm->my_blocks(b);
    // set start/end array indices depending on whether ghost zones are included
    out_is = pmb->is;
    out_ie = pmb->ie;
    out_js = pmb->js;
    out_je = pmb->je;
    out_ks = pmb->ks;
    out_ke = pmb->ke;

    // \todo include_ghost zones probably doesn't work with grids other than CCC
    if (output_params.include_ghost_zones) {
      out_is -= NGHOST;
      out_ie += NGHOST;
      if (out_js != out_je) {
        out_js -= NGHOST;
        out_je += NGHOST;
      }
      if (out_ks != out_ke) {
        out_ks -= NGHOST;
        out_ke += NGHOST;
      }
    }

    LoadOutputData(pmb);

    // 6. each meshblock writes coordinates and variables
    MPI_Offset ncells1 = out_ie - out_is + 1;
    MPI_Offset ncells2 = out_je - out_js + 1;
    MPI_Offset ncells3 = out_ke - out_ks + 1;
    MPI_Offset nfaces1 = ncells1;
    if (ncells1 > 1) nfaces1++;
    MPI_Offset nfaces2 = ncells2;
    if (ncells2 > 1) nfaces2++;
    MPI_Offset nfaces3 = ncells3;
    if (ncells3 > 1) nfaces3++;

    MPI_Offset start[4] = {0, ncells1 * pmb->loc.lx1, ncells2 * pmb->loc.lx2,
                           ncells3 * pmb->loc.lx3};
    MPI_Offset count[4] = {1, ncells1, ncells2, ncells3};
    MPI_Offset count1[4] = {1, nfaces1, ncells2, ncells3};
    MPI_Offset count2[4] = {1, ncells1, nfaces2, ncells3};
    MPI_Offset count3[4] = {1, ncells1, ncells2, nfaces3};
    MPI_Offset startr[4] = {0, 0, ncells2 * pmb->loc.lx2,
                            ncells3 * pmb->loc.lx3};

    MPI_Offset countr[4] = {1, (MPI_Offset)nrays, ncells2, ncells3};

    // MPI_Offset start_ax1[2] = {0, ncells1*pmb->loc.lx1};
    // MPI_Offset count_ax1[2] = {1, ncells1};
    // MPI_Offset count_ax1f[2] = {1, nfaces1};

    if (strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
      for (int i = out_is; i <= out_ie; ++i)
        (*ib)[i - out_is] = (float)(pmb->pcoord->x1v(i)) - radius;
    } else {
      for (int i = out_is; i <= out_ie; ++i)
        (*ib)[i - out_is] = (float)(pmb->pcoord->x1v(i));
    }
    err = ncmpi_iput_vara_float(ifile, ivx1, start + 1, count + 1, *ib++, ir++);
    ERR;

    if (nx1 > 1) {
      if (strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
        for (int i = out_is; i <= out_ie + 1; ++i)
          (*ib)[i - out_is] = (float)(pmb->pcoord->x1f(i)) - radius;
      } else {
        for (int i = out_is; i <= out_ie + 1; ++i)
          (*ib)[i - out_is] = (float)(pmb->pcoord->x1f(i));
      }
      err = ncmpi_iput_vara_float(ifile, ivx1f, start + 1, count1 + 1, *ib++,
                                  ir++);
      ERR;
    }

    if (strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
      for (int j = out_js; j <= out_je; ++j)
        (*ib)[j - out_js] = 90. - (float)rad2deg(pmb->pcoord->x2v(j));
    } else {
      for (int j = out_js; j <= out_je; ++j)
        (*ib)[j - out_js] = (float)(pmb->pcoord->x2v(j));
    }
    err = ncmpi_iput_vara_float(ifile, ivx2, start + 2, count + 2, *ib++, ir++);
    ERR;

    if (nx2 > 1) {
      if (strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
        for (int j = out_js; j <= out_je + 1; ++j)
          (*ib)[j - out_js] = 90. - (float)rad2deg(pmb->pcoord->x2f(j));
      } else {
        for (int j = out_js; j <= out_je + 1; ++j)
          (*ib)[j - out_js] = (float)(pmb->pcoord->x2f(j));
      }
      err = ncmpi_iput_vara_float(ifile, ivx2f, start + 2, count2 + 2, *ib++,
                                  ir++);
      ERR;
    }

    if (strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
      for (int k = out_ks; k <= out_ke; ++k)
        (*ib)[k - out_ks] = (float)rad2deg(pmb->pcoord->x3v(k));
    } else {
      for (int k = out_ks; k <= out_ke; ++k)
        (*ib)[k - out_ks] = (float)(pmb->pcoord->x3v(k));
    }
    err = ncmpi_iput_vara_float(ifile, ivx3, start + 3, count + 3, *ib++, ir++);
    ERR;

    if (nx3 > 1) {
      if (strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
        for (int k = out_ks; k <= out_ke + 1; ++k)
          (*ib)[k - out_ks] = (float)rad2deg(pmb->pcoord->x3f(k));
      } else {
        for (int k = out_ks; k <= out_ke + 1; ++k)
          (*ib)[k - out_ks] = (float)(pmb->pcoord->x3f(k));
      }
      err = ncmpi_iput_vara_float(ifile, ivx3f, start + 3, count3 + 3, *ib++,
                                  ir++);
      ERR;
    }

    if (nrays > 0) {
      auto prad = pmb->pimpl->prad;
      int m = 0;
      for (int b = 0; b < prad->GetNumBands(); ++b) {
        auto p = prad->GetBand(b);
        for (int n = 0; n < p->GetNumOutgoingRays(); ++n)
          (*ib)[m++] = (float)(p->GetCosinePolarAngle(n));
      }
      err = ncmpi_iput_var_float(ifile, imu, *ib++, ir++);
      ERR;

      m = 0;
      for (int b = 0; b < prad->GetNumBands(); ++b) {
        auto p = prad->GetBand(b);
        for (int n = 0; n < p->GetNumOutgoingRays(); ++n)
          (*ib)[m++] = (float)(p->GetAzimuthalAngle(n));
      }
      err = ncmpi_iput_var_float(ifile, iphi, *ib++, ir++);
      ERR;
    }

    ivar = var_ids;
    pdata = pfirst_data_;
    while (pdata != nullptr) {
      auto names = Vectorize<std::string>(pdata->name.c_str(), ",");
      std::string grid = pmeta->GetGridType(names[0]);
      int nvar = get_num_variables(grid, pdata->data);

      if (grid == "RCC") {  // radiation rays
        for (int n = 0; n < nvar; n++) {
          float *it = *ib;
          for (int m = 0; m < nrays; ++m)
            for (int j = out_js; j <= out_je; ++j)
              for (int k = out_ks; k <= out_ke; ++k)
                *it++ = (float)pdata->data(n, m, k, j);
          err = ncmpi_iput_vara_float(ifile, *ivar++, startr, countr, *ib++,
                                      ir++);
          ERR;
        }
      } else if (grid == "CCF") {
        for (int n = 0; n < nvar; n++) {
          float *it = *ib;
          for (int i = out_is; i <= out_ie + 1; ++i)
            for (int j = out_js; j <= out_je; ++j)
              for (int k = out_ks; k <= out_ke; ++k)
                *it++ = (float)pdata->data(n, k, j, i);
          err =
              ncmpi_iput_vara_float(ifile, *ivar++, start, count1, *ib++, ir++);
          ERR;
        }
      } else if ((grid == "CFC") && (nx2 > 1)) {
        for (int n = 0; n < nvar; n++) {
          float *it = *ib;
          for (int i = out_is; i <= out_ie; ++i)
            for (int j = out_js; j <= out_je + 1; ++j)
              for (int k = out_ks; k <= out_ke; ++k)
                *it++ = (float)pdata->data(n, k, j, i);
          err =
              ncmpi_iput_vara_float(ifile, *ivar++, start, count2, *ib++, ir++);
          ERR;
        }
      } else if ((grid == "FCC") && (nx3 > 1)) {
        for (int n = 0; n < nvar; n++) {
          float *it = *ib;
          for (int i = out_is; i <= out_ie; ++i)
            for (int j = out_js; j <= out_je; ++j)
              for (int k = out_ks; k <= out_ke + 1; ++k)
                *it++ = (float)pdata->data(n, k, j, i);
          err =
              ncmpi_iput_vara_float(ifile, *ivar++, start, count3, *ib++, ir++);
          ERR;
        }
      } else if (grid == "--C") {
        for (int n = 0; n < nvar; n++) {
          float *it = *ib;
          for (int i = out_is; i <= out_ie; ++i)
            *it++ = (float)pdata->data(n, i);
          err =
              ncmpi_iput_vara_float(ifile, *ivar++, start, count, *ib++, ir++);
          ERR;
        }
      } else if (grid == "--F") {
        for (int n = 0; n < nvar; n++) {
          float *it = *ib;
          for (int i = out_is; i <= out_ie + 1; ++i)
            *it++ = (float)pdata->data(n, i);
          err =
              ncmpi_iput_vara_float(ifile, *ivar++, start, count1, *ib++, ir++);
          ERR;
        }
      } else if (grid == "---") {
        for (int n = 0; n < nvar; n++) {
          **ib = (float)pdata->data(n, 0);
          err =
              ncmpi_iput_vara_float(ifile, *ivar++, start, count, *ib++, ir++);
          ERR;
        }
      } else {
        for (int n = 0; n < nvar; n++) {
          float *it = *ib;
          for (int i = out_is; i <= out_ie; ++i)
            for (int j = out_js; j <= out_je; ++j)
              for (int k = out_ks; k <= out_ke; ++k)
                *it++ = (float)pdata->data(n, k, j, i);
          err =
              ncmpi_iput_vara_float(ifile, *ivar++, start, count, *ib++, ir++);
          ERR;
        }
      }

      pdata = pdata->pnext;
    }

    ClearOutputData();  // required when LoadOutputData() is used.
  }
  //} while (LoadOutputData(pmb));   // end loop over MeshBLocks

  // 7. wait for all writings to complete
  ncmpi_wait_all(ifile, nbufs, reqs, stts);
  for (int i = 0; i < nbufs; ++i)
    if (stts[i] != NC_NOERR) {
      std::stringstream msg;
      msg << "### FATAL ERROR in PnetcdfOutput::WriteOutputFile" << std::endl
          << "Nonblocking write fails on request " << i << std::endl
          << ncmpi_strerror(stts[i]);
      ATHENA_ERROR(msg);
    }

  // 8. close nc file
  ncmpi_close(ifile);

  // 9. clear data buffers
  delete[] buf[0];
  delete[] buf;
  delete[] var_ids;
  delete[] reqs;
  delete[] stts;

  // 10. increment counters
  output_params.file_number++;
  output_params.next_time += output_params.dt;
  pin->SetInteger(output_params.block_name, "file_number",
                  output_params.file_number);
  pin->SetReal(output_params.block_name, "next_time", output_params.next_time);

  return;
}

#endif  // PNETCDFOUTPUT
