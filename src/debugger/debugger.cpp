// C/C++ headers
#include <algorithm>
#include <cstring>
#include <iomanip>
#include <memory>
#include <sstream>

// debugger headers
#include <configure.hpp>

#include "debugger.hpp"
// #include "../particles/material_point.hpp"

// MPI
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

char const Debugger::cgreen[] = "\033[0;32m";
char const Debugger::cend[] = "\033[0m";

// global debugger
std::unique_ptr<Debugger> pdebug;

Debugger::Debugger(int depth)
    : prev(nullptr),
      next(nullptr),
      fname_("HEAD"),
      depth_(depth),
      current_depth_(0) {
  // data_.NewAthenaArray(2*NHYDRO, pmb->ncells3, pmb->ncells2, pmb->ncells1);
  idstack_next_.push_back("1.");
}

Debugger::~Debugger() {
  if (prev != nullptr) prev->next = next;
  if (next != nullptr) next->prev = prev;
}

/*Debugger* Debugger::StartTracking(std::string name)
{
  Debugger *my = this;
  while (my->next != nullptr) my = my->next;
  my->next = new Debugger(pmy_block);
  my->next->fname_ = name;
  my->next->prev = my;
  my->next->next = nullptr;

  return my->next;
}

void Debugger::Track3D(std::string name, TestFunc_t test, AthenaArray<Real>&
var, int n)
{
  std::stringstream ss;
  bool pass = true;
  MeshBlock *pmb = pmy_block;

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;
  int c1, c2, c3;

  is -= NGHOST; ie += NGHOST;
  if (pmb->pmy_mesh->f2) {
    js -= NGHOST; je += NGHOST;
  }
  if (pmb->pmy_mesh->f3) {
    ks -= NGHOST; ke += NGHOST;
  }

  vnames_.push_back(name);

  if (vnames_.size() >= 2*NHYDRO) {
    ss << "### FATAL ERROR in function [Debugger::Track3D]"
       << std::endl << "Maximum number ("
       << 2*NHYDRO << ") of trakcers reached." << std::endl;
    ATHENA_ERROR(ss);
  }

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        data_(vnames_.size() - 1,k,j,i) = var(n,k,j,i);
        bool inside = (i >= pmb->is) && (i <= pmb->ie) && (j >= pmb->js) && (j
<= pmb->je)
                   && (k >= pmb->ks) && (k <= pmb->ke);
        if (!test(var(n,k,j,i)) && pass && inside) {
          pass = false;
          c1 = i;
          c2 = j;
          c3 = k;
        }
      }

  if (!pass) {
    DumpTracking(name, c1, c2, c3, "w");
    ss << "### FATAL ERROR in function [Debugger::Track3D]"
        << std::endl << "Invalid \033[0;32m" << name << "\033[0m at "
        << "(k=" << c3 << ", " << "j=" << c2 << ", " << "i=" << c1 << ") "
        << "on rank \033[0;32m" << Globals::my_rank << "\033[0m called from
function \033[0;32m"
        << fname_ << "\033[0m" << std::endl;
    ATHENA_ERROR(ss);
  }
}

void Debugger::Track1D(std::string name, TestFunc_t test, AthenaArray<Real>&
var, int n, int k, int j)
{
  std::stringstream ss;
  bool pass = true;
  MeshBlock *pmb = pmy_block;

  int is = pmb->is - NGHOST, ie = pmb->ie + NGHOST;
  int c1, c2, c3;

  std::vector<std::string>::iterator it = std::find(vnames_.begin(),
vnames_.end(), name); if (it == vnames_.end()) vnames_.push_back(name); int m =
std::find(vnames_.begin(), vnames_.end(), name) - vnames_.begin();

  if (m >= 2*NHYDRO) {
    ss << "### FATAL ERROR in function [Debugger::Track1D]"
        << std::endl << "Maximum number ("
        << 2*NHYDRO << ") of trakcers reached." << std::endl;
    ATHENA_ERROR(ss);
  }

  for (int i = is; i <= ie; ++i) {
    data_(m,k,j,i) = var(n,i);
    bool inside = (i >= pmb->is) && (i <= pmb->ie);
    if (!test(var(n,i)) && pass && inside) {
      pass = false;
      c1 = i;
      c2 = j;
      c3 = k;
    }
  }

  if (!pass) {
    DumpTracking(name, c1, c2, c3, "w");
    ss << "### FATAL ERROR in function [Debugger::Track1D]"
        << std::endl << "Invalid \033[0;32m" << name << "\033[0m at "
        << "(k=" << c3 << ", " << "j=" << c2 << ", " << "i=" << c1 << ") "
        << "on rank \033[0;32m" << Globals::my_rank << "\033[0m called from
function \033[0;32m"
        << fname_ << "\033[0m" << std::endl;
    ATHENA_ERROR(ss);
  }
}

void Debugger::DumpTracking(std::string name, int c1, int c2, int c3, char
const* mode)
{
  MeshBlock *pmb = pmy_block;
  Coordinates *pcoord = pmb->pcoord;

  int is = c1 - NGHOST, ie = c1 + NGHOST;
  int js, je, ks, ke;
  if (pmb->pmy_mesh->f2) {
    js = c2 - NGHOST;
    je = c2 + NGHOST;
  } else
    js = je = 0;
  if (pmb->pmy_mesh->f3) {
    ks = c3 - NGHOST;
    ke = c3 + NGHOST;
  } else
    ks = ke = 0;

  // loop over all variables
  for (int n = 0; n < vnames_.size(); ++n) {
    // open file for output
    std::string out;
    char blockid[12];
    std::snprintf(blockid, sizeof(blockid), "block%d", pmb->gid);
    out.assign(vnames_[n]);
    out.append(".");
    out.append(blockid);
    out.append(".txt");

    FILE *pfile;
    std::stringstream ss;
    if ((pfile = std::fopen(out.c_str(), mode)) == nullptr) {
      ss << "### FATAL ERROR in function [Debugger::DumpTracking]"
          << std::endl << "Output file '" << out << "' could not be opened"
          << std::endl;
      ATHENA_ERROR(ss);
    }

    // file header
    std::fprintf(pfile, "# Debugger dump at time = %e\n", pmb->pmy_mesh->time);
    std::fprintf(pfile, "# Variable = %s\n", vnames_[n].c_str());
    std::fprintf(pfile, "# Subdomain extent = [%.3g, %.3g] x [%.3g, %.3g] x
[%.3g, %.3g]\n", pcoord->x1f(pmb->is), pcoord->x1f(pmb->ie+1),
      pcoord->x2f(pmb->js), pcoord->x2f(pmb->je+1),
      pcoord->x3f(pmb->ks), pcoord->x3f(pmb->ke+1));
    std::fprintf(pfile, "# Break at (k=%d, j=%d, i=%d) of variable %s\n", c3,
c2, c1, name.c_str());

    // data table
    std::fprintf(pfile, "# Currently at [%s]\n", fname_.c_str());
    std::fprintf(pfile, "# Previously at [%s]\n", prev->fname_.c_str());
    for (int k = ks; k <= ke; ++k) {
      std::fprintf(pfile, "k=%02d        ", k);
      for (int j = js; j <= je; ++j)
        std::fprintf(pfile, " (k=%02d,j=%02d)", k, j);
      std::fprintf(pfile, "\n");

      for (int i = is; i <= ie; ++i) {
        std::fprintf(pfile, "i=%02d        ", i);
        for (int j = js; j <= je; ++j)
          std::fprintf(pfile, "%12.3g", data_(n,k,j,i));
        std::fprintf(pfile, "\n");
      }
      std::fprintf(pfile, "\n");
    }

    fclose(pfile);
  }

  if (prev != nullptr)
    prev->DumpTracking(name, c1, c2, c3, "a");
}*/

void Debugger::Enter(std::string name, std::string heil) {
  sections_.push_back(name);
  std::string id = idstack_next_.back();

  if (Globals::my_rank == 0 && current_depth_ < depth_) {
    std::cout << msg.str() << id << " " << heil << " " << name << " ..."
              << std::endl;
  }

  idstack_next_.push_back(id + "1.");
  msg.str("");
  current_depth_++;
}

void Debugger::Leave() {
  current_depth_--;
  std::stringstream ss;
  if (sections_.size() == 0) {
    Fatal("Debugger::Leave", "Caller stack is empty");
  }

  std::string name = sections_.back();
  sections_.pop_back();
  idstack_next_.pop_back();
  increment_id(&idstack_next_.back());
  if (Globals::my_rank == 0 && current_depth_ < depth_) {
    std::cout << msg.str() << name << " done." << std::endl;
  }
  msg.str("");
}

/*void Debugger::CheckConservation(std::string name, AthenaArray<Real> const&
var, int is, int ie, int js, int je, int ks, int ke) { int nvar = var.GetDim4();
  Real *sum = new Real [nvar];
  std::fill(sum, sum + nvar, 0.);

  for (int n = 0; n < nvar; ++n)
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i)
          sum[n] += var(n,k,j,i);

#ifdef MPI_PARALLEL
  if (Globals::my_rank == 0) {
    MPI_Reduce(MPI_IN_PLACE, sum, nvar, MPI_ATHENA_REAL, MPI_SUM, 0,
MPI_COMM_WORLD); } else { MPI_Reduce(sum, sum, nvar, MPI_ATHENA_REAL, MPI_SUM,
0, MPI_COMM_WORLD);
  }
#endif

  if (Globals::my_rank == 0) {
    //msg << cgreen << ">> Total " << name << "  = (";
    msg << "\t>> Total " << name << "  = (";
    for (int n = 0; n < nvar; ++n)
      msg << std::setprecision(12) << sum[n] << ", ";
    //msg << ")" << cend << std::endl;
    msg << ")" << std::endl;
  }
  delete[] sum;
}*/

/*void Debugger::CheckParticleConservation(std::vector<std::string> const&
cnames, std::vector<MaterialPoint> const& mp)
{
  Real *sum = new Real [cnames.size()];
  int *num = new int [cnames.size()];

  std::fill(sum, sum + cnames.size(), 0.);
  std::fill(num, num + cnames.size(), 0);

  for (std::vector<MaterialPoint>::const_iterator q = mp.begin(); q != mp.end();
++q) { sum[q->type] += q->rho; num[q->type] += 1;
  }

#ifdef MPI_PARALLEL
  if (Globals::my_rank == 0) {
    MPI_Reduce(MPI_IN_PLACE, sum, cnames.size(), MPI_ATHENA_REAL, MPI_SUM, 0,
MPI_COMM_WORLD); MPI_Reduce(MPI_IN_PLACE, num, cnames.size(), MPI_INT, MPI_SUM,
0, MPI_COMM_WORLD); } else { MPI_Reduce(sum, sum, cnames.size(),
MPI_ATHENA_REAL, MPI_SUM, 0, MPI_COMM_WORLD); MPI_Reduce(num, num,
cnames.size(), MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  }
#endif

  if (Globals::my_rank == 0) {
    //for (int i = 0; i < cnames.size(); ++i)
    //  msg << cgreen << ">> Total " << cnames[i] << " = " << sum[i] << ", num =
" << num[i] << cend << std::endl; for (int i = 0; i < cnames.size(); ++i) msg <<
"\t>> Total " << cnames[i] << " = " << sum[i] << ", num = " << num[i] <<
std::endl;
  }

  delete[] sum;
  delete[] num;
}*/

Debugger* Debugger::Message(std::string str) {
  msg << "- " << str << std::endl;
  return this;
}

void Debugger::Fatal(std::string where, std::string what) {
  std::cout << "### FATAL ERROR in " << where << std::endl;
#ifdef MPI_PARALLEL
  MPI_Finalize();
#endif
  throw std::runtime_error(what.c_str());
}

void Debugger::Fatal(std::string where, std::string str, std::string what) {
  std::cout << "### FATAL ERROR in " << where << " " << str << std::endl;
#ifdef MPI_PARALLEL
  MPI_Finalize();
#endif
  throw std::runtime_error(what.c_str());
}

void Debugger::Print(std::string str) {
  if (Globals::my_rank == 0) {
    std::cout << str << std::endl;
  }
}

void increment_id(std::string* str) {
  int len = str->size();
  int x, i = len - 1;
  while ((i > 0) && ((*str)[i - 1] != '.')) i--;
  x = std::stoi(str->substr(i, len - i - 1));
  x += 1;
  *str = str->substr(0, i) + std::to_string(x) + '.';
}
