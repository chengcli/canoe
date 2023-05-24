#include "inversion_helper.hpp"

#include <vector>

#include "inversion.hpp"

void gather_probability(std::vector<Inversion *> &fitq) {
  auto qlast = fitq[fitq.size() - 1];

  // replace the log probability by the last one
  for (auto q : fitq) {
    int nwalker = q->GetWalkers();

    for (int k = 0; k < nwalker; ++k) {
      q->SetLogProbability(k, qlast->GetLogProbability(k));
    }
  }
}
