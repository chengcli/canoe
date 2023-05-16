// C/C++ headers
#include <cstring>

// Athena++ headers
#include <parameter_input.hpp>
#include <mesh/mesh.hpp>

// harp2 headers
#include "../debugger/debugger.hpp"
#include "../utils/sentinelq.hpp"
#include "../utils/vectorize.hpp"
#include "../mesh/meshblock_impl.hpp"
#include "../mesh/block_index.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "profile_inversion.hpp"

void new_inversion_queue(SentinelQ<Inversion*> &fitq,
  MeshBlock *pmb, ParameterInput *pin,
  BlockIndex *pblock, Thermodynamics *pthermo)
{
  std::string str = pin->GetOrAddString("inversion", "tasks", "");
  std::vector<std::string> task_names = Vectorize<std::string>(str.c_str(), " ,");

  Inversion *pfit;
  for (auto p : task_names) {
    if (p == "VLAProfileInversion") {
      pfit = new VLAProfileInversion(pmb, pin);
    } else if (p == "JunoProfileInversion") {
      pfit = new JunoProfileInversion(pmb, pin);
    } else if (p == "VLACompositionInversion") {
    } else if (p == "JunoCompositionInversion") {
    } else {
      Debugger::Fatal("new_inversion_queue", "task::" + p, "unrecognized");
    }
    fitq.push(pfit);
  }

  auto q = fitq.getNext();
  int jl = pmb->js;
  while (q != nullptr) {
    pfit = q->getData();
    pfit->use(pblock)->use(pthermo);
    pfit->InitializePositions();
    pfit->setX2Indices(jl);
    jl += pfit->getX2Span();
    q = q->getNext();
  }
}
