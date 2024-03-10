// athena
#include <athena/athena.hpp>
#include <athena/mesh/mesh.hpp>

// forcing
#include "forcing.hpp"

BotForcing::BotForcing(MeshBlock *pmb, int nvar) {
  bot_data_.NewAthenaArray(nvar, pmb->ncells3, pmb->ncells2);
}

size_t BotForcing::RestartDataSizeInBytes() const {
  size_t size = 0;
  size += bot_data_.GetSizeInBytes();

  return size;
}

size_t BotForcing::DumpRestartData(char *pdst) const {
  std::memcpy(pdst, bot_data_.data(), bot_data_.GetSizeInBytes());
  return RestartDataSizeInBytes();
}

size_t BotForcing::LoadRestartData(char *psrc) {
  std::memcpy(bot_data_.data(), psrc, bot_data_.GetSizeInBytes());
  return RestartDataSizeInBytes();
}
