// C/C++ headers
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

// external libraries
#include <application/application.hpp>
#include <application/exceptions.hpp>

// Athena++ headers
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// exchanger
#include <exchanger/exchanger.hpp>

// forcing
#include "forcing.hpp"

ForcingContainer ForcingFactory::CreateFrom(MeshBlock *pmb,
                                            ParameterInput *pin) {
  ForcingContainer forcing;

  char package_names[1024];
  std::string str = pin->GetOrAddString("forcing", "packages", "");
  std::strcpy(package_names, str.c_str());
  char *p = std::strtok(package_names, " ,");

  while (p != NULL) {
    if (std::strcmp(p, "relax_bot_temp") == 0) {
      auto nb = ExchangerHelper::find_bot_neighbor(pmb);
      if (nb == nullptr) {  // no bottom neighbor
        forcing.push_back(std::make_shared<RelaxBotTemp>(pmb, pin));
      }
    } else if (std::strcmp(p, "relax_bot_comp") == 0) {
      auto nb = ExchangerHelper::find_bot_neighbor(pmb);
      if (nb == nullptr) {  // no bottom neighbor
        forcing.push_back(std::make_shared<RelaxBotComp>(pmb, pin));
      }
    } else if (std::strcmp(p, "relax_bot_velo") == 0) {
      auto nb = ExchangerHelper::find_bot_neighbor(pmb);
      if (nb == nullptr) {  // no bottom neighbor
        forcing.push_back(std::make_shared<RelaxBotVelo>(pmb, pin));
      }
    } else if (std::strcmp(p, "top_sponge_lyr") == 0) {
      auto nb = ExchangerHelper::find_top_neighbor(pmb);
      if (nb == nullptr) {  // no top neighbor
        forcing.push_back(std::make_shared<TopSpongeLyr>(pmb, pin));
      }
    } else if (std::strcmp(p, "bot_sponge_lyr") == 0) {
      auto nb = ExchangerHelper::find_bot_neighbor(pmb);
      if (nb == nullptr) {  // no bottom neighbor
        forcing.push_back(std::make_shared<BotSpongeLyr>(pmb, pin));
      }
    } else if (std::strcmp(p, "left_sponge_lyr") == 0) {
      auto nb = ExchangerHelper::find_left_neighbor(pmb);
      if (nb == nullptr) {  // no left neighbor
        forcing.push_back(std::make_shared<LeftSpongeLyr>(pmb, pin));
      }
    } else if (std::strcmp(p, "right_sponge_lyr") == 0) {
      auto nb = ExchangerHelper::find_right_neighbor(pmb);
      if (nb == nullptr) {  // no right neighbor
        forcing.push_back(std::make_shared<RightSpongeLyr>(pmb, pin));
      }
    } else if (std::strcmp(p, "top_cooling") == 0) {
      auto nb = ExchangerHelper::find_top_neighbor(pmb);
      if (nb == nullptr) {  // no top neighbor
        forcing.push_back(std::make_shared<TopCooling>(pmb, pin));
      }
    } else if (std::strcmp(p, "bot_heating") == 0) {
      auto nb = ExchangerHelper::find_bot_neighbor(pmb);
      if (nb == nullptr) {  // no bottom neighbor
        forcing.push_back(std::make_shared<BotHeating>(pmb, pin));
      }
    } else if (std::strcmp(p, "body_heating") == 0) {
      forcing.push_back(std::make_shared<BodyHeating>(pmb, pin));
    } else {
      std::stringstream msg;
      msg << "Package '" << p << "' "
          << "unrecognized." << std::endl;
      throw RuntimeError("ForcingFactory::CreateFrom", msg.str());
    }

    p = std::strtok(NULL, " ,");
  }

  return forcing;
}
