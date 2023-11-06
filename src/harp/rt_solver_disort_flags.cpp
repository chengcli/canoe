// external
#include <yaml-cpp/yaml.h>

// canoe
#include <configure.hpp>

// harp
#include "rt_solvers.hpp"

#ifdef RT_DISORT

void RadiationBand::RTSolverDisort::setFlagsFromNode(YAML::Node const& my) {
  if (my["ibcnd"]) {
    ds_.flag.ibcnd = my["ibcnd"].as<bool>();
  } else {
    ds_.flag.ibcnd = 0;
  }

  if (my["usrtau"]) {
    ds_.flag.usrtau = my["usrtau"].as<bool>();
  } else {
    ds_.flag.usrtau = 1;
  }

  if (my["usrang"]) {
    ds_.flag.usrang = my["usrang"].as<bool>();
  } else {
    ds_.flag.usrang = 1;
  }

  if (my["lamber"]) {
    ds_.flag.lamber = my["lamber"].as<bool>();
  } else {
    ds_.flag.lamber = 0;
  }

  if (my["planck"]) {
    ds_.flag.planck = my["planck"].as<bool>();
  } else {
    ds_.flag.planck = 0;
  }

  if (my["spher"]) {
    ds_.flag.spher = my["spher"].as<bool>();
  } else {
    ds_.flag.spher = 0;
  }

  if (my["onlyfl"]) {
    ds_.flag.onlyfl = my["onlyfl"].as<bool>();
  } else {
    ds_.flag.onlyfl = 0;
  }

  if (my["quiet"]) {
    ds_.flag.quiet = my["quiet"].as<bool>();
  } else {
    ds_.flag.quiet = 0;
  }

  if (my["brdf_type"]) {
    ds_.flag.brdf_type = my["brdf_type"].as<int>();
  } else {
    ds_.flag.brdf_type = 0;
  }

  if (my["intensity_correction"]) {
    ds_.flag.intensity_correction = my["intensity_correction"].as<bool>();
  } else {
    ds_.flag.intensity_correction = 0;
  }

  if (my["old_intensity_correction"]) {
    ds_.flag.old_intensity_correction =
        my["old_intensity_correction"].as<bool>();
  } else {
    ds_.flag.old_intensity_correction = 0;
  }

  if (my["general_source"]) {
    ds_.flag.general_source = my["general_source"].as<bool>();
  } else {
    ds_.flag.general_source = 0;
  }

  if (my["output_uum"]) {
    ds_.flag.output_uum = my["output_uum"].as<bool>();
  } else {
    ds_.flag.output_uum = 0;
  }

  if (my["prnt"]) {
    std::vector<int> prnt = my["prnt"].as<std::vector<int>>();
    for (int i = 0; i < 5; ++i) {
      ds_.flag.prnt[i] = prnt[i];
    }
  } else {
    for (int i = 0; i < 5; ++i) {
      ds_.flag.prnt[i] = 0;
    }
  }
}

#endif  // RT_DISORT
