/* Copyright (c) 2012-2022 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#pragma once

#include <functional>

namespace Kiva {

using NaturalConvectionAlgorithm = std::function<double(double, double, double, double, double)>;

using ForcedConvectionTerm = std::function<double(double, double, double, double)>;

typedef std::function<double(double, double, double, double)> ForcedConvectionTerm;

struct SurfaceBoundaryConditions {
  double convective_temperature{293.15};
  double radiant_temperature{293.15};
  double absorbed_radiation{0.};
  NaturalConvectionAlgorithm natural_convection_algorithm;
};

struct CombinedConvectionAlgorithm {
  NaturalConvectionAlgorithm natural_convection_algorithm;
  ForcedConvectionTerm forced_convection_term;
};

struct BoundaryConditions {
  double outdoor_temperature{273.15};
  double local_wind_speed{0.};
  double wind_direction{0.};
  double solar_azimuth{3.14};
  double solar_altitude{0.};
  double direct_normal_solar_flux{0.};
  double diffuse_horizontal_solar_flux{0.};
  double sky_emissivity{0.8};
  double deep_ground_temperature{283.15};
  SurfaceBoundaryConditions slab_interior;
  SurfaceBoundaryConditions wall_interior;
  CombinedConvectionAlgorithm wall_exterior;
  CombinedConvectionAlgorithm grade;
};

} // namespace Kiva
