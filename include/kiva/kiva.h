/* Copyright (c) 2012-2022 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#pragma once

#include "boundary_conditions.h"
#include "foundation.h"

namespace Kiva {

class Kiva {
public:
  void calculate(const BoundaryConditions &boundary_conditions, const double timestep = 0.0);
};

class Aggregator {
public:
};

} // namespace Kiva
