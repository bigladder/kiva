/* Copyright (c) 2012-2016 Big Ladder Software. All rights reserved.
* See the LICENSE file for additional terms and conditions. */

#ifndef BASE_FIXTURE_HPP_
#define BASE_FIXTURE_HPP_

#include <gtest/gtest.h>
#include "Ground.hpp"

using namespace Kiva;

class BaseFixture : public testing::Test {
protected:
  double calcQ(){
    Ground ground(fnd,outputMap);
    ground.buildDomain();
    ground.calculate(bcs);
    ground.calculateSurfaceAverages();
    return ground.getSurfaceAverageValue({Surface::ST_SLAB_CORE,GroundOutput::OT_RATE});

  }
  std::map<Surface::SurfaceType, std::vector<GroundOutput::OutputType>> outputMap;
  BoundaryConditions bcs;
  Foundation fnd;
};

#endif /* BASE_FIXTURE_HPP_ */
