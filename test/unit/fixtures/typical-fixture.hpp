/* Copyright (c) 2012-2016 Big Ladder Software. All rights reserved.
* See the LICENSE file for additional terms and conditions. */

#ifndef TYPICAL_FIXTURE_HPP_
#define TYPICAL_FIXTURE_HPP_

#include "base-fixture.hpp"

using namespace Kiva;

class TypicalFixture : public BaseFixture {
protected:
  void SetUp() {

    fnd.reductionStrategy = Foundation::RS_BOUNDARY;

    fnd.hasSlab = true;
    fnd.slab.emissivity = 0.8;


    Layer tempLayer;
    tempLayer.thickness = 0.24;
    tempLayer.material = soil;

    fnd.wall.layers.push_back(tempLayer);

    fnd.wall.heightAboveGrade = 0.0;
    fnd.wall.depthBelowSlab = 0.0;
    fnd.wall.interiorEmissivity = 0.0;
    fnd.wall.exteriorEmissivity = 0.0;
    fnd.wall.exteriorAbsorptivity = 0.0;

    fnd.numericalScheme = Foundation::NS_STEADY_STATE;

    bcs.localWindSpeed = 0;
    bcs.outdoorTemp = 283.15;
    bcs.indoorTemp = 303.15;


    outputMap[Surface::ST_SLAB_CORE] = {
      GroundOutput::OT_RATE
    };
  }

};

#endif /* TYPICAL_FIXTURE_HPP_ */
