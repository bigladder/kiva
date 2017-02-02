/* Copyright (c) 2012-2016 Big Ladder Software. All rights reserved.
* See the LICENSE file for additional terms and conditions. */

#include <gtest/gtest.h>
#include "fixtures/bestest.hpp"

using namespace Kiva;

TEST_F( BESTEST, GC10a)
{
  fnd.wallTopBoundary = Foundation::WTB_LINEAR_DT;
  fnd.wallTopInteriorTemperature = 303.15;
  fnd.wallTopExteriorTemperature = 283.15;
  EXPECT_NEAR(calcQ(), 2435, 50);
}

TEST_F( BESTEST, GC30a)
{
  fnd.deepGroundDepth = 30.0;
  fnd.farFieldWidth = 20.0;
  EXPECT_NEAR(calcQ(), 2650, 100);
}
