/* Copyright (c) 2012-2016 Big Ladder Software. All rights reserved.
* See the LICENSE file for additional terms and conditions. */

#include "fixtures/bestest-fixture.hpp"
//#include "fixtures/typical-fixture.hpp"

using namespace Kiva;

TEST_F( BESTESTFixture, GC10a)
{
  fnd.wallTopBoundary = Foundation::WTB_LINEAR_DT;
  fnd.wallTopInteriorTemperature = 303.15;
  fnd.wallTopExteriorTemperature = 283.15;
  EXPECT_NEAR(calcQ(), 2435, 55);
}

TEST_F( BESTESTFixture, GC30a)
{
  fnd.deepGroundDepth = 30.0;
  fnd.farFieldWidth = 20.0;
  EXPECT_NEAR(calcQ(), 2650, 100);
}

// Not an actual BESTEST Case, though recommended by TRNSYS report
TEST_F( BESTESTFixture, 1D)
{
  fnd.exposedFraction = 0.0;
  fnd.deepGroundDepth = 1.0;
  fnd.useDetailedExposedPerimeter = false;

  double area = 144; // m2
  double expectedQ = fnd.soil.conductivity/fnd.deepGroundDepth*area
    *(bcs.indoorTemp - fnd.deepGroundTemperature);
  EXPECT_NEAR(calcQ(), expectedQ, 1.0);
}


// Google Test main
int
main( int argc, char **argv )
{
	::testing::InitGoogleTest( &argc, argv );
	return RUN_ALL_TESTS();
}
