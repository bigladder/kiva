/* Copyright (c) 2012-2016 Big Ladder Software. All rights reserved.
* See the LICENSE file for additional terms and conditions. */

#ifndef BoundaryConditions_HPP
#define BoundaryConditions_HPP

namespace Kiva {

class BoundaryConditions {
public:
  double indoorTemp;
  double outdoorTemp;
  double localWindSpeed;
  double solarAzimuth;
  double solarAltitude;
  double directNormalFlux;
  double diffuseHorizontalFlux;
  double skyEmissivity;

  BoundaryConditions() :
    indoorTemp(293.15),
    outdoorTemp(273.15),
    localWindSpeed(0.0),
    solarAzimuth(3.14),
    solarAltitude(0.0),
    directNormalFlux(0.0),
    diffuseHorizontalFlux(0.0),
    skyEmissivity(0.0)
  {}

};

}
#endif // BoundaryConditions_HPP
