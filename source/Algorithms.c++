/* Algorithms.c++ is part of Kiva (Written by Neal Kruis)
 * Copyright (C) 2012-2015 Big Ladder Software <info@bigladdersoftware.com>
 * All rights reserved.
 *
 * Kiva is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Kiva is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Kiva.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ConvectionAlgorithms_CPP
#define ConvectionAlgorithms_CPP

#include "Algorithms.h"

static const double PI = 4.0*atan(1.0);
static const double SIGMA = 5.67*pow(10.0,-8);  // [W/m2-K4]

double getDOE2ConvectionCoeff(double tilt,
                  double azimuth,
                  double windDirection,
                  double Tsurf,
                  double Tamb,
                  double Vair,
                  double roughness)
{
  /* Based on the DOE-2 convection model as used in EnergyPlus
   *
   * Roughness factors:
   * Very Rough = 2.17
   * Rough = 1.67
   * Medium Rough = 1.52
   * Medium Smooth = 1.13
   * Smooth = 1.11
   * Very Smooth = 1.0
  */

  double hn;

  if (cos(tilt) == 0.0)
  {
    hn = 1.31*(pow(fabs(Tsurf - Tamb),1.0/3.0));
  }
  else if ((cos(tilt) < 0.0 && Tsurf < Tamb) ||
       (cos(tilt) > 0.0 && Tsurf > Tamb))
  {
    hn = 9.482*(pow(fabs(Tsurf - Tamb),1.0/3.0))/
        (7.238 - fabs(cos(tilt)));
  }
  else if ((cos(tilt) < 0.0 && Tsurf > Tamb) ||
       (cos(tilt) > 0.0 && Tsurf < Tamb))
  {
    hn = 1.810*(pow(fabs(Tsurf - Tamb),1.0/3.0))/
        (1.382 - fabs(cos(tilt)));
  }

  double hcGlass;

  if (isWindward(tilt, azimuth, windDirection))
  {
    hcGlass = sqrt(pow(hn,2) + pow(pow(3.26*Vair,0.89),2));
  }
  else
  {
    hcGlass = sqrt(pow(hn,2) + pow(pow(3.55*Vair,0.617),2));
  }

  double hf = roughness*(hcGlass - hn);
  return hn + hf;
}

bool isWindward(double tilt, double azimuth, double windDirection)
{
  if (fabs(cos(tilt)) < 0.98)
  {
    double diff = fabs(windDirection - azimuth);
    if ((diff - PI) > 0.001)
    {
      diff -= 2*PI;
    }
    if (fabs((diff) - 100.0*PI/180.0) > 0.001)
    {
      return false;
    }
    else return true;
  }
  else return true;
}

double getExteriorIRCoeff(double eSurf, double Tsurf, double Tamb, double eSky, double tilt)
{

  double F = getEffectiveExteriorViewFactor(eSky, tilt);
  return eSurf*SIGMA*(Tamb*Tamb*pow(F,0.5)+Tsurf*Tsurf)*
      (Tamb*pow(F,0.25)+Tsurf);
}

double getEffectiveExteriorViewFactor(double eSky, double tilt)
{
  double Fsky = 0.5*(1.0 + cos(tilt));
  double beta = cos(tilt*0.5);
  return Fsky*beta*(eSky - 1.0) + 1.0;

}

double getSimpleInteriorIRCoeff(double eSurf, double Tsurf, double Tamb)
{
  return eSurf*SIGMA*(Tamb*Tamb+Tsurf*Tsurf)*(Tamb + Tsurf);
}

#endif
