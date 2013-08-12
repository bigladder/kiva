/* Algorithms.h is part of Kiva (Written by Neal Kruis)
 * Copyright (C) 2012-2013 Big Ladder Software <info@bigladdersoftware.com>
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

#ifndef ConvectionAlgorithms_HPP
#define ConvectionAlgorithms_HPP

#include <cmath>

using namespace std;

// TODO: use defaulting inputs
double getDOE2ConvectionCoeff(double tilt,
		  	  	  	  	  	  double azimuth,
		  	  	  	  	  	  double windDirection,
		  	  	  	  	  	  double Tsurf,
		  	  	  	  	  	  double Tamb,
		  	  	  	  	  	  double Vair,
		  	  	  	  	  	  double roughness);

bool isWindward(double tilt, double azimuth, double windDirection);

double getIRCoeff(double eSurf,
		          double Tsurf,
		          double Tamb,
		          double eSky,
		          double tilt);

#endif
