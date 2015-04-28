/* Mesher.h is part of Kiva (Written by Neal Kruis)
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

#ifndef Mesher_HPP
#define Mesher_HPP

#include <vector>
#include <cmath>
#include <math.h>
#include "Functions.h"

class Interval
{
public:

  double maxGrowthCoeff;
  double minCellDim;
  enum Growth
  {
    FORWARD,
    BACKWARD,
    UNIFORM,
    CENTERED
  };
  Growth growthDir;

};

class MeshData
{
public:
  std::vector<double> points;
  std::vector<Interval> intervals;

};

class Mesher
{
private:
    MeshData data;

public:
    std::vector<double> dividers; // center is between divider[i] and divider[i+1]
    std::vector<double> deltas;
    std::vector<double> centers;

public:

    Mesher();
    Mesher(MeshData &data);
    std::size_t getNearestIndex(double position);
    std::size_t getNextIndex(double position);
    std::size_t getPreviousIndex(double position);

};


#endif
