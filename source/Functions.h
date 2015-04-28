/* Functions.h is part of Kiva (Written by Neal Kruis)
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

#ifndef Functions_HPP
#define Functions_HPP

#include <math.h>
#include <algorithm>
#include <vector>

bool isLessThan(double first, double second);
bool isLessOrEqual(double first, double second);
bool isEqual(double first, double second);
bool isEqual(double first, double second, double epsilon);
bool isGreaterThan(double first, double second);
bool isGreaterOrEqual(double first, double second);
bool isEven(int N);
bool isOdd(int N);
void solveTDM(const std::vector<double>& a1, const std::vector<double>& a2,
          std::vector<double>& a3, std::vector<double>& b,
              std::vector<double>& x);

#endif
