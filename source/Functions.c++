/* Functions.c++ is part of Kiva (Written by Neal Kruis)
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

#ifndef Functions_CPP
#define Functions_CPP

#include "Functions.h"

static const double EPSILON = 1E-7;

bool isLessThan(double first, double second)
{
	if(first - second < -EPSILON)
		return true;
	else
		return false;
}

bool isLessOrEqual(double first, double second)
{
	if(first - second < EPSILON)
		return true;
	else
		return false;
}

bool isEqual(double first, double second)
{
	if(fabs(first - second) < EPSILON)
		return true;
	else
		return false;
}

bool isGreaterThan(double first, double second)
{
	if(first - second > EPSILON)
		return true;
	else
		return false;
}

bool isGreaterOrEqual(double first, double second)
{
	if(first - second > -EPSILON)
		return true;
	else
		return false;
}

bool isOdd(int N)
{
	return (N % 2 != 0);
}

bool isEven(int N)
{
	return (N % 2 == 0);
}

#endif



