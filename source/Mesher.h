/* Mesher.h is part of Kiva (Written by Neal Kruis)
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

#ifndef Mesher_HPP
#define Mesher_HPP

#include <boost/numeric/ublas/vector.hpp>
namespace blas = boost::numeric::ublas;

class Interval
{
public:

	double maxGrowthCoeff;
	double minCellDim;
	enum Growth
	{
		FORWARD,
		BACKWARD,
		UNIFORM
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
		MeshData xData;
		MeshData yData;

public:
		blas::vector<double> xDividers, yDividers;
		blas::vector<double> xDeltas, yDeltas;
		blas::vector<double> xCenters, yCenters;

public:

		Mesher();
		Mesher(MeshData &xData, MeshData &yData);
		void makeMesh(MeshData data, blas::vector<double> &dividers,
				      blas::vector<double> &deltas,
				      blas::vector<double> &centers);

};


#endif
