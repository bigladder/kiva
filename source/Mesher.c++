/* Mesher.c++ is part of Kiva (Written by Neal Kruis)
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

#ifndef Mesher_CPP
#define Mesher_CPP

#include "Mesher.h"

Mesher::Mesher()
{
  dividers.push_back(0.0);
  dividers.push_back(1.0);
  deltas.push_back(1.0);
  centers.push_back(0.5);
}

Mesher::Mesher(MeshData &data) : data(data)
{
	dividers.push_back(data.points[0]);

	// Loop through intervals
	for (std::size_t i = 0; i < data.points.size() -1; ++i)
	{
		double min = data.points[i];
		double max = data.points[i+1];
		double length = max - min;
		double cellWidth;
		int numCells;

		if (isEqual(length, 0.0))
		{
			// Zero width cells (used for boundary conditions)
			dividers.push_back(max);
			deltas.push_back(0.0);
			centers.push_back(max);
		}
		else
		{
			// Uniform Meshing
			if (data.intervals[i].growthDir == Interval::UNIFORM)
			{
				numCells = length / data.intervals[i].minCellDim;

				// Make sure that there is at least one cell
				if (numCells == 0) numCells = 1;

				if (isEqual(length/(numCells + 1),
				data.intervals[i].minCellDim))
				{
					cellWidth = data.intervals[i].minCellDim;
					numCells = numCells + 1;
				}
				else
				{
					cellWidth = length / numCells;
				}

				for (int j = 1; j <= numCells; ++j)
				{
					dividers.push_back(min + j*cellWidth);
					deltas.push_back(cellWidth);
					centers.push_back(min + j*cellWidth - cellWidth/2.0);
				}

			}
			else
			{
				// Geometric Meshing
				std::vector<double> temp;

				// Create temp. delta array
				if (isGreaterOrEqual(data.intervals[i].minCellDim, length))
				{
					// If min. cell width is greater than interval length then
					// set width equal to the length
					temp.push_back(length);

				}
				else
				{
					bool search = true;
					int N = 0;
					double multiplier;

					while (search)
					{
						multiplier = 0.0;
						for (int j = 0; j <= N; j++)
						{
							multiplier +=
								pow(data.intervals[i].maxGrowthCoeff,j);
						}

						if (data.intervals[i].minCellDim*multiplier > length)
						{
							numCells = N;
							multiplier -=
								pow(data.intervals[i].maxGrowthCoeff,N);
							search = false;
						}
						else
						{
							N += 1;
						}
					}
					double initialCellWidth = length / multiplier;
					temp.push_back(initialCellWidth);

					for (int j = 1; j < numCells; j++)
					{
						temp.push_back(temp[j - 1]
						               *data.intervals[i].maxGrowthCoeff);
					}
				}

				// build arrays
				if (data.intervals[i].growthDir == Interval::FORWARD)
				{
					double position = min;
					for (int j = 0; j < numCells; ++j)
					{
						dividers.push_back(position + temp[j]);
						deltas.push_back(temp[j]);
						centers.push_back(position + temp[j]/2.0);
						position += temp[j];
					}

				}
				else
				{
					double position = min;
					for (int j = 1; j <= numCells; ++j)
					{
						dividers.push_back(position +
								           temp[numCells - j]);
						deltas.push_back(temp[numCells - j]);
						centers.push_back(position +
								            temp[numCells - j]/2.0);
						position += temp[numCells - j];
					}
				}
			}
		}
	}
}

std::size_t Mesher::getNearestIndex(double position)
{
	if (isLessOrEqual(position,this->centers[0]))
		return 0;
	else if (isGreaterOrEqual(position,this->centers[this->centers.size() - 1]))
		return this->centers.size() - 1;
	else
	{
		for (std::size_t i = 1; i < this->centers.size(); i++)
		{
			if (isGreaterOrEqual(position,this->centers[i-1]) &&
				isLessOrEqual(position,this->centers[i]))
			{
				double diffDown = position - this->centers[i-1];
				double diffUp = this->centers[i] - position;

				if (isLessOrEqual(diffDown, diffUp))
					return i-1;
				else
					return i;
			}
		}
	}
}

std::size_t Mesher::getNextIndex(double position)
{
	if (isLessThan(position,this->centers[0]))
		return 0;
	else if (isGreaterOrEqual(position,this->centers[this->centers.size() - 1]))
		return this->centers.size() - 1;
	else
	{
		for (std::size_t i = 1; i < this->centers.size(); i++)
		{
			if (isGreaterOrEqual(position,this->centers[i-1]) &&
				isLessThan(position,this->centers[i]))
				return i;
		}
	}
}

std::size_t Mesher::getPreviousIndex(double position)
{
	if (isLessOrEqual(position,this->centers[0]))
		return 0;
	else if (isGreaterThan(position,this->centers[this->centers.size() - 1]))
		return this->centers.size() - 1;
	else
	{
		for (std::size_t i = 1; i < this->centers.size(); i++)
		{
			if (isGreaterThan(position,this->centers[i-1]) &&
				isLessOrEqual(position,this->centers[i]))
				return i-1;
		}
	}
}

#endif
