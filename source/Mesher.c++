// Mesher.cpp
//
// A simple mesher on a 1d domain. We divide
// an interval into J+1 mesh points, J-1 of which
// are internal mesh points.
//
// DD 2005-10-1 Mesh can be in any time interval [t1, T]
// DD 2005-12-1 hpp and cpp versions
// DD 2005-12-1 mesh in T direction
// DD 2010-2-19 Vector<double, long> Mesher::xarr(int J, int start), now start value possible
//
// (C) Datasim Education BV 2006-2010
//

#ifndef Mesher_CPP
#define Mesher_CPP

#include "Mesher.h"
#include <cmath>
#include <math.h>

Mesher::Mesher()
{

}

Mesher::Mesher(MeshData &xData, MeshData &yData) : xData(xData), yData(yData)
{ // Describe the domain of integration

	makeMesh(xData, xDividers, xDeltas, xCenters);
	makeMesh(yData, yDividers, yDeltas, yCenters);

}

void Mesher::makeMesh(MeshData data, blas::vector<double> &dividers,
	      	  	  	  blas::vector<double> &deltas,
	      	  	  	  blas::vector<double> &centers)
{
	std::vector<double> divArray;
	std::vector<double> deltArray;
	std::vector<double> centArray;
	divArray.push_back(data.points[0]);

	// Loop through intervals
	for (std::size_t i = 0; i < data.points.size() -1; ++i)
	{
		double min = data.points[i];
		double max = data.points[i+1];
		double length = max - min;
		double cellWidth;
		int numCells;

		if (length < 0.00000001)
		{
			// Zero width cells (used for boundary conditions)
			divArray.push_back(max);
			deltArray.push_back(0.0);
			centArray.push_back(max);
		}
		else
		{
			// Uniform Meshing
			if (data.intervals[i].growthDir == Interval::UNIFORM)
			{
				numCells = length / data.intervals[i].minCellDim;

				// Make sure that there is at least one cell
				if (numCells == 0) numCells = 1;

				// Check to make sure that an exact fit isn't excluded due to
				// numerical precision
				if (std::abs(length/(numCells + 1) -
				data.intervals[i].minCellDim) < 0.00000001)
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
					divArray.push_back(min + j*cellWidth);
					deltArray.push_back(cellWidth);
					centArray.push_back(min + j*cellWidth - cellWidth/2.0);
				}

			}
			else
			{
				// Geometric Meshing
				std::vector<double> temp;

				// Create temp. delta array
				if (data.intervals[i].minCellDim >= length)
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
						divArray.push_back(position + temp[j]);
						deltArray.push_back(temp[j]);
						centArray.push_back(position + temp[j]/2.0);
						position += temp[j];
					}

				}
				else
				{
					double position = min;
					for (int j = 1; j <= numCells; ++j)
					{
						divArray.push_back(position +
								           temp[numCells - j]);
						deltArray.push_back(temp[numCells - j]);
						centArray.push_back(position +
								            temp[numCells - j]/2.0);
						position += temp[numCells - j];
					}
				}
			}
		}

	}

	// Create BLAS vectors from standard vectors
	blas::vector<double> divTemp(divArray.size());
	copy(divArray.begin(), divArray.end(), divTemp.begin());
	dividers = divTemp;

	blas::vector<double> deltTemp(deltArray.size());
	copy(deltArray.begin(), deltArray.end(), deltTemp.begin());
	deltas = deltTemp;

	blas::vector<double> centTemp(centArray.size());
	copy(centArray.begin(), centArray.end(), centTemp.begin());
	centers = centTemp;


}

#endif
