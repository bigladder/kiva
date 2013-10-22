/* Domain.c++ is part of Kiva (Written by Neal Kruis)
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

#ifndef Domain_CPP
#define Domain_CPP

#include "Domain.h"

Domain::Domain()
{

}

Domain::Domain(Foundation &foundation)
{
	setDomain(foundation);
}

// Having this separate from the constructor allows the correct resizing of
// multidimensional arrays on pre-existing initialized instances.
void Domain::setDomain(Foundation &foundation)
{

	{
		Mesher mX(foundation.xMeshData);
		meshX = mX;

		Mesher mZ(foundation.zMeshData);
		meshZ = mZ;
	}

	nX = meshX.centers.size();
	nY = meshY.centers.size();
	nZ = meshZ.centers.size();

	cell.resize(boost::extents[nX][nY][nZ]);

	for (size_t k = 0; k < nZ; k++)
	{
		for (size_t j = 0; j < nY; j++)
		{
			for (size_t i = 0; i < nX; i++)
			{

				// Set Cell Properties
				cell[i][j][k].density = foundation.soil.density;
				cell[i][j][k].specificHeat = foundation.soil.specificHeat;
				cell[i][j][k].conductivity = foundation.soil.conductivity;

				// Default to normal cells
				cell[i][j][k].cellType = Cell::NORMAL;

				// Next set interior zero-width cells
				if (meshX.deltas[i] < 0.00000001 ||
					meshZ.deltas[k] < 0.00000001)
				{
					cell[i][j][k].cellType = Cell::ZERO_THICKNESS;
				}

				for (size_t b = 0; b < foundation.blocks.size(); b++)
				{
					if (meshX.centers[i] > foundation.blocks[b].xMin &&
						meshX.centers[i] < foundation.blocks[b].xMax &&
						meshZ.centers[k] > foundation.blocks[b].zMin &&
						meshZ.centers[k] < foundation.blocks[b].zMax)
					{
						cell[i][j][k].density = foundation.blocks[b].material.density;
						cell[i][j][k].specificHeat = foundation.blocks[b].material.specificHeat;
						cell[i][j][k].conductivity = foundation.blocks[b].material.conductivity;

						cell[i][j][k].blockNumber = b;
						cell[i][j][k].block = foundation.blocks[b];


						if (foundation.blocks[b].blockType == Block::INTERIOR_AIR)
						{
							cell[i][j][k].cellType = Cell::INTERIOR_AIR;
						}
						else if (foundation.blocks[b].blockType == Block::EXTERIOR_AIR)
						{
							cell[i][j][k].cellType = Cell::EXTERIOR_AIR;
						}
					}
				}

				for (size_t s = 0; s < foundation.surfaces.size(); s++)
				{
					if (meshX.centers[i] - foundation.surfaces[s].xMin > -0.00000001
					&&  meshX.centers[i] - foundation.surfaces[s].xMax < 0.00000001)
					{
						if (meshY.centers[j] - foundation.surfaces[s].yMin > -0.00000001
						&&  meshY.centers[j] - foundation.surfaces[s].yMax < 0.00000001)
						{
							if (meshZ.centers[k] - foundation.surfaces[s].zMin > -0.00000001
							&&  meshZ.centers[k] - foundation.surfaces[s].zMax < 0.00000001)
							{
								cell[i][j][k].cellType = Cell::BOUNDARY;

								cell[i][j][k].surfaceNumber = s;

								cell[i][j][k].surface = foundation.surfaces[s];

								// Point/Line cells not on the boundary should be
								// zero-thickness cells
								if ((meshX.deltas[i] < 0.00000001 &&
									meshZ.deltas[k] < 0.00000001) &&
									i != 0 && i != nX -1 && k != 0 && k != nZ - 1)
								{
									cell[i][j][k].cellType = Cell::ZERO_THICKNESS;
								}
							}
						}
					}
				}
			}
		}

		// Calculate matrix coefficients
		for (size_t k = 0; k < nZ; k++)
		{
			for (size_t j = 0; j < nY; j++)
			{
				for (size_t i = 0; i < nX; i++)
				{

					// First set effective properties of zero-thickness cells
					// based on other cells
					if (i != 0 && i != nX - 1 && k != 0 && k != nZ - 1)
					{
						if (meshX.deltas[i] < 0.00000001 &&
							meshZ.deltas[k] < 0.00000001)
						{
							double volTL = meshX.deltas[i-1]*meshZ.deltas[k+1];
							double volTR = meshX.deltas[i+1]*meshZ.deltas[k+1];
							double volBL = meshX.deltas[i-1]*meshZ.deltas[k-1];
							double volBR = meshX.deltas[i+1]*meshZ.deltas[k-1];
							double volTotal = volTL + volTR + volBL + volBR;

							double densTL = cell[i-1][j][k+1].density;
							double densTR = cell[i+1][j][k+1].density;
							double densBL = cell[i-1][j][k-1].density;
							double densBR = cell[i+1][j][k-1].density;

							cell[i][j][k].density = (volTL*densTL + volTR*densTR +
												volBL*densBL + volBR*densBR) /
												volTotal;

							double cpTL = cell[i-1][j][k+1].specificHeat;
							double cpTR = cell[i+1][j][k+1].specificHeat;
							double cpBL = cell[i-1][j][k-1].specificHeat;
							double cpBR = cell[i+1][j][k-1].specificHeat;

							cell[i][j][k].specificHeat = (volTL*cpTL*densTL +
													 volTR*cpTR*densTR +
													 volBL*cpBL*densBL +
													 volBR*cpBR*densBR) /
													 (volTotal*cell[i][j][k].density);

							double condTL = cell[i-1][j][k+1].conductivity;
							double condTR = cell[i+1][j][k+1].conductivity;
							double condBL = cell[i-1][j][k-1].conductivity;
							double condBR = cell[i+1][j][k-1].conductivity;

							cell[i][j][k].conductivity = (volTL*condTL +
													 volTR*condTR +
													 volBL*condBL +
													 volBR*condBR) /
													 volTotal;
						}
						else if (meshX.deltas[i] < 0.00000001)
						{
							double volL = meshX.deltas[i-1];
							double volR = meshX.deltas[i+1];
							double volTotal = volL + volR;

							double densL = cell[i-1][j][k].density;
							double densR = cell[i+1][j][k].density;

							cell[i][j][k].density = (volL*densL + volR*densR) /
												volTotal;

							double cpL = cell[i-1][j][k].specificHeat;
							double cpR = cell[i+1][j][k].specificHeat;

							cell[i][j][k].specificHeat = (volL*cpL*densL +
													 volR*cpR*densR) /
													 (volTotal*cell[i][j][k].density);

							double condL = cell[i-1][j][k].conductivity;
							double condR = cell[i+1][j][k].conductivity;

							cell[i][j][k].conductivity = (volL*condL +
													 volR*condR) /
													 volTotal;

						}
						else if (meshZ.deltas[k] < 0.00000001)
						{
							double volT = meshZ.deltas[k+1];
							double volB = meshZ.deltas[k-1];
							double volTotal = volT + volB;

							double densT = cell[i][j][k+1].density;
							double densB = cell[i][j][k-1].density;

							cell[i][j][k].density = (volT*densT + volB*densB) /
												volTotal;

							double cpT = cell[i][j][k+1].specificHeat;
							double cpB = cell[i][j][k-1].specificHeat;

							cell[i][j][k].specificHeat = (volT*cpT*densT +
													 volB*cpB*densB) /
													 (volTotal*cell[i][j][k].density);

							double condT = cell[i][j][k+1].conductivity;
							double condB = cell[i][j][k-1].conductivity;

							cell[i][j][k].conductivity = (volT*condT +
													 volB*condB) /
													 volTotal;
						}
					}

					// PDE Coefficients

					// Radial X terms
					if (foundation.coordinateSystem == Foundation::CS_2DAXIAL)
					{
						cell[i][j][k].cxp_c = (getDXM(i)*getKXP(i,j,k))/
								((getDXM(i) + getDXP(i))*getDXP(i));
						cell[i][j][k].cxm_c = (getDXP(i)*getKXM(i,j,k))/
								((getDXM(i) + getDXP(i))*getDXM(i));
					}
					else
					{
						cell[i][j][k].cxp_c = 0.0;
						cell[i][j][k].cxm_c = 0.0;
					}

					// Cartesian X terms
					cell[i][j][k].cxp = (2*getKXP(i,j,k))/
							((getDXM(i) + getDXP(i))*getDXP(i));
					cell[i][j][k].cxm = -1*(2*getKXM(i,j,k))/
							((getDXM(i) + getDXP(i))*getDXM(i));

					// Cartesian Y terms
					cell[i][j][k].czp = (2*getKZP(i,j,k))/
							((getDZM(k) + getDZP(k))*getDZP(k));
					cell[i][j][k].czm = -1*(2*getKZM(i,j,k))/
							((getDZM(k) + getDZP(k))*getDZM(k));

					// Cartesian Y terms
					if (foundation.coordinateSystem == Foundation::CS_3D)
					{
						cell[i][j][k].cyp = (2*getKYP(i,j,k))/
								((getDYM(j) + getDYP(j))*getDYP(j));
						cell[i][j][k].cym = -1*(2*getKYM(i,j,k))/
								((getDYM(j) + getDYP(j))*getDYM(j));
					}
					else
					{
						cell[i][j][k].cyp = 0.0;
						cell[i][j][k].cym = 0.0;
					}

				}
			}
		}
	}
}

double Domain::getDXP(size_t i)
{
	if (i == nX - 1)
	{
		// For boundary cells assume that the cell on the other side of the
		// boundary is the same as the previous cell
		return (meshX.deltas[i] + meshX.deltas[i - 1])/2.0;
	}
	else
	{
		return (meshX.deltas[i] + meshX.deltas[i + 1])/2.0;
	}
}

double Domain::getDXM(size_t i)
{
	if (i == 0)
	{
		// For boundary cells assume that the cell on the other side of the
		// boundary is the same as the previous cell
		return (meshX.deltas[i] + meshX.deltas[i + 1])/2.0;
	}
	else
	{
		return (meshX.deltas[i] + meshX.deltas[i - 1])/2.0;
	}
}

double Domain::getDYP(size_t j)
{
	if (j == nY - 1)
	{
		// For boundary cells assume that the cell on the other side of the
		// boundary is the same as the previous cell
		return (meshY.deltas[j] + meshY.deltas[j - 1])/2.0;
	}
	else
	{
		return (meshY.deltas[j] + meshY.deltas[j + 1])/2.0;
	}
}

double Domain::getDYM(size_t j)
{
	if (j == 0)
	{
		// For boundary cells assume that the cell on the other side of the
		// boundary is the same as the previous cell
		return (meshY.deltas[j] + meshY.deltas[j + 1])/2.0;
	}
	else
	{
		return (meshY.deltas[j] + meshY.deltas[j - 1])/2.0;
	}
}

double Domain::getDZP(size_t k)
{
	if (k == nZ - 1)
	{
		// For boundary cells assume that the cell on the other side of the
		// boundary is the same as the previous cell
		return (meshZ.deltas[k] + meshZ.deltas[k - 1])/2.0;
	}
	else
	{
		return (meshZ.deltas[k] + meshZ.deltas[k + 1])/2.0;
	}
}

double Domain::getDZM(size_t k)
{
	if (k == 0)
	{
		// For boundary cells assume that the cell on the other side of the
		// boundary is the same as the previous cell
		return (meshZ.deltas[k] + meshZ.deltas[k + 1])/2.0;
	}
	else
	{
		return (meshZ.deltas[k] + meshZ.deltas[k - 1])/2.0;
	}
}

double Domain::getKXP(size_t i, size_t j, size_t k)
{
	if (i == nX - 1)
	{
		// For boundary cells assume that the cell on the other side of the
		// boundary is the same as the current cell
		return cell[i][j][k].conductivity;
	}
	else
	{
		return 1/(meshX.deltas[i]/(2*getDXP(i)*cell[i][j][k].conductivity) +
				meshX.deltas[i + 1]/(2*getDXP(i)*cell[i+1][j][k].conductivity));
	}
}

double Domain::getKXM(size_t i, size_t j, size_t k)
{
	if (i == 0)
	{
		// For boundary cells assume that the cell on the other side of the
		// boundary is the same as the current cell
		return cell[i][j][k].conductivity;
	}
	else
	{
		return 1/(meshX.deltas[i]/(2*getDXM(i)*cell[i][j][k].conductivity) +
				meshX.deltas[i - 1]/(2*getDXM(i)*cell[i-1][j][k].conductivity));
	}
}

double Domain::getKYP(size_t i, size_t j, size_t k)
{
	if (j == nY - 1)
	{
		// For boundary cells assume that the cell on the other side of the
		// boundary is the same as the current cell
		return cell[i][j][k].conductivity;
	}
	else
	{
		return 1/(meshY.deltas[j]/(2*getDYP(j)*cell[i][j][k].conductivity) +
				meshY.deltas[j + 1]/(2*getDYP(j)*cell[i][j+1][k].conductivity));
	}
}

double Domain::getKYM(size_t i, size_t j, size_t k)
{
	if (j == 0)
	{
		// For boundary cells assume that the cell on the other side of the
		// boundary is the same as the current cell
		return cell[i][j][k].conductivity;
	}
	else
	{
		return 1/(meshY.deltas[j]/(2*getDYM(j)*cell[i][j][k].conductivity) +
				meshY.deltas[j - 1]/(2*getDYM(j)*cell[i][j-1][k].conductivity));
	}
}

double Domain::getKZP(size_t i, size_t j, size_t k)
{
	if (k == nZ - 1)
	{
		// For boundary cells assume that the cell on the other side of the
		// boundary is the same as the current cell
		return cell[i][j][k].conductivity;
	}
	else
	{
		return 1/(meshZ.deltas[k]/(2*getDZP(k)*cell[i][j][k].conductivity) +
				meshZ.deltas[k + 1]/(2*getDZP(k)*cell[i][j][k+1].conductivity));
	}
}

double Domain::getKZM(size_t i, size_t j, size_t k)
{
	if (k == 0)
	{
		// For boundary cells assume that the cell on the other side of the
		// boundary is the same as the current cell
		return cell[i][j][k].conductivity;
	}
	else
	{
		return 1/(meshZ.deltas[k]/(2*getDZM(k)*cell[i][j][k].conductivity) +
				meshZ.deltas[k - 1]/(2*getDZM(k)*cell[i][j][k-1].conductivity));
	}
}

void Domain::printCellTypes()
{
	// TODO: Make the ability to output a specific slice in i, j, or k
	std::ofstream output;
	output.open("Cells.csv");

	for (size_t i = 0; i < nX; i++)
	{

		output << ", " << i;

	}

	output << std::endl;

	for (size_t k = nZ - 1; k >= 0 && k < nZ; k--)
	{

		output << k;

		for (size_t i = 0; i < nX; i++)
		{

			output << ", " << cell[i][0][k].cellType;

		}

		output << std::endl;
	}
	output.close();

}

#endif
