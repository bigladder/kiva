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

	density.resize(boost::extents[nX][nY][nZ]);
	specificHeat.resize(boost::extents[nX][nY][nZ]);
	conductivity.resize(boost::extents[nX][nY][nZ]);

	// PDE Coefficients
	cxp_c.resize(boost::extents[nX][nY][nZ]);
	cxm_c.resize(boost::extents[nX][nY][nZ]);
	cxp.resize(boost::extents[nX][nY][nZ]);
	cxm.resize(boost::extents[nX][nY][nZ]);
	cyp.resize(boost::extents[nX][nY][nZ]);
	cym.resize(boost::extents[nX][nY][nZ]);
	czp.resize(boost::extents[nX][nY][nZ]);
	czm.resize(boost::extents[nX][nY][nZ]);

	cellType.resize(boost::extents[nX][nY][nZ]);

	// Slab and Wall indices
	bool slabKSet = false;
	bool slabIminSet = false;
	
	// z dimensions
	double zWallTop = foundation.wall.height - foundation.wall.depth;
	double zGrade = 0.0;
	// double zInsEdge = zWallTop - foundation.interiorVerticalInsulation.depth;
	double zSlab = zWallTop - foundation.excavationDepth;
	double zDeepGround = -foundation.deepGroundDepth;

	// r dimensions
	double xAxis = 0.0;
	double xIntIns = foundation.effectiveLength - foundation.interiorVerticalInsulation.layer.thickness;
	double xExtIns = foundation.effectiveLength +
					 foundation.wall.totalWidth() +
					 foundation.exteriorVerticalInsulation.layer.thickness;
	double xFarField = foundation.effectiveLength + foundation.farFieldWidth;

	for (size_t k = 0; k < nZ; k++)
	{
		for (size_t j = 0; j < nY; j++)
		{
			for (size_t i = 0; i < nX; i++)
			{

				// Set Cell Properties
				density[i][j][k] = foundation.soil.density;
				specificHeat[i][j][k] = foundation.soil.specificHeat;
				conductivity[i][j][k] = foundation.soil.conductivity;

				for (size_t b = 0; b < foundation.blocks.size(); b++)
				{
					if (meshX.centers[i] > foundation.blocks[b].xMin &&
						meshX.centers[i] < foundation.blocks[b].xMax &&
						meshZ.centers[k] > foundation.blocks[b].zMin &&
						meshZ.centers[k] < foundation.blocks[b].zMax)
					{
						density[i][j][k] = foundation.blocks[b].material.density;
						specificHeat[i][j][k] = foundation.blocks[b].material.specificHeat;
						conductivity[i][j][k] = foundation.blocks[b].material.conductivity;
					}
				}

				// Set Cell Types

				// Default to normal cells
				cellType[i][j][k] = NORMAL;

				// Next set interior zero-width cells
				if (meshX.deltas[i] < 0.00000001 ||
					meshZ.deltas[k] < 0.00000001)
				{
					cellType[i][j][k] = ZERO_THICKNESS;
				}
				// Interior Air
				if (meshX.centers[i] - xIntIns < -0.00000001 &&
					meshZ.centers[k] - zSlab > 0.00000001)
				{
					cellType[i][j][k] = INTERIOR_AIR;
				}

				// Exterior Air
				if (meshX.centers[i] - xExtIns > 0.00000001 &&
					meshZ.centers[k] - zGrade > 0.00000001)
				{
					cellType[i][j][k] = EXTERIOR_AIR;
				}

				// Top of Wall
				if (fabs(meshZ.centers[k] - zWallTop) < 0.00000001)
				{
					if (meshX.centers[i] - xIntIns > -0.00000001 &&
						meshX.centers[i] - xExtIns < 0.00000001)
					{
						cellType[i][j][k] = WALL_TOP;
					}
				}

				// Exterior Grade
				if (fabs(meshZ.centers[k] - zGrade) < 0.00000001)
				{
					if (meshX.centers[i] - xExtIns > 0.00000001)
					{
						cellType[i][j][k] = EXTERIOR_GRADE;
					}
				}

				// Exterior Wall
				if (fabs(meshX.centers[i] - xExtIns) < 0.00000001)
				{
					if (meshZ.centers[k] - zGrade > 0.00000001 &&
						meshZ.centers[k] - zWallTop < -0.00000001)
					{
						cellType[i][j][k] = EXTERIOR_WALL;
					}
				}

				// Interior Slab
				if (fabs(meshZ.centers[k] - zSlab) < 0.00000001)
				{
					if (! slabKSet)
					{
						slabK = k;
						slabKSet = true;
					}

					if (meshX.centers[i] - xIntIns < -0.00000001)
					{
						if (! slabIminSet)
						{
							slabImin = i;
							slabIminSet = true;
						}

						slabImax = i;

						cellType[i][j][k] = INTERIOR_SLAB;
					}
				}

				// Interior Wall
				if (fabs(meshX.centers[i] - xIntIns) < 0.00000001)
				{
					if (meshZ.centers[k] - zSlab > 0.00000001 &&
						meshZ.centers[k] - zWallTop < -0.00000001)
					{
						cellType[i][j][k] = INTERIOR_WALL;
					}
				}

				// Axis
				if (fabs(meshX.centers[i] - xAxis) < 0.00000001)
				{
					if (meshZ.centers[k] - zSlab < -0.00000001)
					{
						cellType[i][j][k] = SYMMETRY;
					}
				}

				// Far Field
				if (fabs(meshX.centers[i] - xFarField) < 0.00000001)
				{
					if (meshZ.centers[k] - zGrade < -0.00000001)
					{
						cellType[i][j][k] = FAR_FIELD;
					}
				}

				// Deep Ground
				if (fabs(meshZ.centers[k] - zDeepGround) < 0.00000001)
				{
					cellType[i][j][k] = DEEP_GROUND;
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

							double densTL = density[i-1][j][k+1];
							double densTR = density[i+1][j][k+1];
							double densBL = density[i-1][j][k-1];
							double densBR = density[i+1][j][k-1];

							density[i][j][k] = (volTL*densTL + volTR*densTR +
												volBL*densBL + volBR*densBR) /
												volTotal;

							double cpTL = specificHeat[i-1][j][k+1];
							double cpTR = specificHeat[i+1][j][k+1];
							double cpBL = specificHeat[i-1][j][k-1];
							double cpBR = specificHeat[i+1][j][k-1];

							specificHeat[i][j][k] = (volTL*cpTL*densTL +
													 volTR*cpTR*densTR +
													 volBL*cpBL*densBL +
													 volBR*cpBR*densBR) /
													 (volTotal*density[i][j][k]);

							double condTL = conductivity[i-1][j][k+1];
							double condTR = conductivity[i+1][j][k+1];
							double condBL = conductivity[i-1][j][k-1];
							double condBR = conductivity[i+1][j][k-1];

							conductivity[i][j][k] = (volTL*condTL +
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

							double densL = density[i-1][j][k];
							double densR = density[i+1][j][k];

							density[i][j][k] = (volL*densL + volR*densR) /
												volTotal;

							double cpL = specificHeat[i-1][j][k];
							double cpR = specificHeat[i+1][j][k];

							specificHeat[i][j][k] = (volL*cpL*densL +
													 volR*cpR*densR) /
													 (volTotal*density[i][j][k]);

							double condL = conductivity[i-1][j][k];
							double condR = conductivity[i+1][j][k];

							conductivity[i][j][k] = (volL*condL +
													 volR*condR) /
													 volTotal;

						}
						else if (meshZ.deltas[k] < 0.00000001)
						{
							double volT = meshZ.deltas[k+1];
							double volB = meshZ.deltas[k-1];
							double volTotal = volT + volB;

							double densT = density[i][j][k+1];
							double densB = density[i][j][k-1];

							density[i][j][k] = (volT*densT + volB*densB) /
												volTotal;

							double cpT = specificHeat[i][j][k+1];
							double cpB = specificHeat[i][j][k-1];

							specificHeat[i][j][k] = (volT*cpT*densT +
													 volB*cpB*densB) /
													 (volTotal*density[i][j][k]);

							double condT = conductivity[i][j][k+1];
							double condB = conductivity[i][j][k-1];

							conductivity[i][j][k] = (volT*condT +
													 volB*condB) /
													 volTotal;
						}
					}

					double dxp;  // Delta x+
					double dxm;  // Delta x-
					double dyp;  // Delta y+
					double dym;  // Delta y-
					double dzp;  // Delta z+
					double dzm;  // Delta z-
					double kxp;  // kx+
					double kxm;  // kx-
					double kyp;  // ky+
					double kym;  // ky-
					double kzp;  // kz+
					double kzm;  // kz-

					// Intermediate Values
					dxp = getDXP(i);
					dxm = getDXM(i);
					dyp = getDYP(j);
					dym = getDYM(j);
					dzp = getDZP(k);
					dzm = getDZM(k);
					kxp = getKXP(i,j,k);
					kxm = getKXM(i,j,k);
					kyp = getKYP(i,j,k);
					kym = getKYM(i,j,k);
					kzp = getKZP(i,j,k);
					kzm = getKZM(i,j,k);

					// PDE Coefficients
					if (dxp < 0.00000001)
					{
						cxp_c[i][j][k] = 0.0;
						cxp[i][j][k] = 0.0;
					}
					else
					{
						cxp_c[i][j][k] = (dxm*kxp)/((dxm + dxp)*dxp);
						cxp[i][j][k] = (2*kxp)/((dxm + dxp)*dxp);
					}

					if (dxm < 0.00000001)
					{
						cxm_c[i][j][k] = 0.0;
						cxm[i][j][k] = 0.0;
					}
					else
					{
						cxm_c[i][j][k] = (dxp*kxm)/((dxm + dxp)*dxm);
						cxm[i][j][k] = -1*(2*kxm)/((dxm + dxp)*dxm);
					}

					if (dyp < 0.00000001)
					{
						cyp[i][j][k] = 0.0;
					}
					else
					{
						cyp[i][j][k] = (2*kyp)/((dym + dyp)*dyp);
					}

					if (dym < 0.00000001)
					{
						cym[i][j][k] = 0.0;
					}
					else
					{
						cym[i][j][k] = -1*(2*kym)/((dym + dyp)*dym);
					}

					if (dzp < 0.00000001)
					{
						czp[i][j][k] = 0.0;
					}
					else
					{
						czp[i][j][k] = (2*kzp)/((dzm + dzp)*dzp);
					}

					if (dzm < 0.00000001)
					{
						czm[i][j][k] = 0.0;
					}
					else
					{
						czm[i][j][k] = -1*(2*kzm)/((dzm + dzp)*dzm);
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
		return conductivity[i][j][k];
	}
	else
	{
		return 1/(meshX.deltas[i]/(2*getDXP(i)*conductivity[i][j][k]) +
				meshX.deltas[i + 1]/(2*getDXP(i)*conductivity[i+1][j][k]));
	}
}

double Domain::getKXM(size_t i, size_t j, size_t k)
{
	if (i == 0)
	{
		// For boundary cells assume that the cell on the other side of the
		// boundary is the same as the current cell
		return conductivity[i][j][k];
	}
	else
	{
		return 1/(meshX.deltas[i]/(2*getDXM(i)*conductivity[i][j][k]) +
				meshX.deltas[i - 1]/(2*getDXM(i)*conductivity[i-1][j][k]));
	}
}

double Domain::getKYP(size_t i, size_t j, size_t k)
{
	if (j == nY - 1)
	{
		// For boundary cells assume that the cell on the other side of the
		// boundary is the same as the current cell
		return conductivity[i][j][k];
	}
	else
	{
		return 1/(meshY.deltas[j]/(2*getDYP(j)*conductivity[i][j][k]) +
				meshY.deltas[j + 1]/(2*getDYP(j)*conductivity[i][j+1][k]));
	}
}

double Domain::getKYM(size_t i, size_t j, size_t k)
{
	if (j == 0)
	{
		// For boundary cells assume that the cell on the other side of the
		// boundary is the same as the current cell
		return conductivity[i][j][k];
	}
	else
	{
		return 1/(meshY.deltas[j]/(2*getDYM(j)*conductivity[i][j][k]) +
				meshY.deltas[j - 1]/(2*getDYM(j)*conductivity[i][j-1][k]));
	}
}

double Domain::getKZP(size_t i, size_t j, size_t k)
{
	if (k == nZ - 1)
	{
		// For boundary cells assume that the cell on the other side of the
		// boundary is the same as the current cell
		return conductivity[i][j][k];
	}
	else
	{
		return 1/(meshZ.deltas[k]/(2*getDZP(k)*conductivity[i][j][k]) +
				meshZ.deltas[k + 1]/(2*getDZP(k)*conductivity[i][j][k+1]));
	}
}

double Domain::getKZM(size_t i, size_t j, size_t k)
{
	if (k == 0)
	{
		// For boundary cells assume that the cell on the other side of the
		// boundary is the same as the current cell
		return conductivity[i][j][k];
	}
	else
	{
		return 1/(meshZ.deltas[k]/(2*getDZM(k)*conductivity[i][j][k]) +
				meshZ.deltas[k - 1]/(2*getDZM(k)*conductivity[i][j][k-1]));
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

			output << ", " << cellType[i][0][k];

		}

		output << std::endl;
	}
	output.close();

}

#endif
