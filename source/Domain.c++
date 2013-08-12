#ifndef Domain_CPP
#define Domain_CPP

#include "Domain.h"
#include <fstream>


Domain::Domain()
{

}

Domain::Domain(Mesher &mesher, Foundation &foundation, double &timestep)  :
			   mesher(mesher)
{

	nR = mesher.xCenters.size();
	nZ = mesher.yCenters.size();

	rho = blas::matrix<double>(nR, nZ);
	cp = blas::matrix<double>(nR, nZ);
	k = blas::matrix<double>(nR, nZ);

	theta = blas::matrix<double>(nR, nZ);
	a = blas::matrix<double>(nR, nZ);
	b = blas::matrix<double>(nR, nZ);
	c = blas::matrix<double>(nR, nZ);
	d = blas::matrix<double>(nR, nZ);
	e = blas::matrix<double>(nR, nZ);
	f = blas::matrix<double>(nR, nZ);

	cellType = blas::matrix<CellType>(nR, nZ);

	// Slab and Wall indices
	bool slabJSet = false;
	bool slabIminSet = false;

	// z dimensions
	double zWallTop = foundation.wall.height - foundation.wall.depth;
	double zGrade = 0.0;
	// double zInsEdge = zWallTop - foundation.interiorVerticalInsulation.depth;
	double zSlab = zWallTop - foundation.excavationDepth;
	double zDeepGround = -foundation.deepGroundDepth;

	// r dimensions
	double rAxis = 0.0;
	double rIntIns = foundation.radius - foundation.interiorVerticalInsulation.layer.d;
	double rExtIns = foundation.radius +
					 foundation.wall.totalWidth() +
					 foundation.exteriorVerticalInsulation.layer.d;
	double rFarField = foundation.radius + foundation.farFieldWidth;


	for (size_t j = 0; j < nZ; j++)
	{
		for (size_t i = 0; i < nR; i++)
		{

			rho(i, j) = foundation.soil.rho;
			cp(i, j) = foundation.soil.cp;
			k(i, j) = foundation.soil.k;

			for (size_t b = 0; b < foundation.blocks.size(); b++)
			{
				if (mesher.xCenters[i] > foundation.blocks[b].rMin &&
					mesher.xCenters[i] < foundation.blocks[b].rMax &&
					mesher.yCenters[j] > foundation.blocks[b].zMin &&
					mesher.yCenters[j] < foundation.blocks[b].zMax)
				{
					rho(i, j) = foundation.blocks[b].material.rho;
					cp(i, j) = foundation.blocks[b].material.cp;
					k(i, j) = foundation.blocks[b].material.k;
				}
			}

			// Set Cell Types

			// Default to normal cells
			cellType(i, j) = NORMAL;

			// Next set zero-width cells
			if (mesher.xDeltas[i] < 0.00000001)
			{
				cellType(i, j) = ZERO_WIDTH_R;
			}

			if (mesher.yDeltas[j] < 0.00000001)
			{
				cellType(i, j) = ZERO_WIDTH_Z;
			}

			if (mesher.xDeltas[i] < 0.00000001 &&
				mesher.yDeltas[j] < 0.00000001)
			{
				cellType(i, j) = ZERO_WIDTH_RZ;
			}

			// Interior Air
			if (mesher.xCenters[i] - rIntIns < -0.00000001 &&
				mesher.yCenters[j] - zSlab > 0.00000001)
			{
				cellType(i, j) = INTERIOR_AIR;
			}

			// Exterior Air
			if (mesher.xCenters[i] - rExtIns > 0.00000001 &&
				mesher.yCenters[j] - zGrade > 0.00000001)
			{
				cellType(i, j) = EXTERIOR_AIR;
			}

			// Top of Wall
			if (fabs(mesher.yCenters[j] - zWallTop) < 0.00000001)
			{
				if (mesher.xCenters[i] - rIntIns > -0.00000001 &&
					mesher.xCenters[i] - rExtIns < 0.00000001)
				{
					cellType(i, j) = WALL_TOP;
				}
			}

			// Exterior Grade
			if (fabs(mesher.yCenters[j] - zGrade) < 0.00000001)
			{
				if (mesher.xCenters[i] - rExtIns > 0.00000001)
				{
					cellType(i, j) = EXTERIOR_GRADE;
				}
			}

			// Exterior Wall
			if (fabs(mesher.xCenters[i] - rExtIns) < 0.00000001)
			{
				if (mesher.yCenters[j] - zGrade > 0.00000001 &&
					mesher.yCenters[j] - zWallTop < -0.00000001)
				{
					cellType(i, j) = EXTERIOR_WALL;
				}
			}

			// Interior Slab
			if (fabs(mesher.yCenters[j] - zSlab) < 0.00000001)
			{
				if (! slabJSet)
				{
					slabJ = j;
					slabJSet = true;
				}

				if (mesher.xCenters[i] - rIntIns < -0.00000001)
				{
					if (! slabIminSet)
					{
						slabImin = i;
						slabIminSet = true;
					}

					slabImax = i;

					cellType(i, j) = INTERIOR_SLAB;
				}
			}

			// Interior Wall
			if (fabs(mesher.xCenters[i] - rIntIns) < 0.00000001)
			{
				if (mesher.yCenters[j] - zSlab > 0.00000001 &&
					mesher.yCenters[j] - zWallTop < -0.00000001)
				{
					cellType(i, j) = INTERIOR_WALL;
				}
			}

			// Axis
			if (fabs(mesher.xCenters[i] - rAxis) < 0.00000001)
			{
				if (mesher.yCenters[j] - zSlab < -0.00000001)
				{
					cellType(i, j) = AXIS;
				}
			}

			// Far Field
			if (fabs(mesher.xCenters[i] - rFarField) < 0.00000001)
			{
				if (mesher.yCenters[j] - zGrade < -0.00000001)
				{
					cellType(i, j) = FAR_FIELD;
				}
			}

			// Deep Ground
			if (fabs(mesher.yCenters[j] - zDeepGround) < 0.00000001)
			{
				cellType(i, j) = DEEP_GROUND;
			}
		}
	}

	// Calculate matrix coefficients
	for (size_t j = 0; j < nZ; j++)
	{
		for (size_t i = 0; i < nR; i++)
		{

			double drp;  // Delta r+
			double drm;  // Delta r-
			double dzp;  // Delta z+
			double dzm;  // Delta z-
			double krp;  // kr+
			double krm;  // kr-
			double kzp;  // kz+
			double kzm;  // kz-

			// Intermediate Values
			drp = getDRP(i);
			drm = getDRM(i);
			dzp = getDZP(j);
			dzm = getDZM(j);
			krp = getKRP(i,j);
			krm = getKRM(i,j);
			kzp = getKZP(i,j);
			kzm = getKZM(i,j);

			theta(i, j) = timestep/(rho(i, j)*cp(i, j));

			// PDE Coefficients
			if (drp < 0.00000001)
			{
				a(i, j) = 0.0;
				c(i, j) = 0.0;
			}
			else
			{
				a(i, j) = (drm*krp)/((drm + drp)*drp);
				c(i, j) = (2*krp)/((drm + drp)*drp);
			}

			if (drm < 0.00000001)
			{
				b(i, j) = 0.0;
				d(i, j) = 0.0;
			}
			else
			{
				b(i, j) = (drp*krm)/((drm + drp)*drm);
				d(i, j) = -1*(2*krm)/((drm + drp)*drm);
			}

			if (dzp < 0.00000001)
			{
				e(i, j) = 0.0;
			}
			else
			{
				e(i, j) = (2*kzp)/((dzm + dzp)*dzp);
			}

			if (dzm < 0.00000001)
			{
				f(i, j) = 0.0;
			}
			else
			{
				f(i, j) = -1*(2*kzm)/((dzm + dzp)*dzm);
			}
		}
	}
}

double Domain::getDRP(size_t i)
{
	if (i == nR - 1)
	{
		// For boundary cells assume that the cell on the other side of the
		// boundary is the same as the previous cell
		return (mesher.xDeltas[i] + mesher.xDeltas[i - 1])/2.0;
	}
	else
	{
		return (mesher.xDeltas[i] + mesher.xDeltas[i + 1])/2.0;
	}
}

double Domain::getDRM(size_t i)
{
	if (i == 0)
	{
		// For boundary cells assume that the cell on the other side of the
		// boundary is the same as the previous cell
		return (mesher.xDeltas[i] + mesher.xDeltas[i + 1])/2.0;
	}
	else
	{
		return (mesher.xDeltas[i] + mesher.xDeltas[i - 1])/2.0;
	}
}

double Domain::getDZP(size_t j)
{
	if (j == nZ - 1)
	{
		// For boundary cells assume that the cell on the other side of the
		// boundary is the same as the previous cell
		return (mesher.yDeltas[j] + mesher.yDeltas[j - 1])/2.0;
	}
	else
	{
		return (mesher.yDeltas[j] + mesher.yDeltas[j + 1])/2.0;
	}
}

double Domain::getDZM(size_t j)
{
	if (j == 0)
	{
		// For boundary cells assume that the cell on the other side of the
		// boundary is the same as the previous cell
		return (mesher.yDeltas[j] + mesher.yDeltas[j + 1])/2.0;
	}
	else
	{
		return (mesher.yDeltas[j] + mesher.yDeltas[j - 1])/2.0;
	}
}

double Domain::getKRP(size_t i, size_t j)
{
	if (i == nR - 1)
	{
		// For boundary cells assume that the cell on the other side of the
		// boundary is the same as the current cell
		return k(i,j);
	}
	else
	{
		return 1/(mesher.xDeltas[i]/(2*getDRP(i)*k(i,j)) +
				  mesher.xDeltas[i + 1]/(2*getDRP(i)*k(i + 1,j)));
	}
}

double Domain::getKRM(size_t i, size_t j)
{
	if (i == 0)
	{
		// For boundary cells assume that the cell on the other side of the
		// boundary is the same as the current cell
		return k(i,j);
	}
	else
	{
		return 1/(mesher.xDeltas[i]/(2*getDRM(i)*k(i,j)) +
				  mesher.xDeltas[i - 1]/(2*getDRM(i)*k(i - 1,j)));
	}
}

double Domain::getKZP(size_t i, size_t j)
{
	if (j == nZ - 1)
	{
		// For boundary cells assume that the cell on the other side of the
		// boundary is the same as the current cell
		return k(i,j);
	}
	else
	{
		return 1/(mesher.yDeltas[j]/(2*getDZP(j)*k(i,j)) +
				  mesher.yDeltas[j + 1]/(2*getDZP(j)*k(i,j + 1)));
	}
}

double Domain::getKZM(size_t i, size_t j)
{
	if (j == 0)
	{
		// For boundary cells assume that the cell on the other side of the
		// boundary is the same as the current cell
		return k(i,j);
	}
	else
	{
		return 1/(mesher.yDeltas[j]/(2*getDZM(j)*k(i,j)) +
				  mesher.yDeltas[j - 1]/(2*getDZM(j)*k(i,j - 1)));
	}
}

void Domain::printCellTypes()
{
	ofstream output;
	output.open("Cells.csv");

	for (size_t i = 0; i < nR; i++)
	{

		output << ", " << i;

	}

	output << endl;

	for (size_t j = nZ - 1; j >= 0 && j < nZ; j--)
	{

		output << j;

		for (size_t i = 0; i < nR; i++)
		{

			output << ", " << cellType(i,j);

		}

		output << endl;
	}
	output.close();

}

#endif
