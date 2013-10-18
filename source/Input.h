/* Input.h is part of Kiva (Written by Neal Kruis)
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

#ifndef INPUT_HPP_
#define INPUT_HPP_

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include "Mesher.h"

class SimulationControl
{
public:
	// Simulation Control
	boost::gregorian::date startDate;
	boost::gregorian::date endDate;
	boost::posix_time::time_duration timestep;
	std::string weatherFile;
	boost::posix_time::ptime startTime;

	void setStartTime()
	{
		boost::posix_time::ptime st(startDate,boost::posix_time::hours(0));
		startTime = st;
	}
};

class Material
{
public:

	double conductivity;  // [W/m-K] conductivity (boost function of z, t?)
	double density;  // [kg/m3] density
	double specificHeat;  // [J/kg-K] specific heat
};

class Layer
{
public:

	Material material;
	double thickness;  // [m] thickness
};

class HorizontalInsulation
{
public:

	double depth;  // [m] depth from top of wall
	double width;  // [m] width from side of wall
	Layer layer;


};

class VerticalInsulation
{
public:

	double depth; // [m] depth from top of wall
	Layer layer;

};

class Wall
{
public:

	double interiorEmissivity;
	double exteriorEmissivity;
	double exteriorAbsorptivity;
	double depth;  // [m] below grade depth
	double height;  // [m] total height
	std::vector <Layer> layers;

	double totalWidth()
	{
		double width = 0.0;

		for (size_t n = 0; n < layers.size(); n++) width += layers[n].thickness;

		return width;
	}
};

class Slab
{
public:

	double emissivity;
	std::vector <Layer> layers;

	double totalWidth()
	{
		double width = 0.0;

		for (size_t n = 0; n < layers.size(); n++) width += layers[n].thickness;

		return width;
	}

};

class Mesh
{
public:

	double maxRadialGrowthCoeff;
	double maxInteriorGrowthCoeff;
	double maxDepthGrowthCoeff;
	double minCellDim;  // [m]

};

class OutputAnimation
{
public:

	std::string name;
	boost::posix_time::time_duration frequency;
	bool grid;
	bool contours;
	bool gradients;
	int size;
};

class Block
{
public:

	double xMin, xMax, yMin, yMax, zMin, zMax;
	Material material;

	Block()
	{

	}

	Block(Material mat, double xmin, double xmax, double zmin, double zmax)
	{
		material = mat;
		xMin = xmin;
		xMax = xmax;
		yMin = 0.5;
		yMax = 0.5;
		zMin = zmin;
		zMax = zmax;
	}
};

class Surface
{
public:
	double xMin, xMax, yMin, yMax, zMin, zMax;
	std::string name;

	Surface()
	{

	}

	Surface(std::string name_, double xmin, double xmax, double zmin, double zmax)
	{
		name = name_;
		xMin = xmin;
		xMax = xmax;
		yMin = 0.5;
		yMax = 0.5;
		zMin = zmin;
		zMax = zmax;
	}

	bool constX()
	{
		if(xMin == xMax)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	bool constY()
	{
		if(yMin == yMax)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	bool constZ()
	{
		if(zMin == zMax)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

};

class Foundation
{
public:

	// Inputs

	// Site
	double deepGroundDepth;  // [m]
	double farFieldWidth;  // [m] distance from outside of wall to the edge
						   // of the domain
	double deepGroundTemperature;  // [K]
	double excavationDepth; // [m] below top of wall

	enum DeepGroundBoundary
	{
		DGB_AUTO,
		DGB_CONSTANT_TEMPERATURE,
		DGB_ZERO_FLUX
	};

	DeepGroundBoundary deepGroundBoundary;

	double indoorAirTemperature; // [K]

	Material soil;
	double soilAbsorptivity;  // [frac]
	double soilEmissivity;  // [frac]

	// Local wind speed characteristics
	double vegetationHeight;  // [m]
	double deltaLocal;  // [m]
	double alphaLocal;  // [-]


	// Geometry
	enum CoordinateSystem
	{
		CS_2DAXIAL,
		CS_2DLINEAR,
		CS_3D
	};
	CoordinateSystem coordinateSystem;
	double length;  // [m]
	double width;  // [m]

	// Constructions
	Wall wall;
	bool hasWall;
	Slab slab;
	bool hasSlab;
	HorizontalInsulation interiorHorizontalInsulation;
	bool hasInteriorHorizontalInsulation;
	HorizontalInsulation exteriorHorizontalInsulation;
	bool hasExteriorHorizontalInsulation;
	VerticalInsulation interiorVerticalInsulation;
	bool hasInteriorVerticalInsulation;
	VerticalInsulation exteriorVerticalInsulation;
	bool hasExteriorVerticalInsulation;


	// Meshing
	Mesh mesh;


	// Simulation Control
	enum NumericalScheme
	{
		NS_ADE,
		NS_EXPLICIT,
		NS_IMPLICIT,
		NS_CRANK_NICOLSON,
		NS_STEADY_STATE
	};

	NumericalScheme numericalScheme;

	enum GeometryTransform
	{
		GT_CYLINDRICAL,
		GT_CARTESIAN
	};

	GeometryTransform geometryTransform;

	double initialTemperature;
	long implicitAccelTimestep;
	long implicitAccelPeriods;
	enum InitializationMethod
	{
		IM_KUSUDA,
		IM_CONSTANT_TEMPERATURE,
		IM_IMPLICIT_ACCELERATION,
		IM_STEADY_STATE
	};

	InitializationMethod initializationMethod;

	double interiorConvectiveCoefficient;
	double exteriorConvectiveCoefficient;
	enum ConvectionCalculationMethod
	{
		CCM_AUTO,
		CCM_CONSTANT_COEFFICIENT
	};

	ConvectionCalculationMethod convectionCalculationMethod;

	double outdoorDryBulbTemperature;
	enum OutdoorTemperatureMethod
	{
		OTM_WEATHER_FILE,
		OTM_CONSTANT_TEMPERATURE
	};

	OutdoorTemperatureMethod outdoorTemperatureMethod;


	// Output Animations
	OutputAnimation outputAnimation;


	// Output Reports


	// Derived variables
	double area;  // [m2] Area of foundation
	double perimeter;  // [m] Perimeter of foundation
	double effectiveLength;

	MeshData xMeshData;
	MeshData zMeshData;
	std::vector<Block> blocks;
	std::vector<Surface> surfaces;

	void setMeshData()
	{
		area = length*width;  // [m2] Area of foundation
		perimeter = 2*length + 2*width;  // [m] Perimeter of foundation
		effectiveLength = 2.0*area/perimeter;

		// points are the coordinate values of the interval boundaries.
		// Intervals explain how each interval between a set of points is
		// discretized.

		Material air;
		air.conductivity = 0.02587;
		air.density = 1.275;
		air.specificHeat = 1007;

		// Meshing info
		Interval zeroThickness;
		zeroThickness.maxGrowthCoeff = 1.0;
		zeroThickness.minCellDim = 1.0;
		zeroThickness.growthDir = Interval::UNIFORM;

		Interval near;
		near.maxGrowthCoeff = 1.0;
		near.minCellDim = mesh.minCellDim;
		near.growthDir = Interval::UNIFORM;

		Interval deep;
		deep.maxGrowthCoeff = mesh.maxDepthGrowthCoeff;
		deep.minCellDim = mesh.minCellDim;
		deep.growthDir = Interval::BACKWARD;

		Interval interior;
		interior.maxGrowthCoeff = mesh.maxInteriorGrowthCoeff;
		interior.minCellDim = mesh.minCellDim;
		interior.growthDir = Interval::BACKWARD;

		Interval exterior;
		exterior.maxGrowthCoeff = mesh.maxRadialGrowthCoeff;
		exterior.minCellDim = mesh.minCellDim;
		exterior.growthDir = Interval::FORWARD;


		if (coordinateSystem == CS_2DAXIAL ||
			coordinateSystem == CS_2DLINEAR)
		{
			// Set misc. dimensions
			double zTop;
			if (hasWall)
				zTop = wall.height - wall.depth;
			else
				zTop = 0.0;

			double zBottom = -deepGroundDepth;

			double xLeft = 0.0;
			double xRight = effectiveLength + farFieldWidth;

			double zGrade = 0.0;
			double zSlab = zTop - excavationDepth;

			double zNearBottom = std::min(zSlab, zGrade);  // This will change depending on configuration

    		double xWallExterior;
    		if (hasWall)
    		{
        		if (hasExteriorVerticalInsulation)
        		{
        			xWallExterior = effectiveLength
            				+ wall.totalWidth()
            				+ exteriorVerticalInsulation.layer.thickness;

        		}
        		else
        		{
        			xWallExterior = effectiveLength
            				+ wall.totalWidth();
        		}
    		}
    		else
    		{
    			xWallExterior = effectiveLength;
    		}

    		double xWallInterior;
    		if (hasInteriorVerticalInsulation)
    		{
    			xWallInterior = effectiveLength
    					- interiorVerticalInsulation.layer.thickness;

    		}
    		else
    		{
    			xWallInterior = effectiveLength;
    		}

    		double xNearLeft = effectiveLength;  // This will change depending on configuration

    		double xNearRight = effectiveLength;  // This will change depending on configuration

			// Deep ground surface
			Surface deepGround("Deep Ground",
					xLeft,
					xRight,
					zBottom,
					zBottom);
			surfaces.push_back(deepGround);

			// Slab
			Surface slabInterior("Slab Interior",
					xLeft,
					xWallInterior,
					zSlab,
					zSlab);
			surfaces.push_back(slabInterior);

			// Grade
			Surface grade("Grade",
					xWallExterior,
					xRight,
					zGrade,
					zGrade);
			surfaces.push_back(grade);

			// Wall Top
			if(hasWall)
			{
				Surface wallTop("Wall Top",
						xWallInterior,
						xWallExterior,
						zTop,
						zTop);
				surfaces.push_back(wallTop);
			}

			// Symmetry Surface
			Surface symmetry("Symmetry",
					xLeft,
					xLeft,
					zBottom,
					zSlab);
			surfaces.push_back(symmetry);

			if(excavationDepth > 0.0)
			{
				// Interior Wall Surface
				Surface interiorWall("Interior Wall",
						xWallInterior,
						xWallInterior,
						zSlab,
						zTop);
				surfaces.push_back(interiorWall);
			}

			if(zTop > 0.0)
			{
				// Exterior Wall Surface
				Surface exteriorWall("Exterior Wall",
						xWallExterior,
						xWallExterior,
						zGrade,
						zTop);
				surfaces.push_back(exteriorWall);
			}

			// Far Field Surface
			Surface farField("Far Field",
					xRight,
					xRight,
					zBottom,
					zGrade);
			surfaces.push_back(farField);

			// Interior Horizontal Insulation
			if (hasInteriorHorizontalInsulation)
			{
				Block block(interiorHorizontalInsulation.layer.material,
						    xNearLeft,
						    effectiveLength,
						    zTop - interiorHorizontalInsulation.depth - interiorHorizontalInsulation.layer.thickness,
						    zTop - interiorHorizontalInsulation.depth);
				blocks.push_back(block);

				if(block.zMin < zNearBottom)
				{
					zNearBottom = block.zMin;
				}
				if(block.xMin < xNearLeft)
				{
					xNearLeft = block.xMin;
				}
			}

			if (hasSlab)
			{
				// Foundation Slab
				double zPosition = zSlab - slab.totalWidth();

				for (size_t n = 0; n < slab.layers.size(); n++)
				{
					Block block(slab.layers[n].material,
							    xLeft,
							    effectiveLength,
							    zPosition,
							    zPosition + slab.layers[n].thickness);
					blocks.push_back(block);
					zPosition = block.zMax;

					if(block.zMin < zNearBottom)
					{
						zNearBottom = block.zMin;
					}
				}

			}

			// Interior Vertical Insulation
			if (hasInteriorVerticalInsulation)
			{
				Block block(interiorVerticalInsulation.layer.material,
				  			xWallInterior,
				  			effectiveLength,
				  			zTop - interiorVerticalInsulation.depth,
				  			zTop);
				blocks.push_back(block);

				if(block.zMin < zNearBottom)
				{
					zNearBottom = block.zMin;
				}
				if(block.xMin < xNearLeft)
				{
					xNearLeft = block.xMin;
				}

			}

			// Indoor Air
			{
				Block block(air,
					  		xLeft,
					  		xWallInterior,
					  		zTop,
					  		zSlab);
				blocks.push_back(block);
			}

			if (hasWall)
			{
				double xPosition = effectiveLength;

				// Foundation Wall
				for (int n = wall.layers.size() - 1; n >= 0; n--)
				{
					Block block(wall.layers[n].material,
								xPosition,
								xPosition + wall.layers[n].thickness,
								-wall.depth,
						  		zTop);
					xPosition = block.xMax;
					blocks.push_back(block);

					if(block.zMin < zNearBottom)
					{
						zNearBottom = block.zMin;
					}
					if(block.xMax > xNearRight)
					{
						xNearRight = block.xMax;
					}

				}
			}

			// Exterior Vertical Insulation
			if (hasExteriorVerticalInsulation)
			{
				Block block(exteriorVerticalInsulation.layer.material,
						    xWallExterior - exteriorVerticalInsulation.layer.thickness,
					  		xWallExterior,
					  		zTop - exteriorVerticalInsulation.depth,
					  		zTop);
				blocks.push_back(block);

				if(block.zMin < zNearBottom)
				{
					zNearBottom = block.zMin;
				}
				if(block.xMax > xNearRight)
				{
					xNearRight = block.xMax;
				}

			}

			// Exterior Horizontal Insulation
			if (hasExteriorHorizontalInsulation)
			{
				Block block(exteriorHorizontalInsulation.layer.material,
						    effectiveLength + wall.totalWidth(),
						    xWallExterior,
					  		zTop - exteriorHorizontalInsulation.depth - exteriorHorizontalInsulation.layer.thickness,
					  		zTop - exteriorHorizontalInsulation.depth);
				blocks.push_back(block);

				if(block.zMin < zNearBottom)
				{
					zNearBottom = block.zMin;
				}
				if(block.xMax > xNearRight)
				{
					xNearRight = block.xMax;
				}

			}

			// Exterior Air
			{
				Block block(air,
							xWallExterior,
					  		xRight,
					  		zGrade,
					  		zTop);
				blocks.push_back(block);
			}

			std::vector<double> xPoints;
			std::vector<double> zPoints;

			std::vector<double> xSurfaces;
			std::vector<double> zSurfaces;

			for (size_t s = 0; s < surfaces.size(); s++)
			{
				xPoints.push_back(surfaces[s].xMin);
				xPoints.push_back(surfaces[s].xMax);
				zPoints.push_back(surfaces[s].zMax);
				zPoints.push_back(surfaces[s].zMin);

				if (surfaces[s].constX())
				{
					xSurfaces.push_back(surfaces[s].xMin);
				}
				if (surfaces[s].constZ())
				{
					zSurfaces.push_back(surfaces[s].zMin);
				}

			}

			for (size_t b = 0; b < blocks.size(); b++)
			{
				xPoints.push_back(blocks[b].xMin);
				xPoints.push_back(blocks[b].xMax);
				zPoints.push_back(blocks[b].zMax);
				zPoints.push_back(blocks[b].zMin);
			}

			// Sort the points vectors
			sort(xPoints.begin(), xPoints.end());
			sort(zPoints.begin(), zPoints.end());

			// erase (approximately) duplicate elements
			for (size_t i = 1; i < xPoints.size(); i++)
			{
				if (std::fabs(xPoints[i] - xPoints[i-1]) < 0.0000001)
				{
					xPoints.erase(xPoints.begin() + i-1);
					i -= 1;
				}
			}

			for (size_t k = 1; k < zPoints.size(); k++)
			{
				if (std::fabs(zPoints[k] - zPoints[k-1]) < 0.0000001)
				{
					zPoints.erase(zPoints.begin() + k-1);
					k -= 1;
				}
			}

			// Sort the surfaces vectors
			sort(xSurfaces.begin(), xSurfaces.end());
			sort(zSurfaces.begin(), zSurfaces.end());

			// erase (approximately) duplicate elements
			for (size_t i = 1; i < xSurfaces.size(); i++)
			{
				if (std::fabs(xSurfaces[i] - xSurfaces[i-1]) < 0.0000001)
				{
					xSurfaces.erase(xSurfaces.begin() + i-1);
					i -= 1;
				}
			}

			for (size_t k = 1; k < zSurfaces.size(); k++)
			{
				if (std::fabs(zSurfaces[k] - zSurfaces[k-1]) < 0.0000001)
				{
					zSurfaces.erase(zSurfaces.begin() + k-1);
					k -= 1;
				}
			}

			// re-add the extra surface elements to create zero-thickness cells
			for (size_t i = 0; i < xSurfaces.size(); i++)
			{
				xPoints.push_back(xSurfaces[i]);
			}

			for (size_t k = 0; k < zSurfaces.size(); k++)
			{
				zPoints.push_back(zSurfaces[k]);
			}

			// Sort the range vectors again this time it includes doubles for zero thickness cells
			sort(xPoints.begin(), xPoints.end());
			sort(zPoints.begin(), zPoints.end());


			xMeshData.points = xPoints;
			zMeshData.points = zPoints;

			std::vector<Interval> xIntervals;
			std::vector<Interval> zIntervals;

			for (size_t i = 1; i < xPoints.size(); i++)
			{
				if (std::fabs(xPoints[i] - xPoints[i-1]) < 0.0000001)
					xIntervals.push_back(zeroThickness);
				else if (xPoints[i] - xNearLeft < 0.0000001)
					xIntervals.push_back(interior);
				else if (xPoints[i] - xNearRight < 0.0000001)
					xIntervals.push_back(near);
				else if (xPoints[i] > xNearRight)
					xIntervals.push_back(exterior);
			}

			for (size_t k = 1; k < zPoints.size(); k++)
			{
				if (std::fabs(zPoints[k] - zPoints[k-1]) < 0.0000001)
					zIntervals.push_back(zeroThickness);
				else if (zPoints[k] - zNearBottom < 0.0000001)
					zIntervals.push_back(deep);
				else if (zPoints[k] > zNearBottom)
					zIntervals.push_back(near);
			}

			xMeshData.intervals = xIntervals;
			zMeshData.intervals = zIntervals;

		}
	}

};


class Input
{
public:
	SimulationControl simulationControl;
	std::vector <Foundation> foundations;
};

#endif /* INPUT_HPP_ */
