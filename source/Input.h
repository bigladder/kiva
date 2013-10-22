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
#include "Functions.h"

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

	enum BlockType
	{
		SOLID,
		INTERIOR_AIR,
		EXTERIOR_AIR
	};

	BlockType blockType;

	Block()
	{

	}

	Block(Material mat,
		  BlockType bt,
		  double xmin, double xmax, double zmin, double zmax)
	{
		material = mat;
		blockType = bt;
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
	double emissivity, absorptivity, temperature;


	enum BoundaryConditionType
	{
		ZERO_FLUX,
		INTERIOR_FLUX,
		EXTERIOR_FLUX,
		CONSTANT_TEMPERATURE,
		INTERIOR_TEMPERATURE,
		EXTERIOR_TEMPERATURE
	};
	BoundaryConditionType boundaryConditionType;

	enum Orientation
	{
		X_POS,
		X_NEG,
		Y_POS,
		Y_NEG,
		Z_POS,
		Z_NEG
	};
	Orientation orientation;

	Surface()
	{

	}

	Surface(std::string name_,
			double xmin, double xmax, double zmin, double zmax,
			BoundaryConditionType bc, Orientation o,
			double emiss, double abs)
	{
		name = name_;
		xMin = xmin;
		xMax = xmax;
		yMin = 0.5;
		yMax = 0.5;
		zMin = zmin;
		zMax = zmax;
		boundaryConditionType = bc;
		orientation = o;
		emissivity = emiss;
		absorptivity = abs;
	}

	bool isConstX()
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
	bool isConstY()
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
	bool isConstZ()
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

		// Set misc. "Z" dimensions (relative to grade)
		double zMax;
		if (hasWall)
			zMax = wall.height - wall.depth;
		else
			zMax = 0.0;

		double zMin = -deepGroundDepth;

		double zGrade = 0.0;
		double zSlab = zMax - excavationDepth;

		double zNearDeep = std::min(zSlab, zGrade);  // This will change depending on configuration

		// Set misc. "X/Y" dimensions (relative to foundation outline --
		// currently the interior of the wall)

		double xyWallExterior;
		if (hasWall)
		{
    		if (hasExteriorVerticalInsulation)
    		{
    			xyWallExterior = wall.totalWidth()
        				+ exteriorVerticalInsulation.layer.thickness;

    		}
    		else
    		{
    			xyWallExterior = wall.totalWidth();
    		}
		}
		else
		{
			xyWallExterior = 0.0;
		}

		double xyWallInterior;
		if (hasInteriorVerticalInsulation)
		{
			xyWallInterior = -interiorVerticalInsulation.layer.thickness;

		}
		else
		{
			xyWallInterior = 0.0;
		}

		double xyNearInt = xyWallInterior;  // This will change depending on configuration
		double xyNearExt = xyWallExterior;  // This will change depending on configuration

		double xyIntHIns = 0.0;
		double zIntHIns = 0.0;
		if (hasInteriorHorizontalInsulation)
		{
			xyIntHIns = -interiorHorizontalInsulation.width;
			zIntHIns = zMax - interiorHorizontalInsulation.depth - interiorHorizontalInsulation.layer.thickness;

			if(zIntHIns < zNearDeep)
			{
				zNearDeep = zIntHIns;
			}
			if(xyIntHIns < xyNearInt)
			{
				xyNearInt = xyIntHIns;
			}

		}

		double zSlabBottom = zSlab;
		if (hasSlab)
		{
			zSlabBottom = zSlab - slab.totalWidth();

			if(zSlabBottom < zNearDeep)
			{
				zNearDeep = zSlabBottom;
			}
		}

		double zIntVIns = 0.0;
		if (hasInteriorVerticalInsulation)
		{
			zIntVIns = zMax -interiorVerticalInsulation.depth;

			if(zIntVIns < zNearDeep)
			{
				zNearDeep = zIntVIns;
			}
		}

		double zWall = 0.0;
		if (hasWall)
		{
			zWall = -wall.depth;
			if(zWall < zNearDeep)
			{
				zNearDeep = zWall;
			}
		}

		double zExtVIns = 0.0;
		if (hasExteriorVerticalInsulation)
		{
			zExtVIns = zMax - exteriorVerticalInsulation.depth;

			if(zExtVIns < zNearDeep)
			{
				zNearDeep = zExtVIns;
			}
		}

		double xyExtHIns = 0.0;
		double zExtHIns = 0.0;
		if (hasExteriorHorizontalInsulation)
		{
			xyExtHIns = exteriorHorizontalInsulation.width;
			zExtHIns = zMax - exteriorHorizontalInsulation.depth - exteriorHorizontalInsulation.layer.thickness;

			if(zExtHIns < zNearDeep)
			{
				zNearDeep = zExtHIns;
			}
			if(xyExtHIns > xyNearExt)
			{
				xyNearExt = xyExtHIns;
			}


		}

		if (coordinateSystem == CS_2DAXIAL ||
			coordinateSystem == CS_2DLINEAR)
		{

			double xMin = 0.0;
			double xMax = effectiveLength + farFieldWidth;

			double xRef = effectiveLength;

			// Symmetry Surface
			Surface symmetry("Symmetry",
					xMin,
					xMin,
					zMin,
					zSlab,
					Surface::ZERO_FLUX,
					Surface::X_NEG,
					0.8,
					0.8);
			surfaces.push_back(symmetry);

			if(excavationDepth > 0.0)
			{
				// Interior Wall Surface
				Surface interiorWall("Interior Wall",
						xRef + xyWallInterior,
						xRef + xyWallInterior,
						zSlab,
						zMax,
						Surface::INTERIOR_FLUX,
						Surface::X_NEG,
						wall.interiorEmissivity,
						0.8);
				surfaces.push_back(interiorWall);

				// Interior Air Left Temperature
				Surface interiorAirLeft("Interior Air Left",
						xMin,
						xMin,
						zSlab,
						zMax,
						Surface::INTERIOR_TEMPERATURE,
						Surface::X_NEG,
						0.8,
						0.8);
				surfaces.push_back(interiorAirLeft);

			}

			if(zMax > 0.0)
			{
				// Exterior Wall Surface
				Surface exteriorWall("Exterior Wall",
						xRef + xyWallExterior,
						xRef + xyWallExterior,
						zGrade,
						zMax,
						Surface::EXTERIOR_FLUX,
						Surface::X_POS,
						wall.exteriorEmissivity,
						wall.exteriorAbsorptivity);
				surfaces.push_back(exteriorWall);

				// Exterior Air Right Surface
				Surface exteriorAirRight("Exterior Air Right",
						xMax,
						xMax,
						zGrade,
						zMax,
						Surface::EXTERIOR_TEMPERATURE,
						Surface::X_POS,
						0.8,
						0.8);
				surfaces.push_back(exteriorAirRight);
			}

			// Far Field Surface
			Surface farField("Far Field",
					xMax,
					xMax,
					zMin,
					zGrade,
					Surface::ZERO_FLUX,
					Surface::X_POS,
					0.8,
					0.8);
			surfaces.push_back(farField);

			// Deep ground surface
    		if (deepGroundBoundary == DGB_CONSTANT_TEMPERATURE ||
    			deepGroundBoundary == DGB_AUTO)
    		{
    			Surface deepGround("Deep Ground",
    					xMin,
    					xMax,
    					zMin,
    					zMin,
    					Surface::CONSTANT_TEMPERATURE,
    					Surface::Z_NEG,
    					0.8,
    					0.8);
   				deepGround.temperature = deepGroundTemperature;
    			surfaces.push_back(deepGround);

    		}
    		else if (deepGroundBoundary == DGB_ZERO_FLUX)
    		{
    			Surface deepGround("Deep Ground",
    					xMin,
    					xMax,
    					zMin,
    					zMin,
    					Surface::ZERO_FLUX,
    					Surface::Z_NEG,
    					0.8,
    					0.8);
    			surfaces.push_back(deepGround);
    		}

			// Slab
			Surface slabInterior("Slab Interior",
					xMin,
					xRef + xyWallInterior,
					zSlab,
					zSlab,
					Surface::INTERIOR_FLUX,
					Surface::Z_POS,
					wall.interiorEmissivity,
					0.8);
			surfaces.push_back(slabInterior);

			// Grade
			Surface grade("Grade",
					xRef + xyWallExterior,
					xMax,
					zGrade,
					zGrade,
					Surface::EXTERIOR_FLUX,
					Surface::Z_POS,
					soilEmissivity,
					soilAbsorptivity);
			surfaces.push_back(grade);

			if(excavationDepth > 0.0)
			{
				// Interior Air Top Surface
				Surface interiorAirTop("Interior Air Top",
						xMin,
						xRef + xyWallInterior,
						zMax,
						zMax,
						Surface::INTERIOR_TEMPERATURE,
						Surface::Z_POS,
						0.8,
						0.8);
				surfaces.push_back(interiorAirTop);
			}

			if(zMax > 0.0)
			{
				// Exterior Air Top Surface
				Surface exteriorAirTop("Exterior Air Top",
						xRef + xyWallExterior,
						xMax,
						zMax,
						zMax,
						Surface::EXTERIOR_TEMPERATURE,
						Surface::Z_POS,
						0.8,
						0.8);
				surfaces.push_back(exteriorAirTop);
			}

			// Wall Top
			if(hasWall)
			{
				Surface wallTop("Wall Top",
						xRef + xyWallInterior,
						xRef + xyWallExterior,
						zMax,
						zMax,
						Surface::ZERO_FLUX,
						Surface::Z_POS,
    					0.8,
    					0.8);
				surfaces.push_back(wallTop);
			}

			// Interior Horizontal Insulation
			if (hasInteriorHorizontalInsulation)
			{
				Block block(interiorHorizontalInsulation.layer.material,
						    Block::SOLID,
						    xRef + xyIntHIns,
						    effectiveLength,
						    zIntHIns,
						    zIntHIns + interiorHorizontalInsulation.layer.thickness);
				blocks.push_back(block);
			}

			if (hasSlab)
			{
				// Foundation Slab
				double zPosition = zSlabBottom;

				for (size_t n = 0; n < slab.layers.size(); n++)
				{
					Block block(slab.layers[n].material,
								Block::SOLID,
								xMin,
							    xRef,
							    zPosition,
							    zPosition + slab.layers[n].thickness);
					blocks.push_back(block);
					zPosition = block.zMax;
				}

			}

			// Interior Vertical Insulation
			if (hasInteriorVerticalInsulation)
			{
				Block block(interiorVerticalInsulation.layer.material,
							Block::SOLID,
							xRef + xyWallInterior,
				  			xRef,
				  			zIntVIns,
				  			zMax);
				blocks.push_back(block);
			}

			// Indoor Air
			{
				Block block(air,
						    Block::INTERIOR_AIR,
						    xMin,
						    xRef + xyWallInterior,
					  		zSlab,
					  		zMax);
				blocks.push_back(block);
			}

			if (hasWall)
			{
				double xPosition = xRef;

				// Foundation Wall
				for (int n = wall.layers.size() - 1; n >= 0; n--)
				{
					Block block(wall.layers[n].material,
								Block::SOLID,
								xPosition,
								xPosition + wall.layers[n].thickness,
								zWall,
								zMax);
					xPosition = block.xMax;
					blocks.push_back(block);

				}
			}

			// Exterior Vertical Insulation
			if (hasExteriorVerticalInsulation)
			{
				Block block(exteriorVerticalInsulation.layer.material,
							Block::SOLID,
							xRef + xyWallExterior - exteriorVerticalInsulation.layer.thickness,
							xRef + xyWallExterior,
							zExtVIns,
					  		zMax);
				blocks.push_back(block);
			}

			// Exterior Horizontal Insulation
			if (hasExteriorHorizontalInsulation)
			{
				Block block(exteriorHorizontalInsulation.layer.material,
							Block::SOLID,
							xRef + wall.totalWidth(),
							xRef + xyWallExterior,
							zExtHIns,
							zExtHIns + exteriorHorizontalInsulation.layer.thickness);
				blocks.push_back(block);
			}

			// Exterior Air
			{
				Block block(air,
							Block::EXTERIOR_AIR,
							xRef + xyWallExterior,
							xMax,
					  		zGrade,
					  		zMax);
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

				if (surfaces[s].isConstX())
				{
					xSurfaces.push_back(surfaces[s].xMin);
				}
				if (surfaces[s].isConstZ())
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

			// erase duplicate elements
			for (size_t i = 1; i < xPoints.size(); i++)
			{
				if (isEqual(xPoints[i], xPoints[i-1]))
				{
					xPoints.erase(xPoints.begin() + i-1);
					i -= 1;
				}
			}

			for (size_t k = 1; k < zPoints.size(); k++)
			{
				if (isEqual(zPoints[k], zPoints[k-1]))
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
				if (isEqual(xSurfaces[i], xSurfaces[i-1]))
				{
					xSurfaces.erase(xSurfaces.begin() + i-1);
					i -= 1;
				}
			}

			for (size_t k = 1; k < zSurfaces.size(); k++)
			{
				if (isEqual(zSurfaces[k], zSurfaces[k-1]))
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
				if (isEqual(xPoints[i], xPoints[i-1]))
					xIntervals.push_back(zeroThickness);
				else if (isLessOrEqual(xPoints[i], xRef + xyNearInt))
					xIntervals.push_back(interior);
				else if (isLessOrEqual(xPoints[i], xRef + xyNearExt))
					xIntervals.push_back(near);
				else if (isGreaterThan(xPoints[i], xRef + xyNearExt))
					xIntervals.push_back(exterior);
			}

			for (size_t k = 1; k < zPoints.size(); k++)
			{
				if (isEqual(zPoints[k], zPoints[k-1]))
					zIntervals.push_back(zeroThickness);
				else if (isLessOrEqual(zPoints[k], zNearDeep))
					zIntervals.push_back(deep);
				else if (isGreaterThan(zPoints[k], zNearDeep))
					zIntervals.push_back(near);
			}

			xMeshData.intervals = xIntervals;
			zMeshData.intervals = zIntervals;

		}
		else  // if(coordinateSystem == CS_3D)
		{

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
