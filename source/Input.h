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

	double rMin, rMax, zMin, zMax;
	Material material;

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
	double area;  // [m2] Area of foundation
	double perimeter;  // [m] Perimeter of foundation

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
	double radius;
	MeshData rMeshData;
	MeshData zMeshData;
	std::vector<Block> blocks;

	void setMeshData()
	{
		radius = 2.0*area/perimeter;

		// points are the coordinate values of the interval boundaries.
		// Intervals explain how each interval between a set of points is
		// discretized.

		// Meshing info
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

		std::vector<double> rRange;
		std::vector<double> zRange;

		double rPosition;
		double zPosition;

		Material null;
		null.conductivity = 1;
		null.density = 1;
		null.specificHeat = 1;

		// Cylinder Axis
		rRange.push_back(0.0);

		// Deep ground
		zRange.push_back(-deepGroundDepth);

		// Grade
		zRange.push_back(0.0);

		// Foundation Radius
		rRange.push_back(radius);
		rPosition = radius;

		// Top of domain
		double zTop;
		if (hasWall)
			zTop = wall.height - wall.depth;
		else
			zTop = 0.0;

		zRange.push_back(zTop);

		// Interior Horizontal Insulation
		if (hasInteriorHorizontalInsulation)
		{
			Block block;
			block.rMin = rPosition - interiorHorizontalInsulation.width;
			block.rMax = rPosition;
			block.zMax = zTop - interiorHorizontalInsulation.depth;
			block.zMin = block.zMax - interiorHorizontalInsulation.layer.thickness;
			block.material = interiorHorizontalInsulation.layer.material;
			rRange.push_back(block.rMin);
			rRange.push_back(block.rMax);
			zRange.push_back(block.zMax);
			zRange.push_back(block.zMin);
			blocks.push_back(block);
		}

		double zSlab;
		if (hasSlab)
		{
			// Foundation Slab
			zSlab = zTop - excavationDepth;
			zPosition = zSlab - slab.totalWidth();

			for (size_t n = 0; n < slab.layers.size(); n++)
			{
				Block block;
				block.rMin = 0.0;
				block.rMax = rPosition;
				block.zMin = zPosition;
				block.zMax = block.zMin + slab.layers[n].thickness;
				zPosition = block.zMax;
				block.material = slab.layers[n].material;
				rRange.push_back(block.rMin);
				rRange.push_back(block.rMax);
				zRange.push_back(block.zMax);
				zRange.push_back(block.zMin);
				blocks.push_back(block);

			}

		}
		else
		{
			zSlab = zTop - excavationDepth;
			zPosition = zSlab;
			zRange.push_back(zPosition);
		}

		// Interior Vertical Insulation
		if (hasInteriorVerticalInsulation)
		{
			Block block;
			block.rMin = rPosition - interiorVerticalInsulation.layer.thickness;
			block.rMax = rPosition;
			block.zMax = zTop;
			block.zMin = block.zMax - interiorVerticalInsulation.depth;
			block.material = interiorVerticalInsulation.layer.material;
			rRange.push_back(block.rMin);
			rRange.push_back(block.rMax);
			zRange.push_back(block.zMax);
			zRange.push_back(block.zMin);
			blocks.push_back(block);
		}

		// Indoor Air
		{
			Block block;
			block.rMin = 0.0;
			block.rMax = radius;
			block.zMax = zTop;
			block.zMin = zSlab;
			block.material = null;
			rRange.push_back(block.rMin);
			rRange.push_back(block.rMax);
			zRange.push_back(block.zMax);
			zRange.push_back(block.zMin);
			blocks.push_back(block);
		}

		if (hasWall)
		{
			// Foundation Wall
			for (int n = wall.layers.size() - 1; n >= 0; n--)
			{
				Block block;
				block.rMin = rPosition;
				block.rMax = block.rMin + wall.layers[n].thickness;
				rPosition = block.rMax;
				block.zMax = zTop;
				block.zMin = -1*wall.depth;
				block.material = wall.layers[n].material;
				rRange.push_back(block.rMin);
				rRange.push_back(block.rMax);
				zRange.push_back(block.zMax);
				zRange.push_back(block.zMin);
				blocks.push_back(block);
			}
		}
		else if (excavationDepth > 0.0)
			rRange.push_back(radius);

		// Exterior Vertical Insulation
		if (hasExteriorVerticalInsulation)
		{
			Block block;
			block.rMin = rPosition;
			block.rMax = rPosition + exteriorVerticalInsulation.layer.thickness;
			block.zMax = zTop;
			block.zMin = block.zMax - exteriorVerticalInsulation.depth;
			block.material = exteriorVerticalInsulation.layer.material;
			rRange.push_back(block.rMin);
			rRange.push_back(block.rMax);
			zRange.push_back(block.zMax);
			zRange.push_back(block.zMin);
			blocks.push_back(block);
		}

		// Interior Horizontal Insulation
		if (hasExteriorHorizontalInsulation)
		{
			Block block;
			block.rMin = rPosition;
			block.rMax = rPosition + exteriorHorizontalInsulation.width;
			block.zMax = zTop - exteriorHorizontalInsulation.depth;
			block.zMin = block.zMax - exteriorHorizontalInsulation.layer.thickness;
			block.material = exteriorHorizontalInsulation.layer.material;
			rRange.push_back(block.rMin);
			rRange.push_back(block.rMax);
			zRange.push_back(block.zMax);
			zRange.push_back(block.zMin);
			blocks.push_back(block);
		}

		double rWallExt = rPosition;
		double rFarField = radius + farFieldWidth;

		// Exterior Air
		{
			Block block;
			block.rMin = rWallExt;
			block.rMax = rFarField;
			block.zMax = zTop;
			block.zMin = 0.0;
			block.material = null;
			rRange.push_back(block.rMin);
			rRange.push_back(block.rMax);
			zRange.push_back(block.zMax);
			zRange.push_back(block.zMin);
			blocks.push_back(block);
		}

		// Far field
		rRange.push_back(rFarField);

		zRange.push_back(0.0);

		// Sort the range vectors
		sort(rRange.begin(), rRange.end());
		sort(zRange.begin(), zRange.end());

		// erase (approximately) duplicate elements
		for (size_t i = 1; i < rRange.size(); i++)
		{
			if (std::fabs(rRange[i] - rRange[i-1]) < 0.0000001)
			{
				rRange.erase(rRange.begin() + i-1);
				i -= 1;
			}
		}

		for (size_t j = 1; j < zRange.size(); j++)
		{
			if (std::fabs(zRange[j] - zRange[j-1]) < 0.0000001)
			{
				zRange.erase(zRange.begin() + j-1);
				j -= 1;
			}
		}

		// insert zero-thickness cells for boundary conditions

		// Deep ground
		zRange.push_back(-deepGroundDepth);

		// Top of domain
		if (zTop > 0.0)
			zRange.push_back(zTop);

		// Grade
		zRange.push_back(0.0);

		// Axis
		rRange.push_back(0.0);

		// Wall
		if (excavationDepth > 0.0)
		{
			rRange.push_back(radius);

			if (hasWall)
				rRange.push_back(rWallExt);

			// Slab
			if (fabs(zSlab) > 0.000001)
				zRange.push_back(zSlab);
		}

		// Far field
		rRange.push_back(radius + farFieldWidth);

		// Sort the range vectors again this time it includes doubles for zero thickness cells
		sort(rRange.begin(), rRange.end());
		sort(zRange.begin(), zRange.end());


		rMeshData.points = rRange;
		zMeshData.points = zRange;

		std::vector<Interval> rintervals;
		rintervals.push_back(interior);
		rintervals.push_back(interior);

		for (size_t i = 2; i < rRange.size() - 3; i++)
		{
			rintervals.push_back(near);
		}
		rintervals.push_back(exterior);

		// Zero-thickness far-field boundary
		rintervals.push_back(exterior);

		std::vector<Interval> zintervals;

		// Zero-thickness lower boundary
		zintervals.push_back(deep);

		// Main "deep" interval
		zintervals.push_back(deep);

		for (size_t j = 2; j < zRange.size() - 1; j++)
		{
			zintervals.push_back(near);
		}

		rMeshData.intervals = rintervals;
		zMeshData.intervals = zintervals;

	}

};


class Input
{
public:
	SimulationControl simulationControl;
	std::vector <Foundation> foundations;
};

#endif /* INPUT_HPP_ */
