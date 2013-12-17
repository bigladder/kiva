/* Input.c++ is part of Kiva (Written by Neal Kruis)
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

#ifndef INPUT_CPP_
#define INPUT_CPP_

#include "Input.h"

void SimulationControl::setStartTime()
{
	boost::posix_time::ptime st(startDate,boost::posix_time::hours(0));
	startTime = st;
}



double Wall::totalWidth()
{
	double width = 0.0;

	for (size_t n = 0; n < layers.size(); n++) width += layers[n].thickness;

	return width;
}

double Wall::totalResistance()
{
	double R = 0.0;

	for (size_t n = 0; n < layers.size(); n++) R += (layers[n].thickness/layers[n].material.conductivity);

	return R;
}


double Slab::totalWidth()
{
	double width = 0.0;

	for (size_t n = 0; n < layers.size(); n++) width += layers[n].thickness;

	return width;
}

double Slab::totalResistance()
{
	double R = 0.0;

	for (size_t n = 0; n < layers.size(); n++) R += (layers[n].thickness/layers[n].material.conductivity);

	return R;
}

OutputVariable::OutputVariable(int varID)
{
	headers.resize(14);

	headers[0] = "Slab Core Average Heat Flux [W/m2]";
	headers[1] = "Slab Core Average Temperature [K]";
	headers[2] = "Slab Core Average Effective Temperature [C]";
	headers[3] = "Slab Core Total Heat Transfer Rate [W]";
	headers[4] = "Slab Perimeter Average Heat Flux [W/m2]";
	headers[5] = "Slab Perimeter Average Temperature [K]";
	headers[6] = "Slab Perimeter Average Effective Temperature [C]";
	headers[7] = "Slab Perimeter Total Heat Transfer Rate [W]";
	headers[8] = "Slab Total Heat Transfer Rate [W]";
	headers[9] = "Wall Average Heat Flux [W/m2]";
	headers[10] = "Wall Average Temperature [K]";
	headers[11] = "Wall Average Effective Temperature [C]";
	headers[12] = "Wall Total Heat Transfer Rate [W]";
	headers[13] = "Total Heat Transfer Rate [W]";

	variableID = varID;
	headerText = headers[varID];
}

void Block::setSquarePolygon()
{
	polygon.outer().push_back(Point(xMin,yMin));
	polygon.outer().push_back(Point(xMin,yMax));
	polygon.outer().push_back(Point(xMax,yMax));
	polygon.outer().push_back(Point(xMax,yMin));
}

void Surface::setSquarePolygon()
{
	polygon.outer().push_back(Point(xMin,yMin));
	polygon.outer().push_back(Point(xMin,yMax));
	polygon.outer().push_back(Point(xMax,yMax));
	polygon.outer().push_back(Point(xMax,yMin));
}

inline bool compareRanges(RangeType first,	RangeType second)
{
	return (first.range.first < second.range.first);
}

bool Ranges::isType(double position,RangeType::Type type)
{
	// find specific Range
	for (std::size_t r = 0; r < ranges.size(); r++)
	{
		if (isGreaterThan(position,ranges[r].range.first) &&
			isLessOrEqual(position,ranges[r].range.second))
		{
			if (ranges[r].type == type)
				return true;
		}
	}
	return false;
}


void Foundation::setMeshData()
{
	area = boost::geometry::area(polygon);  // [m2] Area of foundation
	perimeter = boost::geometry::perimeter(polygon);  // [m] Perimeter of foundation
	effectiveLength = 2.0*area/perimeter;

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

	Interval minInterior;
	minInterior.maxGrowthCoeff = mesh.maxInteriorGrowthCoeff;
	minInterior.minCellDim = mesh.minCellDim;
	minInterior.growthDir = Interval::BACKWARD;

	Interval midInterior;
	midInterior.maxGrowthCoeff = mesh.maxInteriorGrowthCoeff;
	midInterior.minCellDim = mesh.minCellDim;
	midInterior.growthDir = Interval::BACKWARD;

	Interval minExterior;
	minExterior.maxGrowthCoeff = mesh.maxExteriorGrowthCoeff;
	minExterior.minCellDim = mesh.minCellDim;
	minExterior.growthDir = Interval::BACKWARD;

	Interval maxExterior;
	maxExterior.maxGrowthCoeff = mesh.maxExteriorGrowthCoeff;
	maxExterior.minCellDim = mesh.minCellDim;
	maxExterior.growthDir = Interval::FORWARD;

	// Set misc. "Z" dimensions (relative to grade)
	double zMax;
	if (hasWall)
		zMax = wall.heightAboveGrade;
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

	double xyPerimeterSurface;
	if (hasPerimeterSurface)
	{
		xyPerimeterSurface = -perimeterSurfaceWidth;
		if(xyPerimeterSurface < xyNearInt)
		{
			xyNearInt = xyPerimeterSurface;
		}
	}
	else
		xyPerimeterSurface = 0.0;


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
		zWall = wall.heightAboveGrade - wall.height;
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


	Ranges xRanges;
	Ranges yRanges;
	Ranges zRanges;

	RangeType zDeepRange;
	zDeepRange.range.first = zMin;
	zDeepRange.range.second = zNearDeep;
	zDeepRange.type = RangeType::DEEP;

	zRanges.ranges.push_back(zDeepRange);

	RangeType zNearRange;
	zNearRange.range.first = zNearDeep;
	zNearRange.range.second = zMax;
	zNearRange.type = RangeType::NEAR;

	zRanges.ranges.push_back(zNearRange);

	if (coordinateSystem == CS_2DAXIAL ||
		coordinateSystem == CS_2DLINEAR)
	{

		// TODO: 2D
		double xMin = 0.0;
		double xMax = effectiveLength + farFieldWidth;

		double xRef = effectiveLength;

		// Symmetry Surface
		{
			Surface surface;
			surface.name = "Symmetry";
			surface.xMin = xMin;
			surface.xMax = xMin;
			surface.yMin = 0.0;
			surface.yMax = 1.0;
			surface.setSquarePolygon();
			surface.zMin = zMin;
			surface.zMax = zSlab;
			surface.boundaryConditionType = Surface::ZERO_FLUX;
			surface.orientation = Surface::X_NEG;
			surfaces.push_back(surface);
		}

		if(excavationDepth > 0.0)
		{
			// Interior Wall Surface
			{
				Surface surface;
				surface.name = "Interior Wall";
				surface.xMin = xRef + xyWallInterior;
				surface.xMax = xRef + xyWallInterior;
				surface.yMin = 0.0;
				surface.yMax = 1.0;
				surface.setSquarePolygon();
				surface.zMin = zSlab;
				surface.zMax = zMax;
				surface.boundaryConditionType = Surface::INTERIOR_FLUX;
				surface.orientation = Surface::X_NEG;
				surface.emissivity = wall.interiorEmissivity;
				surfaces.push_back(surface);
			}

			// Interior Air Left Temperature
			{
			Surface surface;
			surface.name = "Interior Air Left";
			surface.xMin = xMin;
			surface.xMax = xMin;
			surface.yMin = 0.0;
			surface.yMax = 1.0;
			surface.setSquarePolygon();
			surface.zMin = zSlab;
			surface.zMax = zMax;
			surface.boundaryConditionType = Surface::INTERIOR_TEMPERATURE;
			surface.orientation = Surface::X_NEG;
			surfaces.push_back(surface);
			}
		}

		if(zMax > 0.0)
		{
			// Exterior Wall Surface
			{
				Surface surface;
				surface.name = "Exterior Wall";
				surface.xMin = xRef + xyWallExterior;
				surface.xMax = xRef + xyWallExterior;
				surface.yMin = 0.0;
				surface.yMax = 1.0;
				surface.setSquarePolygon();
				surface.zMin = zGrade;
				surface.zMax = zMax;
				surface.boundaryConditionType = Surface::EXTERIOR_FLUX;
				surface.orientation = Surface::X_POS;
				surface.emissivity = wall.exteriorEmissivity;
				surface.absorptivity = wall.exteriorAbsorptivity;
				surfaces.push_back(surface);
			}

			// Exterior Air Right Surface
			{
				Surface surface;
				surface.name = "Exterior Air Right";
				surface.xMin = xMax;
				surface.xMax = xMax;
				surface.yMin = 0.0;
				surface.yMax = 1.0;
				surface.setSquarePolygon();
				surface.zMin = zGrade;
				surface.zMax = zMax;
				surface.boundaryConditionType = Surface::EXTERIOR_TEMPERATURE;
				surface.orientation = Surface::X_POS;
				surfaces.push_back(surface);
			}
		}

		// Far Field Surface
		{
			Surface surface;
			surface.name = "Far Field";
			surface.xMin = xMax;
			surface.xMax = xMax;
			surface.yMin = 0.0;
			surface.yMax = 1.0;
			surface.setSquarePolygon();
			surface.zMin = zMin;
			surface.zMax = zGrade;
			surface.boundaryConditionType = Surface::ZERO_FLUX;
			surface.orientation = Surface::X_POS;
			surfaces.push_back(surface);
		}

		// Deep ground surface
		if (deepGroundBoundary == DGB_CONSTANT_TEMPERATURE ||
			deepGroundBoundary == DGB_AUTO)
		{
			Surface surface;
			surface.name = "Deep Ground";
			surface.xMin = xMin;
			surface.xMax = xMax;
			surface.yMin = 0.0;
			surface.yMax = 1.0;
			surface.setSquarePolygon();
			surface.zMin = zMin;
			surface.zMax = zMin;
			surface.boundaryConditionType = Surface::CONSTANT_TEMPERATURE;
			surface.orientation = Surface::Z_NEG;
			surface.temperature = deepGroundTemperature;
			surfaces.push_back(surface);
		}
		else if (deepGroundBoundary == DGB_ZERO_FLUX)
		{
			Surface surface;
			surface.name = "Deep Ground";
			surface.xMin = xMin;
			surface.xMax = xMax;
			surface.yMin = 0.0;
			surface.yMax = 1.0;
			surface.setSquarePolygon();
			surface.zMin = zMin;
			surface.zMax = zMin;
			surface.boundaryConditionType = Surface::ZERO_FLUX;
			surface.orientation = Surface::Z_NEG;
			surfaces.push_back(surface);
		}

		// Slab
		{
			Surface surface;
			surface.name = "Slab Interior";
			surface.xMin = xMin;
			surface.xMax = xRef + xyPerimeterSurface;
			surface.yMin = 0.0;
			surface.yMax = 1.0;
			surface.setSquarePolygon();
			surface.zMin = zSlab;
			surface.zMax = zSlab;
			surface.boundaryConditionType = Surface::INTERIOR_FLUX;
			surface.orientation = Surface::Z_POS;
			surface.emissivity = wall.interiorEmissivity;
			surfaces.push_back(surface);
		}
		if (hasPerimeterSurface)
		{
			Surface surface;
			surface.name = "Slab Perimeter";
			surface.xMin = xRef + xyPerimeterSurface;
			surface.xMax = xRef + xyWallInterior;
			surface.yMin = 0.0;
			surface.yMax = 1.0;
			surface.setSquarePolygon();
			surface.zMin = zSlab;
			surface.zMax = zSlab;
			surface.boundaryConditionType = Surface::INTERIOR_FLUX;
			surface.orientation = Surface::Z_POS;
			surface.emissivity = wall.interiorEmissivity;
			surfaces.push_back(surface);
		}

		// Grade
		{
			Surface surface;
			surface.name = "Grade";
			surface.xMin = xRef + xyWallExterior;
			surface.xMax = xMax;
			surface.yMin = 0.0;
			surface.yMax = 1.0;
			surface.setSquarePolygon();
			surface.zMin = zGrade;
			surface.zMax = zGrade;
			surface.boundaryConditionType = Surface::EXTERIOR_FLUX;
			surface.orientation = Surface::Z_POS;
			surface.emissivity = soilEmissivity;
			surface.absorptivity = soilAbsorptivity;
			surfaces.push_back(surface);
		}

		if(excavationDepth > 0.0)
		{
			// Interior Air Top Surface
			Surface surface;
			surface.name = "Interior Air Top";
			surface.xMin = xMin;
			surface.xMax = xRef + xyWallInterior;
			surface.yMin = 0.0;
			surface.yMax = 1.0;
			surface.setSquarePolygon();
			surface.zMin = zMax;
			surface.zMax = zMax;
			surface.boundaryConditionType = Surface::INTERIOR_TEMPERATURE;
			surface.orientation = Surface::Z_POS;
			surfaces.push_back(surface);
		}

		if(zMax > 0.0)
		{
			// Exterior Air Top Surface
			Surface surface;
			surface.name = "Exterior Air Top";
			surface.xMin = xRef + xyWallExterior;
			surface.xMax = xMax;
			surface.yMin = 0.0;
			surface.yMax = 1.0;
			surface.setSquarePolygon();
			surface.zMin = zMax;
			surface.zMax = zMax;
			surface.boundaryConditionType = Surface::EXTERIOR_TEMPERATURE;
			surface.orientation = Surface::Z_POS;
			surfaces.push_back(surface);
		}

		if (wallTopBoundary == WTB_LINEAR_DT)
		{
			// Wall Top
			if(hasWall)
			{
				double position = 0.0;
				double Tin = indoorAirTemperature;
				std::size_t N = xyWallExterior/mesh.minCellDim;
				double temperature = Tin - (1.0/N)/2*wallTopTemperatureDifference;

				for (std::size_t n = 1; n <= N; n++)
				{
					Surface surface;
					surface.name = "Wall Top";
					surface.xMin = xRef + position;
					surface.xMax = xRef + position + xyWallExterior/N;
					surface.yMin = 0.0;
					surface.yMax = 1.0;
					surface.setSquarePolygon();
					surface.zMin = zMax;
					surface.zMax = zMax;
					surface.boundaryConditionType = Surface::CONSTANT_TEMPERATURE;
					surface.orientation = Surface::Z_POS;
					surface.temperature = temperature;
					surfaces.push_back(surface);

					position += xyWallExterior/N;
					temperature -= (1.0/N)*wallTopTemperatureDifference;

				}
			}
		}
		else
		{
			// Wall Top
			if(hasWall)
			{
				Surface surface;
				surface.name = "Wall Top";
				surface.xMin = xRef + xyWallInterior;
				surface.xMax = xRef + xyWallExterior;
				surface.yMin = 0.0;
				surface.yMax = 1.0;
				surface.setSquarePolygon();
				surface.zMin = zMax;
				surface.zMax = zMax;
				surface.boundaryConditionType = Surface::ZERO_FLUX;
				surface.orientation = Surface::Z_POS;
				surfaces.push_back(surface);
			}
		}

		// Interior Horizontal Insulation
		if (hasInteriorHorizontalInsulation)
		{
			Block block;
			block.material = interiorHorizontalInsulation.layer.material;
			block.blockType = Block::SOLID;
			block.xMin = xRef + xyIntHIns;
			block.xMax = xRef;
			block.yMin = 0.0;
			block.yMax = 1.0;
			block.setSquarePolygon();
			block.zMin = zIntHIns;
			block.zMax = zIntHIns + interiorHorizontalInsulation.layer.thickness;
			blocks.push_back(block);
		}

		if (hasSlab)
		{
			// Foundation Slab
			double zPosition = zSlabBottom;

			for (size_t n = 0; n < slab.layers.size(); n++)
			{
				Block block;
				block.material = slab.layers[n].material;
				block.blockType = Block::SOLID;
				block.xMin = xMin;
				block.xMax = xRef;
				block.yMin = 0.0;
				block.yMax = 1.0;
				block.setSquarePolygon();
				block.zMin = zPosition;
				block.zMax = zPosition + slab.layers[n].thickness;
				blocks.push_back(block);
				zPosition = block.zMax;
			}

		}

		// Interior Vertical Insulation
		if (hasInteriorVerticalInsulation)
		{
			Block block;
			block.material = interiorVerticalInsulation.layer.material;
			block.blockType = Block::SOLID;
			block.xMin = xRef + xyWallInterior;
			block.xMax = xRef;
			block.yMin = 0.0;
			block.yMax = 1.0;
			block.setSquarePolygon();
			block.zMin = zIntVIns;
			block.zMax = zMax;
			blocks.push_back(block);
		}

		// Indoor Air
		{
			Block block;
			block.material = air;
			block.blockType = Block::INTERIOR_AIR;
			block.xMin = xMin;
			block.xMax = xRef + xyWallInterior;
			block.yMin = 0.0;
			block.yMax = 1.0;
			block.setSquarePolygon();
			block.zMin = zSlab;
			block.zMax = zMax;
			blocks.push_back(block);
		}

		if (hasWall)
		{
			double xPosition = xRef;

			// Foundation Wall
			for (int n = wall.layers.size() - 1; n >= 0; n--)
			{
				Block block;
				block.material = wall.layers[n].material;
				block.blockType = Block::SOLID;
				block.xMin = xPosition;
				block.xMax = xPosition + wall.layers[n].thickness;
				block.yMin = 0.0;
				block.yMax = 1.0;
				block.setSquarePolygon();
				block.zMin = zWall;
				block.zMax = zMax;
				xPosition = block.xMax;
				blocks.push_back(block);
			}
		}

		// Exterior Vertical Insulation
		if (hasExteriorVerticalInsulation)
		{
			Block block;
			block.material = exteriorVerticalInsulation.layer.material;
			block.blockType = Block::SOLID;
			block.xMin = xRef + wall.totalWidth();
			block.xMax = xRef + xyWallExterior;
			block.yMin = 0.0;
			block.yMax = 1.0;
			block.setSquarePolygon();
			block.zMin = zExtVIns;
			block.zMax = zMax;
			blocks.push_back(block);
		}

		// Exterior Horizontal Insulation
		if (hasExteriorHorizontalInsulation)
		{
			Block block;
			block.material = exteriorHorizontalInsulation.layer.material;
			block.blockType = Block::SOLID;
			block.xMin = xRef + wall.totalWidth();
			block.xMax = xRef + xyExtHIns;
			block.yMin = 0.0;
			block.yMax = 1.0;
			block.setSquarePolygon();
			block.zMin = zExtHIns;
			block.zMax = zExtHIns + exteriorHorizontalInsulation.layer.thickness;
			blocks.push_back(block);
		}

		// Exterior Air
		{
			Block block;
			block.material = air;
			block.blockType = Block::EXTERIOR_AIR;
			block.xMin = xRef + xyWallExterior;
			block.xMax = xMax;
			block.yMin = 0.0;
			block.yMax = 1.0;
			block.setSquarePolygon();
			block.zMin = zGrade;
			block.zMax = zMax;
			blocks.push_back(block);
		}

		// Set range types

		RangeType xInteriorRange;
		xInteriorRange.range.first = xMin;
		xInteriorRange.range.second = xRef + xyNearInt;
		xInteriorRange.type = RangeType::MIN_INTERIOR;
		xRanges.ranges.push_back(xInteriorRange);

		RangeType xNearRange;
		xNearRange.range.first = xRef + xyNearInt;
		xNearRange.range.second = xRef + xyNearExt;
		xNearRange.type = RangeType::NEAR;
		xRanges.ranges.push_back(xNearRange);

		RangeType xExteriorRange;
		xExteriorRange.range.first = xRef + xyNearExt;
		xExteriorRange.range.second = xMax;
		xExteriorRange.type = RangeType::MAX_EXTERIOR;
		xRanges.ranges.push_back(xExteriorRange);

	}
	else if(coordinateSystem == CS_3D)
	{
		// TODO 3D
		Box boundingBox;
		boost::geometry::envelope(polygon, boundingBox);

		double xMinBB = boundingBox.min_corner().get<0>();
		double yMinBB = boundingBox.min_corner().get<1>();

		double xMaxBB = boundingBox.max_corner().get<0>();
		double yMaxBB = boundingBox.max_corner().get<1>();

		std::size_t nV = polygon.outer().size();

		double xMin = xMinBB - farFieldWidth;
		double yMin = yMinBB - farFieldWidth;

		double xMax = xMaxBB + farFieldWidth;
		double yMax = yMaxBB + farFieldWidth;

		if(excavationDepth > 0.0)
		{
			// Interior Wall Surface
			{
				Polygon poly;
				poly = offset(polygon, xyWallInterior);

				for (std::size_t v = 0; v < nV; v++)
				{
					Surface surface;
					surface.name = "Interior Wall";
					surface.xMin = getXmin(poly,v);
					surface.xMax = getXmax(poly,v);
					surface.yMin = getYmin(poly,v);
					surface.yMax = getYmax(poly,v);
					surface.setSquarePolygon();
					surface.zMin = zSlab;
					surface.zMax = zMax;
					surface.boundaryConditionType = Surface::INTERIOR_FLUX;
					switch (getDirectionOut(poly,v))
					{
					case geom::Y_POS:
						surface.orientation = Surface::X_POS;
						break;
					case geom::X_POS:
						surface.orientation = Surface::Y_NEG;
						break;
					case geom::Y_NEG:
						surface.orientation = Surface::X_NEG;
						break;
					case geom::X_NEG:
						surface.orientation = Surface::Y_POS;
						break;
					}
					surface.emissivity = wall.interiorEmissivity;
					surfaces.push_back(surface);
				}
			}

		}

		if(zMax > 0.0)
		{
			// Exterior Wall Surface
			{
				Polygon poly;
				poly = offset(polygon, xyWallExterior);

				for (std::size_t v = 0; v < nV; v++)
				{
					Surface surface;
					surface.name = "Exterior Wall";
					surface.xMin = getXmin(poly,v);
					surface.xMax = getXmax(poly,v);
					surface.yMin = getYmin(poly,v);
					surface.yMax = getYmax(poly,v);
					surface.setSquarePolygon();
					surface.boundaryConditionType = Surface::EXTERIOR_FLUX;
					switch (getDirectionOut(poly,v))
					{
					case geom::Y_POS:
						surface.orientation = Surface::X_NEG;
						break;
					case geom::X_POS:
						surface.orientation = Surface::Y_POS;
						break;
					case geom::Y_NEG:
						surface.orientation = Surface::X_POS;
						break;
					case geom::X_NEG:
						surface.orientation = Surface::Y_NEG;
						break;
					}
					surface.emissivity = wall.exteriorEmissivity;
					surface.absorptivity = wall.exteriorAbsorptivity;
					surfaces.push_back(surface);
				}
			}
		}

		// Far Field Surfaces
		{
			// X Min
			Surface surface;
			surface.name = "Far Field";
			surface.xMin =	xMin;
			surface.xMax =	xMin;
			surface.yMin =	yMin;
			surface.yMax =	yMax;
			surface.setSquarePolygon();
			surface.zMin =	zMin;
			surface.zMax =	zGrade;
			surface.boundaryConditionType = Surface::ZERO_FLUX;
			surface.orientation =	Surface::X_NEG;
			surfaces.push_back(surface);
		}

		{
			// X Max
			Surface surface;
			surface.name = "Far Field";
			surface.xMin =	xMax;
			surface.xMax =	xMax;
			surface.yMin =	yMin;
			surface.yMax =	yMax;
			surface.setSquarePolygon();
			surface.zMin =	zMin;
			surface.zMax =	zGrade;
			surface.boundaryConditionType = Surface::ZERO_FLUX;
			surface.orientation =	Surface::X_POS;
			surfaces.push_back(surface);
		}

		{
			// Y Min
			Surface surface;
			surface.name = "Far Field";
			surface.xMin =	xMin;
			surface.xMax =	xMax;
			surface.yMin =	yMin;
			surface.yMax =	yMin;
			surface.setSquarePolygon();
			surface.zMin =	zMin;
			surface.zMax =	zGrade;
			surface.boundaryConditionType = Surface::ZERO_FLUX;
			surface.orientation =	Surface::Y_NEG;
			surfaces.push_back(surface);
		}

		{
			// Y Max
			Surface surface;
			surface.name = "Far Field";
			surface.xMin =	xMin;
			surface.xMax =	xMax;
			surface.yMin =	yMax;
			surface.yMax =	yMax;
			surface.setSquarePolygon();
			surface.zMin =	zMin;
			surface.zMax =	zGrade;
			surface.boundaryConditionType = Surface::ZERO_FLUX;
			surface.orientation =	Surface::Y_POS;
			surfaces.push_back(surface);
		}

		// Deep ground surface
		if (deepGroundBoundary == DGB_CONSTANT_TEMPERATURE ||
			deepGroundBoundary == DGB_AUTO)
		{
			Surface surface;
			surface.name = "Deep Ground";
			surface.xMin = xMin;
			surface.xMax = xMax;
			surface.yMin = yMin;
			surface.yMax = yMax;
			surface.setSquarePolygon();
			surface.zMin = zMin;
			surface.zMax = zMin;
			surface.boundaryConditionType = Surface::CONSTANT_TEMPERATURE;
			surface.orientation = Surface::Z_NEG;
			surface.temperature = deepGroundTemperature;
			surfaces.push_back(surface);
		}
		else if (deepGroundBoundary == DGB_ZERO_FLUX)
		{
			Surface surface;
			surface.name = "Deep Ground";
			surface.xMin = xMin;
			surface.xMax = xMax;
			surface.yMin = yMin;
			surface.yMax = yMax;
			surface.setSquarePolygon();
			surface.zMin = zMin;
			surface.zMax = zMin;
			surface.boundaryConditionType = Surface::ZERO_FLUX;
			surface.orientation = Surface::Z_NEG;
			surfaces.push_back(surface);
		}

		// Slab
		{
			Polygon poly;
			poly = offset(polygon, xyPerimeterSurface);

			Surface surface;
			surface.name = "Slab Interior";
			surface.polygon = poly;
			surface.zMin = zSlab;
			surface.zMax = zSlab;
			surface.boundaryConditionType = Surface::INTERIOR_FLUX;
			surface.orientation = Surface::Z_POS;
			surface.emissivity = wall.interiorEmissivity;
			surfaces.push_back(surface);
		}
		if (hasPerimeterSurface)
		{
			Polygon poly;
			poly = offset(polygon, xyWallInterior);

			Polygon temp;
			temp = offset(polygon, xyPerimeterSurface);
			Ring ring;
			boost::geometry::convert(temp, ring);
			boost::geometry::reverse(ring);

			poly.inners().push_back(ring);

			Surface surface;
			surface.name = "Slab Perimeter";
			surface.polygon = poly;
			surface.zMin = zSlab;
			surface.zMax = zSlab;
			surface.boundaryConditionType = Surface::INTERIOR_FLUX;
			surface.orientation = Surface::Z_POS;
			surface.emissivity = wall.interiorEmissivity;
			surfaces.push_back(surface);
		}

		// Grade
		{
			Polygon poly;
			poly = offset(polygon, xyWallExterior);

			Ring ring;
			boost::geometry::convert(poly, ring);
			boost::geometry::reverse(ring);

			Surface surface;
			surface.name = "Grade";
			surface.xMin = xMin;
			surface.xMax = xMax;
			surface.yMin = yMin;
			surface.yMax = yMax;
			surface.setSquarePolygon();
			surface.polygon.inners().push_back(ring);
			surface.zMin = zGrade;
			surface.zMax = zGrade;
			surface.boundaryConditionType = Surface::EXTERIOR_FLUX;
			surface.orientation = Surface::Z_POS;
			surface.emissivity = soilEmissivity;
			surface.absorptivity = soilAbsorptivity;
			surfaces.push_back(surface);
		}

		if(excavationDepth > 0.0)
		{
			// Interior Air Top Surface
			Polygon poly;
			poly = offset(polygon, xyWallInterior);

			Surface surface;
			surface.name = "Interior Air Top";
			surface.polygon = poly;
			surface.zMin = zMax;
			surface.zMax = zMax;
			surface.boundaryConditionType = Surface::INTERIOR_TEMPERATURE;
			surface.orientation = Surface::Z_POS;
			surfaces.push_back(surface);
		}

		if(zMax > 0.0)
		{
			// Exterior Air Top Surface
			Polygon poly;
			poly = offset(polygon, xyWallExterior);

			Ring ring;
			boost::geometry::convert(poly, ring);
			boost::geometry::reverse(ring);

			Surface surface;
			surface.name = "Exterior Air Top";
			surface.xMin = xMin;
			surface.xMax = xMax;
			surface.yMin = yMin;
			surface.yMax = yMax;
			surface.setSquarePolygon();
			surface.polygon.inners().push_back(ring);
			surface.zMin = zMax;
			surface.zMax = zMax;
			surface.boundaryConditionType = Surface::EXTERIOR_TEMPERATURE;
			surface.orientation = Surface::Z_POS;
			surfaces.push_back(surface);
		}


		if (wallTopBoundary == WTB_LINEAR_DT)
		{
			// Wall Top
			if(hasWall)
			{
				double position = 0.0;
				double Tin = indoorAirTemperature;
				std::size_t N = xyWallExterior/mesh.minCellDim;
				double temperature = Tin - (1.0/N)/2*wallTopTemperatureDifference;

				for (std::size_t n = 1; n <= N; n++)
				{

					Polygon poly;
					poly = offset(polygon, position + xyWallExterior/N);

					Polygon temp;
					temp = offset(polygon, position);
					Ring ring;
					boost::geometry::convert(temp, ring);
					boost::geometry::reverse(ring);

					poly.inners().push_back(ring);

					Surface surface;
					surface.name = "Wall Top";
					surface.polygon = poly;
					surface.zMin = zMax;
					surface.zMax = zMax;
					surface.boundaryConditionType = Surface::CONSTANT_TEMPERATURE;
					surface.orientation = Surface::Z_POS;
					surface.temperature = temperature;
					surfaces.push_back(surface);

					position += xyWallExterior/N;
					temperature -= (1.0/N)*wallTopTemperatureDifference;
				}
			}
		}
		else
		{
			// Wall Top
			if(hasWall)
			{
				Polygon poly;
				poly = offset(polygon, xyWallExterior);

				Polygon temp;
				temp = offset(polygon, xyWallInterior);
				Ring ring;
				boost::geometry::convert(temp, ring);
				boost::geometry::reverse(ring);

				poly.inners().push_back(ring);

				Surface surface;
				surface.name = "Wall Top";
				surface.polygon = poly;
				surface.zMin = zMax;
				surface.zMax = zMax;
				surface.boundaryConditionType = Surface::ZERO_FLUX;
				surface.orientation = Surface::Z_POS;
				surfaces.push_back(surface);
			}
		}

		// Interior Horizontal Insulation
		if (hasInteriorHorizontalInsulation)
		{
			Polygon poly;
			poly = polygon;

			Polygon temp;
			temp = offset(polygon, xyIntHIns);
			Ring ring;
			boost::geometry::convert(temp, ring);
			boost::geometry::reverse(ring);

			poly.inners().push_back(ring);

			Block block;
			block.material = interiorHorizontalInsulation.layer.material;
			block.blockType = Block::SOLID;
			block.polygon = poly;
			block.zMin = zIntHIns;
			block.zMax = zIntHIns + interiorHorizontalInsulation.layer.thickness;
			blocks.push_back(block);
		}

		if (hasSlab)
		{
			// Foundation Slab
			double zPosition = zSlabBottom;

			for (size_t n = 0; n < slab.layers.size(); n++)
			{
				Block block;
				block.material = slab.layers[n].material;
				block.blockType = Block::SOLID;
				block.polygon = polygon;
				block.zMin = zPosition;
				block.zMax = zPosition + slab.layers[n].thickness;
				blocks.push_back(block);
				zPosition = block.zMax;
			}
		}

		// Interior Vertical Insulation
		if (hasInteriorVerticalInsulation)
		{
			Polygon poly;
			poly = polygon;

			Polygon temp;
			temp = offset(polygon, xyWallInterior);
			Ring ring;
			boost::geometry::convert(temp, ring);
			boost::geometry::reverse(ring);

			poly.inners().push_back(ring);
			Block block;
			block.material = interiorVerticalInsulation.layer.material;
			block.blockType = Block::SOLID;
			block.polygon = poly;
			block.zMin = zIntVIns;
			block.zMax = zMax;
			blocks.push_back(block);
		}

		// Indoor Air
		{
			Polygon poly;
			poly = offset(polygon, xyWallInterior);

			Block block;
			block.material = air;
			block.blockType = Block::INTERIOR_AIR;
			block.polygon = poly;
			block.zMin = zSlab;
			block.zMax = zMax;
			blocks.push_back(block);
		}

		if (hasWall)
		{
			double xyPosition = 0.0;

			// Foundation Wall
			for (int n = wall.layers.size() - 1; n >= 0; n--)
			{
				Polygon poly;
				poly = offset(polygon, xyPosition + wall.layers[n].thickness);

				Polygon temp;
				temp = offset(polygon, xyPosition);
				Ring ring;
				boost::geometry::convert(temp, ring);
				boost::geometry::reverse(ring);

				poly.inners().push_back(ring);

				Block block;
				block.material = wall.layers[n].material;
				block.blockType = Block::SOLID;
				block.polygon = poly;
				block.zMin = zWall;
				block.zMax = zMax;
				xyPosition += wall.layers[n].thickness;
				blocks.push_back(block);
			}
		}

		// Exterior Vertical Insulation
		if (hasExteriorVerticalInsulation)
		{
			Polygon poly;
			poly = offset(polygon, xyWallExterior);

			Polygon temp;
			temp = offset(polygon, wall.totalWidth());
			Ring ring;
			boost::geometry::convert(temp, ring);
			boost::geometry::reverse(ring);

			poly.inners().push_back(ring);

			Block block;
			block.material = exteriorVerticalInsulation.layer.material;
			block.blockType = Block::SOLID;
			block.polygon = poly;
			block.zMin = zExtVIns;
			block.zMax = zMax;
			blocks.push_back(block);
		}

		// Exterior Horizontal Insulation
		if (hasExteriorHorizontalInsulation)
		{
			Polygon poly;
			poly = offset(polygon, xyExtHIns);

			Polygon temp;
			temp = offset(polygon, wall.totalWidth());
			Ring ring;
			boost::geometry::convert(temp, ring);
			boost::geometry::reverse(ring);

			poly.inners().push_back(ring);

			Block block;
			block.material = exteriorHorizontalInsulation.layer.material;
			block.blockType = Block::SOLID;
			block.polygon = poly;
			block.zMin = zExtHIns;
			block.zMax = zExtHIns + exteriorHorizontalInsulation.layer.thickness;
			blocks.push_back(block);
		}

		// Exterior Air
		{
			Polygon poly;
			poly = offset(polygon, xyWallExterior);

			Ring ring;
			boost::geometry::convert(poly, ring);
			boost::geometry::reverse(ring);

			Block block;
			block.material = air;
			block.blockType = Block::EXTERIOR_AIR;
			block.xMin = xMin;
			block.xMax = xMax;
			block.yMin = yMin;
			block.yMax = yMax;
			block.setSquarePolygon();
			block.polygon.inners().push_back(ring);
			block.zMin = zGrade;
			block.zMax = zMax;
			blocks.push_back(block);
		}

		// Set x and y near ranges
		std::vector<RangeType> xNearRanges;
		std::vector<RangeType> yNearRanges;

		for (std::size_t v = 0; v < nV; v++)
		{
			double x = polygon.outer()[v].get<0>();
			double y = polygon.outer()[v].get<1>();

			switch (getDirectionOut(polygon,v))
			{
			case geom::Y_POS:
				{
					RangeType range;
					range.range.first = x - xyNearExt;
					range.range.second = x - xyNearInt;
					range.type = RangeType::NEAR;
					xNearRanges.push_back(range);
				}
				break;

			case geom::Y_NEG:
				{
					RangeType range;
					range.range.first = x + xyNearInt;
					range.range.second = x + xyNearExt;
					range.type = RangeType::NEAR;
					xNearRanges.push_back(range);
				}
				break;
			case geom::X_POS:
				{
					RangeType range;
					range.range.first = y + xyNearInt;
					range.range.second = y + xyNearExt;
					range.type = RangeType::NEAR;
					yNearRanges.push_back(range);
				}
				break;
			case geom::X_NEG:
				{
					RangeType range;
					range.range.first = y - xyNearExt;
					range.range.second = y - xyNearInt;
					range.type = RangeType::NEAR;
					yNearRanges.push_back(range);
				}
				break;
			}
		}

		// Set X range types
		sort(xNearRanges.begin(), xNearRanges.end(), compareRanges);

		// Merge overlapping ranges
		for (std::size_t r = 1; r < xNearRanges.size(); r++)
		{
			if (isLessOrEqual(xNearRanges[r].range.first,xNearRanges[r-1].range.second))
			{
				xNearRanges[r-1].range.second = xNearRanges[r].range.second;
				xNearRanges.erase(xNearRanges.begin() + r-1);
				r -= 1;
			}
		}

		RangeType xMinExteriorRange;
		xMinExteriorRange.range.first = xMin;
		xMinExteriorRange.range.second = xNearRanges[0].range.first;
		xMinExteriorRange.type = RangeType::MIN_EXTERIOR;
		xRanges.ranges.push_back(xMinExteriorRange);

		for (std::size_t r = 0; r < xNearRanges.size(); r++)
		{
			if (r == 0)
			{
				RangeType xNearRange;
				xNearRange.range.first = xNearRanges[r].range.first;
				xNearRange.range.second = xNearRanges[r].range.second;
				xNearRange.type = RangeType::NEAR;
				xRanges.ranges.push_back(xNearRange);
			}
			else
			{
				RangeType xInteriorRange;
				xInteriorRange.range.first = xNearRanges[r-1].range.second;
				xInteriorRange.range.second = xNearRanges[r].range.first;
				xInteriorRange.type = RangeType::MID_INTERIOR;
				xRanges.ranges.push_back(xInteriorRange);

				RangeType xNearRange;
				xNearRange.range.first = xNearRanges[r].range.first;
				xNearRange.range.second = xNearRanges[r].range.second;
				xNearRange.type = RangeType::NEAR;
				xRanges.ranges.push_back(xNearRange);
			}
		}

		RangeType xMaxExteriorRange;
		xMaxExteriorRange.range.first = xNearRanges[xNearRanges.size() - 1].range.second;
		xMaxExteriorRange.range.second = xMax;
		xMaxExteriorRange.type = RangeType::MAX_EXTERIOR;
		xRanges.ranges.push_back(xMaxExteriorRange);


		// Set Y range types
		sort(yNearRanges.begin(), yNearRanges.end(), compareRanges);

		// Merge overlapping ranges
		for (std::size_t r = 1; r < yNearRanges.size(); r++)
		{
			if (isLessOrEqual(yNearRanges[r].range.first,yNearRanges[r-1].range.second))
			{
				yNearRanges[r-1].range.second = yNearRanges[r].range.second;
				yNearRanges.erase(yNearRanges.begin() + r-1);
				r -= 1;
			}
		}

		RangeType yMinExteriorRange;
		yMinExteriorRange.range.first = yMin;
		yMinExteriorRange.range.second = yNearRanges[0].range.first;
		yMinExteriorRange.type = RangeType::MIN_EXTERIOR;
		yRanges.ranges.push_back(yMinExteriorRange);

		for (std::size_t r = 0; r < yNearRanges.size(); r++)
		{
			if (r == 0)
			{
				RangeType yNearRange;
				yNearRange.range.first = yNearRanges[r].range.first;
				yNearRange.range.second = yNearRanges[r].range.second;
				yNearRange.type = RangeType::NEAR;
				yRanges.ranges.push_back(yNearRange);
			}
			else
			{
				RangeType yInteriorRange;
				yInteriorRange.range.first = yNearRanges[r-1].range.second;
				yInteriorRange.range.second = yNearRanges[r].range.first;
				yInteriorRange.type = RangeType::MID_INTERIOR;
				yRanges.ranges.push_back(yInteriorRange);

				RangeType yNearRange;
				yNearRange.range.first = yNearRanges[r].range.first;
				yNearRange.range.second = yNearRanges[r].range.second;
				yNearRange.type = RangeType::NEAR;
				yRanges.ranges.push_back(yNearRange);
			}
		}

		RangeType yMaxExteriorRange;
		yMaxExteriorRange.range.first = yNearRanges[yNearRanges.size() - 1].range.second;
		yMaxExteriorRange.range.second = yMax;
		yMaxExteriorRange.type = RangeType::MAX_EXTERIOR;
		yRanges.ranges.push_back(yMaxExteriorRange);

	}
	else if(coordinateSystem == CS_3D_SYMMETRY)
	{
		// TODO: 3D Symmetric
		bool isXSymm = isXSymmetric(polygon);
		bool isYSymm = isYSymmetric(polygon);

		// Change polygon to minimum symmetric unit
		Polygon symmetricPoly;
		symmetricPoly = symmetricUnit(polygon);

		Box boundingBox;
		boost::geometry::envelope(symmetricPoly, boundingBox);

		double xMinBB = boundingBox.min_corner().get<0>();
		double yMinBB = boundingBox.min_corner().get<1>();

		double xMaxBB = boundingBox.max_corner().get<0>();
		double yMaxBB = boundingBox.max_corner().get<1>();

		double xMin;
		double yMin;

		if (isXSymm)
			xMin = xMinBB;
		else
			xMin = xMinBB - farFieldWidth;

		if (isYSymm)
			yMin = yMinBB;
		else
			yMin = yMinBB - farFieldWidth;

		double xMax = xMaxBB + farFieldWidth;
		double yMax = yMaxBB + farFieldWidth;

		Box domainBox(Point(xMin,yMin),Point(xMax,yMax));

		std::size_t nV = polygon.outer().size();

		if(excavationDepth > 0.0)
		{
			// Interior Wall Surface
			{
				Polygon tempPoly;
				tempPoly = offset(polygon, xyWallInterior);
				MultiPolygon poly;
				boost::geometry::intersection(domainBox,tempPoly,poly);

				for (std::size_t v = 0; v < nV; v++)
				{

					if (!((isEqual(getXmin(poly[0],v),xMin) && isEqual(getXmax(poly[0],v),xMin)) ||
						  (isEqual(getYmin(poly[0],v),yMin) && isEqual(getYmax(poly[0],v),yMin))))
					{
						Surface surface;
						surface.name = "Interior Wall";
						surface.xMin = getXmin(poly[0],v);
						surface.xMax = getXmax(poly[0],v);
						surface.yMin = getYmin(poly[0],v);
						surface.yMax = getYmax(poly[0],v);
						surface.setSquarePolygon();
						surface.zMin = zSlab;
						surface.zMax = zMax;
						surface.boundaryConditionType = Surface::INTERIOR_FLUX;
						switch (getDirectionOut(poly[0],v))
						{
						case geom::Y_POS:
							surface.orientation = Surface::X_POS;
							break;
						case geom::X_POS:
							surface.orientation = Surface::Y_NEG;
							break;
						case geom::Y_NEG:
							surface.orientation = Surface::X_NEG;
							break;
						case geom::X_NEG:
							surface.orientation = Surface::Y_POS;
							break;
						}
						surface.emissivity = wall.interiorEmissivity;
						surfaces.push_back(surface);
					}
				}
			}

		}

		if(zMax > 0.0)
		{
			// Exterior Wall Surface
			{
				Polygon tempPoly;
				tempPoly = offset(polygon, xyWallExterior);
				MultiPolygon poly;
				boost::geometry::intersection(domainBox,tempPoly,poly);

				for (std::size_t v = 0; v < nV; v++)
				{
					if (!((isEqual(getXmin(poly[0],v),xMin) && isEqual(getXmax(poly[0],v),xMin)) ||
						  (isEqual(getYmin(poly[0],v),yMin) && isEqual(getYmax(poly[0],v),yMin))))
					{
						Surface surface;
						surface.name = "Exterior Wall";
						surface.xMin = getXmin(poly[0],v);
						surface.xMax = getXmax(poly[0],v);
						surface.yMin = getYmin(poly[0],v);
						surface.yMax = getYmax(poly[0],v);
						surface.setSquarePolygon();
						surface.boundaryConditionType = Surface::EXTERIOR_FLUX;
						switch (getDirectionOut(poly[0],v))
						{
						case geom::Y_POS:
							surface.orientation = Surface::X_NEG;
							break;
						case geom::X_POS:
							surface.orientation = Surface::Y_POS;
							break;
						case geom::Y_NEG:
							surface.orientation = Surface::X_POS;
							break;
						case geom::X_NEG:
							surface.orientation = Surface::Y_NEG;
							break;
						}
						surface.emissivity = wall.exteriorEmissivity;
						surface.absorptivity = wall.exteriorAbsorptivity;
						surfaces.push_back(surface);
					}
				}
			}
		}

		// Symmetry & Far Field Surfaces
		{
			// X Min
			Surface surface;
			surface.name = "Far Field";
			surface.xMin =	xMin;
			surface.xMax =	xMin;
			surface.yMin =	yMin;
			surface.yMax =	yMax;
			surface.setSquarePolygon();
			surface.zMin =	zMin;
			surface.zMax =	zGrade;
			surface.boundaryConditionType = Surface::ZERO_FLUX;
			surface.orientation = Surface::X_NEG;
			surfaces.push_back(surface);
		}

		{
			// X Max
			Surface surface;
			surface.name = "Far Field";
			surface.xMin =	xMax;
			surface.xMax =	xMax;
			surface.yMin =	yMin;
			surface.yMax =	yMax;
			surface.setSquarePolygon();
			surface.zMin =	zMin;
			surface.zMax =	zGrade;
			surface.boundaryConditionType = Surface::ZERO_FLUX;
			surface.orientation = Surface::X_POS;
			surfaces.push_back(surface);
		}

		{
			// Y Min
			Surface surface;
			surface.name = "Far Field";
			surface.xMin =	xMin;
			surface.xMax =	xMax;
			surface.yMin =	yMin;
			surface.yMax =	yMin;
			surface.setSquarePolygon();
			surface.zMin =	zMin;
			surface.zMax =	zGrade;
			surface.boundaryConditionType = Surface::ZERO_FLUX;
			surface.orientation =	Surface::Y_NEG;
			surfaces.push_back(surface);
		}

		{
			// Y Max
			Surface surface;
			surface.name = "Far Field";
			surface.xMin =	xMin;
			surface.xMax =	xMax;
			surface.yMin =	yMax;
			surface.yMax =	yMax;
			surface.setSquarePolygon();
			surface.zMin =	zMin;
			surface.zMax =	zGrade;
			surface.boundaryConditionType = Surface::ZERO_FLUX;
			surface.orientation =	Surface::Y_POS;
			surfaces.push_back(surface);
		}

		// Deep ground surface
		if (deepGroundBoundary == DGB_CONSTANT_TEMPERATURE ||
			deepGroundBoundary == DGB_AUTO)
		{
			Surface surface;
			surface.name = "Deep Ground";
			surface.xMin = xMin;
			surface.xMax = xMax;
			surface.yMin = yMin;
			surface.yMax = yMax;
			surface.setSquarePolygon();
			surface.zMin = zMin;
			surface.zMax = zMin;
			surface.boundaryConditionType = Surface::CONSTANT_TEMPERATURE;
			surface.orientation = Surface::Z_NEG;
			surface.temperature = deepGroundTemperature;
			surfaces.push_back(surface);
		}
		else if (deepGroundBoundary == DGB_ZERO_FLUX)
		{
			Surface surface;
			surface.name = "Deep Ground";
			surface.xMin = xMin;
			surface.xMax = xMax;
			surface.yMin = yMin;
			surface.yMax = yMax;
			surface.setSquarePolygon();
			surface.zMin = zMin;
			surface.zMax = zMin;
			surface.boundaryConditionType = Surface::ZERO_FLUX;
			surface.orientation = Surface::Z_NEG;
			surfaces.push_back(surface);
		}

		// Slab
		{
			Polygon tempPoly;
			tempPoly = offset(polygon, xyPerimeterSurface);
			MultiPolygon poly;
			boost::geometry::intersection(domainBox,tempPoly,poly);

			Surface surface;
			surface.name = "Slab Interior";
			surface.polygon = poly[0];
			surface.zMin = zSlab;
			surface.zMax = zSlab;
			surface.boundaryConditionType = Surface::INTERIOR_FLUX;
			surface.orientation = Surface::Z_POS;
			surface.emissivity = wall.interiorEmissivity;
			surfaces.push_back(surface);
		}
		if (hasPerimeterSurface)
		{
			Polygon tempPoly;
			tempPoly = offset(polygon, xyWallInterior);

			Polygon temp;
			temp = offset(polygon, xyPerimeterSurface);
			Ring ring;
			boost::geometry::convert(temp, ring);
			boost::geometry::reverse(ring);

			tempPoly.inners().push_back(ring);

			MultiPolygon poly;
			boost::geometry::intersection(domainBox,tempPoly,poly);

			Surface surface;
			surface.name = "Slab Perimeter";
			surface.polygon = poly[0];
			surface.zMin = zSlab;
			surface.zMax = zSlab;
			surface.boundaryConditionType = Surface::INTERIOR_FLUX;
			surface.orientation = Surface::Z_POS;
			surface.emissivity = wall.interiorEmissivity;
			surfaces.push_back(surface);
		}

		// Grade
		{
			Polygon tempPoly;
			tempPoly = offset(polygon, xyWallExterior);

			MultiPolygon poly;
			boost::geometry::difference(domainBox,tempPoly,poly);

			Surface surface;
			surface.name = "Grade";
			surface.polygon = poly[0];
			surface.zMin = zGrade;
			surface.zMax = zGrade;
			surface.boundaryConditionType = Surface::EXTERIOR_FLUX;
			surface.orientation = Surface::Z_POS;
			surface.emissivity = soilEmissivity;
			surface.absorptivity = soilAbsorptivity;
			surfaces.push_back(surface);
		}

		if(excavationDepth > 0.0)
		{
			// Interior Air Top Surface
			Polygon tempPoly;
			tempPoly = offset(polygon, xyWallInterior);
			MultiPolygon poly;
			boost::geometry::intersection(domainBox,tempPoly,poly);

			Surface surface;
			surface.name = "Interior Air Top";
			surface.polygon = poly[0];
			surface.zMin = zMax;
			surface.zMax = zMax;
			surface.boundaryConditionType = Surface::INTERIOR_TEMPERATURE;
			surface.orientation = Surface::Z_POS;
			surfaces.push_back(surface);
		}

		if(zMax > 0.0)
		{
			// Exterior Air Top Surface
			Polygon tempPoly;
			tempPoly = offset(polygon, xyWallExterior);

			MultiPolygon poly;
			boost::geometry::difference(domainBox,tempPoly,poly);

			Surface surface;
			surface.name = "Exterior Air Top";
			surface.polygon = poly[0];
			surface.zMin = zMax;
			surface.zMax = zMax;
			surface.boundaryConditionType = Surface::EXTERIOR_TEMPERATURE;
			surface.orientation = Surface::Z_POS;
			surfaces.push_back(surface);
		}


		if (wallTopBoundary == WTB_LINEAR_DT)
		{
			// Wall Top
			if(hasWall)
			{
				double position = 0.0;
				double Tin = indoorAirTemperature;
				std::size_t N = xyWallExterior/mesh.minCellDim;
				double temperature = Tin - (1.0/N)/2*wallTopTemperatureDifference;

				for (std::size_t n = 1; n <= N; n++)
				{

					Polygon tempPoly;
					tempPoly = offset(polygon, position + xyWallExterior/N);

					Polygon temp;
					temp = offset(polygon, position);
					Ring ring;
					boost::geometry::convert(temp, ring);
					boost::geometry::reverse(ring);

					tempPoly.inners().push_back(ring);

					MultiPolygon poly;
					boost::geometry::intersection(domainBox,tempPoly,poly);

					Surface surface;
					surface.name = "Wall Top";
					surface.polygon = poly[0];
					surface.zMin = zMax;
					surface.zMax = zMax;
					surface.boundaryConditionType = Surface::CONSTANT_TEMPERATURE;
					surface.orientation = Surface::Z_POS;
					surface.temperature = temperature;
					surfaces.push_back(surface);

					position += xyWallExterior/N;
					temperature -= (1.0/N)*wallTopTemperatureDifference;
				}
			}
		}
		else
		{
			// Wall Top
			if(hasWall)
			{
				Polygon tempPoly;
				tempPoly = offset(polygon, xyWallExterior);

				Polygon temp;
				temp = offset(polygon, xyWallInterior);
				Ring ring;
				boost::geometry::convert(temp, ring);
				boost::geometry::reverse(ring);

				tempPoly.inners().push_back(ring);

				MultiPolygon poly;
				boost::geometry::intersection(domainBox,tempPoly,poly);

				Surface surface;
				surface.name = "Wall Top";
				surface.polygon = poly[0];
				surface.zMin = zMax;
				surface.zMax = zMax;
				surface.boundaryConditionType = Surface::ZERO_FLUX;
				surface.orientation = Surface::Z_POS;
				surfaces.push_back(surface);
			}
		}

		// Interior Horizontal Insulation
		if (hasInteriorHorizontalInsulation)
		{
			Polygon tempPoly;
			tempPoly = polygon;

			Polygon temp;
			temp = offset(polygon, xyIntHIns);
			Ring ring;
			boost::geometry::convert(temp, ring);
			boost::geometry::reverse(ring);

			tempPoly.inners().push_back(ring);

			MultiPolygon poly;
			boost::geometry::intersection(domainBox,tempPoly,poly);

			Block block;
			block.material = interiorHorizontalInsulation.layer.material;
			block.blockType = Block::SOLID;
			block.polygon = poly[0];
			block.zMin = zIntHIns;
			block.zMax = zIntHIns + interiorHorizontalInsulation.layer.thickness;
			blocks.push_back(block);
		}

		if (hasSlab)
		{
			// Foundation Slab
			double zPosition = zSlabBottom;

			Polygon tempPoly;
			tempPoly = polygon;

			MultiPolygon poly;
			boost::geometry::intersection(domainBox,tempPoly,poly);

			for (size_t n = 0; n < slab.layers.size(); n++)
			{
				Block block;
				block.material = slab.layers[n].material;
				block.blockType = Block::SOLID;
				block.polygon = poly[0];
				block.zMin = zPosition;
				block.zMax = zPosition + slab.layers[n].thickness;
				blocks.push_back(block);
				zPosition = block.zMax;
			}
		}

		// Interior Vertical Insulation
		if (hasInteriorVerticalInsulation)
		{
			Polygon tempPoly;
			tempPoly = polygon;

			Polygon temp;
			temp = offset(polygon, xyWallInterior);
			Ring ring;
			boost::geometry::convert(temp, ring);
			boost::geometry::reverse(ring);

			tempPoly.inners().push_back(ring);

			MultiPolygon poly;
			boost::geometry::intersection(domainBox,tempPoly,poly);

			Block block;
			block.material = interiorVerticalInsulation.layer.material;
			block.blockType = Block::SOLID;
			block.polygon = poly[0];
			block.zMin = zIntVIns;
			block.zMax = zMax;
			blocks.push_back(block);
		}

		// Indoor Air
		{
			Polygon tempPoly;
			tempPoly = offset(polygon, xyWallInterior);

			MultiPolygon poly;
			boost::geometry::intersection(domainBox,tempPoly,poly);

			Block block;
			block.material = air;
			block.blockType = Block::INTERIOR_AIR;
			block.polygon = poly[0];
			block.zMin = zSlab;
			block.zMax = zMax;
			blocks.push_back(block);
		}

		if (hasWall)
		{
			double xyPosition = 0.0;

			// Foundation Wall
			for (int n = wall.layers.size() - 1; n >= 0; n--)
			{
				Polygon tempPoly;
				tempPoly = offset(polygon, xyPosition + wall.layers[n].thickness);

				Polygon temp;
				temp = offset(polygon, xyPosition);
				Ring ring;
				boost::geometry::convert(temp, ring);
				boost::geometry::reverse(ring);

				tempPoly.inners().push_back(ring);

				MultiPolygon poly;
				boost::geometry::intersection(domainBox,tempPoly,poly);

				Block block;
				block.material = wall.layers[n].material;
				block.blockType = Block::SOLID;
				block.polygon = poly[0];
				block.zMin = zWall;
				block.zMax = zMax;
				xyPosition += wall.layers[n].thickness;
				blocks.push_back(block);
			}
		}

		// Exterior Vertical Insulation
		if (hasExteriorVerticalInsulation)
		{
			Polygon tempPoly;
			tempPoly = offset(polygon, xyWallExterior);

			Polygon temp;
			temp = offset(polygon, wall.totalWidth());
			Ring ring;
			boost::geometry::convert(temp, ring);
			boost::geometry::reverse(ring);

			tempPoly.inners().push_back(ring);

			MultiPolygon poly;
			boost::geometry::intersection(domainBox,tempPoly,poly);

			Block block;
			block.material = exteriorVerticalInsulation.layer.material;
			block.blockType = Block::SOLID;
			block.polygon = poly[0];
			block.zMin = zExtVIns;
			block.zMax = zMax;
			blocks.push_back(block);
		}

		// Exterior Horizontal Insulation
		if (hasExteriorHorizontalInsulation)
		{
			Polygon tempPoly;
			tempPoly = offset(polygon, xyExtHIns);

			Polygon temp;
			temp = offset(polygon, wall.totalWidth());
			Ring ring;
			boost::geometry::convert(temp, ring);
			boost::geometry::reverse(ring);

			tempPoly.inners().push_back(ring);

			MultiPolygon poly;
			boost::geometry::intersection(domainBox,tempPoly,poly);

			Block block;
			block.material = exteriorHorizontalInsulation.layer.material;
			block.blockType = Block::SOLID;
			block.polygon = poly[0];
			block.zMin = zExtHIns;
			block.zMax = zExtHIns + exteriorHorizontalInsulation.layer.thickness;
			blocks.push_back(block);
		}

		// Exterior Air
		{
			Polygon tempPoly;
			tempPoly = offset(polygon, xyWallExterior);

			MultiPolygon poly;
			boost::geometry::difference(domainBox,tempPoly,poly);

			Block block;
			block.material = air;
			block.blockType = Block::EXTERIOR_AIR;
			block.polygon = poly[0];
			block.zMin = zGrade;
			block.zMax = zMax;
			blocks.push_back(block);
		}

		// Set x and y near ranges
		std::vector<RangeType> xNearRanges;
		std::vector<RangeType> yNearRanges;

		for (std::size_t v = 0; v < symmetricPoly.outer().size(); v++)
		{
			if (!((isEqual(getXmin(symmetricPoly,v),xMin) && isEqual(getXmax(symmetricPoly,v),xMin)) ||
				  (isEqual(getYmin(symmetricPoly,v),yMin) && isEqual(getYmax(symmetricPoly,v),yMin))))
			{
				double x = symmetricPoly.outer()[v].get<0>();
				double y = symmetricPoly.outer()[v].get<1>();

				switch (getDirectionOut(symmetricPoly,v))
				{
				case geom::Y_POS:
					{
						RangeType range;
						range.range.first = x - xyNearExt;
						range.range.second = x - xyNearInt;
						range.type = RangeType::NEAR;
						xNearRanges.push_back(range);
					}
					break;

				case geom::Y_NEG:
					{
						RangeType range;
						range.range.first = x + xyNearInt;
						range.range.second = x + xyNearExt;
						range.type = RangeType::NEAR;
						xNearRanges.push_back(range);
					}
					break;
				case geom::X_POS:
					{
						RangeType range;
						range.range.first = y + xyNearInt;
						range.range.second = y + xyNearExt;
						range.type = RangeType::NEAR;
						yNearRanges.push_back(range);
					}
					break;
				case geom::X_NEG:
					{
						RangeType range;
						range.range.first = y - xyNearExt;
						range.range.second = y - xyNearInt;
						range.type = RangeType::NEAR;
						yNearRanges.push_back(range);
					}
					break;
				}
			}
		}

		// Set X range types
		sort(xNearRanges.begin(), xNearRanges.end(), compareRanges);

		// Merge overlapping ranges
		for (std::size_t r = 1; r < xNearRanges.size(); r++)
		{
			if (isLessOrEqual(xNearRanges[r].range.first,xNearRanges[r-1].range.second))
			{
				xNearRanges[r-1].range.second = xNearRanges[r].range.second;
				xNearRanges.erase(xNearRanges.begin() + r-1);
				r -= 1;
			}
		}

		if (isXSymm)
		{
			RangeType xMinInteriorRange;
			xMinInteriorRange.range.first = xMin;
			xMinInteriorRange.range.second = xNearRanges[0].range.first;
			xMinInteriorRange.type = RangeType::MIN_INTERIOR;
			xRanges.ranges.push_back(xMinInteriorRange);
		}
		else
		{
			RangeType xMinExteriorRange;
			xMinExteriorRange.range.first = xMin;
			xMinExteriorRange.range.second = xNearRanges[0].range.first;
			xMinExteriorRange.type = RangeType::MIN_EXTERIOR;
			xRanges.ranges.push_back(xMinExteriorRange);
		}

		for (std::size_t r = 0; r < xNearRanges.size(); r++)
		{
			if (r == 0)
			{
				RangeType xNearRange;
				xNearRange.range.first = xNearRanges[r].range.first;
				xNearRange.range.second = xNearRanges[r].range.second;
				xNearRange.type = RangeType::NEAR;
				xRanges.ranges.push_back(xNearRange);
			}
			else
			{
				RangeType xInteriorRange;
				xInteriorRange.range.first = xNearRanges[r-1].range.second;
				xInteriorRange.range.second = xNearRanges[r].range.first;
				xInteriorRange.type = RangeType::MID_INTERIOR;
				xRanges.ranges.push_back(xInteriorRange);

				RangeType xNearRange;
				xNearRange.range.first = xNearRanges[r].range.first;
				xNearRange.range.second = xNearRanges[r].range.second;
				xNearRange.type = RangeType::NEAR;
				xRanges.ranges.push_back(xNearRange);
			}
		}

		RangeType xMaxExteriorRange;
		xMaxExteriorRange.range.first = xNearRanges[xNearRanges.size() - 1].range.second;
		xMaxExteriorRange.range.second = xMax;
		xMaxExteriorRange.type = RangeType::MAX_EXTERIOR;
		xRanges.ranges.push_back(xMaxExteriorRange);


		// Set Y range types
		sort(yNearRanges.begin(), yNearRanges.end(), compareRanges);

		// Merge overlapping ranges
		for (std::size_t r = 1; r < yNearRanges.size(); r++)
		{
			if (isLessOrEqual(yNearRanges[r].range.first,yNearRanges[r-1].range.second))
			{
				yNearRanges[r-1].range.second = yNearRanges[r].range.second;
				yNearRanges.erase(yNearRanges.begin() + r-1);
				r -= 1;
			}
		}

		if (isYSymm)
		{
			RangeType yMinInteriorRange;
			yMinInteriorRange.range.first = yMin;
			yMinInteriorRange.range.second = yNearRanges[0].range.first;
			yMinInteriorRange.type = RangeType::MIN_INTERIOR;
			yRanges.ranges.push_back(yMinInteriorRange);
		}
		else
		{
			RangeType yMinExteriorRange;
			yMinExteriorRange.range.first = yMin;
			yMinExteriorRange.range.second = yNearRanges[0].range.first;
			yMinExteriorRange.type = RangeType::MIN_EXTERIOR;
			yRanges.ranges.push_back(yMinExteriorRange);
		}

		for (std::size_t r = 0; r < yNearRanges.size(); r++)
		{
			if (r == 0)
			{
				RangeType yNearRange;
				yNearRange.range.first = yNearRanges[r].range.first;
				yNearRange.range.second = yNearRanges[r].range.second;
				yNearRange.type = RangeType::NEAR;
				yRanges.ranges.push_back(yNearRange);
			}
			else
			{
				RangeType yInteriorRange;
				yInteriorRange.range.first = yNearRanges[r-1].range.second;
				yInteriorRange.range.second = yNearRanges[r].range.first;
				yInteriorRange.type = RangeType::MID_INTERIOR;
				yRanges.ranges.push_back(yInteriorRange);

				RangeType yNearRange;
				yNearRange.range.first = yNearRanges[r].range.first;
				yNearRange.range.second = yNearRanges[r].range.second;
				yNearRange.type = RangeType::NEAR;
				yRanges.ranges.push_back(yNearRange);
			}
		}

		RangeType yMaxExteriorRange;
		yMaxExteriorRange.range.first = yNearRanges[yNearRanges.size() - 1].range.second;
		yMaxExteriorRange.range.second = yMax;
		yMaxExteriorRange.type = RangeType::MAX_EXTERIOR;
		yRanges.ranges.push_back(yMaxExteriorRange);

	}

	std::vector<double> xPoints;
	std::vector<double> yPoints;
	std::vector<double> zPoints;

	std::vector<double> xSurfaces;
	std::vector<double> ySurfaces;
	std::vector<double> zSurfaces;

	// Create points for mesh data
	for (size_t s = 0; s < surfaces.size(); s++)
	{
		for (std::size_t v = 0; v < surfaces[s].polygon.outer().size(); v++)
		{
			xPoints.push_back(surfaces[s].polygon.outer()[v].get<0>());
			yPoints.push_back(surfaces[s].polygon.outer()[v].get<1>());
		}
		zPoints.push_back(surfaces[s].zMax);
		zPoints.push_back(surfaces[s].zMin);

		if (surfaces[s].orientation == Surface::X_POS ||
			surfaces[s].orientation == Surface::X_NEG)
		{
			xSurfaces.push_back(surfaces[s].polygon.outer()[0].get<0>());
		}
		else if (surfaces[s].orientation == Surface::Y_POS ||
				 surfaces[s].orientation == Surface::Y_NEG)
		{
			ySurfaces.push_back(surfaces[s].polygon.outer()[0].get<1>());
		}
		else // if (surfaces[s].orientation == Surface::Z_POS ||
			 //     surfaces[s].orientation == Surface::Z_NEG)
		{
			zSurfaces.push_back(surfaces[s].zMin);
		}

	}

	for (size_t b = 0; b < blocks.size(); b++)
	{
		for (std::size_t v = 0; v < surfaces[b].polygon.outer().size(); v++)
		{
			xPoints.push_back(blocks[b].polygon.outer()[v].get<0>());
			yPoints.push_back(blocks[b].polygon.outer()[v].get<1>());
		}
		zPoints.push_back(blocks[b].zMax);
		zPoints.push_back(blocks[b].zMin);
	}

	// Sort the points vectors
	sort(xPoints.begin(), xPoints.end());
	sort(yPoints.begin(), yPoints.end());
	sort(zPoints.begin(), zPoints.end());

	// erase duplicate elements
	for (size_t i = 1; i < xPoints.size(); i++)
	{
		if (isEqual(xPoints[i], xPoints[i-1]))
		{
			xPoints.erase(xPoints.begin() + i);
			i -= 1;
		}
	}

	for (size_t j = 1; j < yPoints.size(); j++)
	{
		if (isEqual(yPoints[j], yPoints[j-1]))
		{
			yPoints.erase(yPoints.begin() + j);
			j -= 1;
		}
	}

	for (size_t k = 1; k < zPoints.size(); k++)
	{
		if (isEqual(zPoints[k], zPoints[k-1]))
		{
			zPoints.erase(zPoints.begin() + k);
			k -= 1;
		}
	}

	// Sort the surfaces vectors
	sort(xSurfaces.begin(), xSurfaces.end());
	sort(ySurfaces.begin(), ySurfaces.end());
	sort(zSurfaces.begin(), zSurfaces.end());

	// erase (approximately) duplicate elements
	for (size_t i = 1; i < xSurfaces.size(); i++)
	{
		if (isEqual(xSurfaces[i], xSurfaces[i-1]))
		{
			xSurfaces.erase(xSurfaces.begin() + i);
			i -= 1;
		}
	}

	for (size_t j = 1; j < ySurfaces.size(); j++)
	{
		if (isEqual(ySurfaces[j], ySurfaces[j-1]))
		{
			ySurfaces.erase(ySurfaces.begin() + j);
			j -= 1;
		}
	}

	for (size_t k = 1; k < zSurfaces.size(); k++)
	{
		if (isEqual(zSurfaces[k], zSurfaces[k-1]))
		{
			zSurfaces.erase(zSurfaces.begin() + k);
			k -= 1;
		}
	}

	// re-add the extra surface elements to create zero-thickness cells
	for (size_t i = 0; i < xSurfaces.size(); i++)
	{
		xPoints.push_back(xSurfaces[i]);
	}

	for (size_t j = 0; j < ySurfaces.size(); j++)
	{
		yPoints.push_back(ySurfaces[j]);
	}

	for (size_t k = 0; k < zSurfaces.size(); k++)
	{
		zPoints.push_back(zSurfaces[k]);
	}

	// Sort the range vectors again this time it includes doubles for zero thickness cells
	sort(xPoints.begin(), xPoints.end());
	sort(yPoints.begin(), yPoints.end());
	sort(zPoints.begin(), zPoints.end());


	xMeshData.points = xPoints;
	yMeshData.points = yPoints;
	zMeshData.points = zPoints;

	std::vector<Interval> xIntervals;
	std::vector<Interval> yIntervals;
	std::vector<Interval> zIntervals;

	for (size_t i = 1; i < xPoints.size(); i++)
	{
		if (isEqual(xPoints[i], xPoints[i-1]))
			xIntervals.push_back(zeroThickness);
		else if (xRanges.isType(xPoints[i],RangeType::MIN_INTERIOR))
			xIntervals.push_back(minInterior);
		else if (xRanges.isType(xPoints[i],RangeType::MID_INTERIOR))
			xIntervals.push_back(midInterior);
		else if (xRanges.isType(xPoints[i],RangeType::NEAR))
			xIntervals.push_back(near);
		else if (xRanges.isType(xPoints[i],RangeType::MIN_EXTERIOR))
			xIntervals.push_back(minExterior);
		else if (xRanges.isType(xPoints[i],RangeType::MAX_EXTERIOR))
			xIntervals.push_back(maxExterior);
	}

	for (size_t j = 1; j < yPoints.size(); j++)
	{
		if (isEqual(yPoints[j], yPoints[j-1]))
			yIntervals.push_back(zeroThickness);
		else if (yRanges.isType(yPoints[j],RangeType::MIN_INTERIOR))
			yIntervals.push_back(minInterior);
		else if (yRanges.isType(yPoints[j],RangeType::MID_INTERIOR))
			yIntervals.push_back(midInterior);
		else if (yRanges.isType(yPoints[j],RangeType::NEAR))
			yIntervals.push_back(near);
		else if (yRanges.isType(yPoints[j],RangeType::MIN_EXTERIOR))
			yIntervals.push_back(minExterior);
		else if (yRanges.isType(yPoints[j],RangeType::MAX_EXTERIOR))
			yIntervals.push_back(maxExterior);
	}

	for (size_t k = 1; k < zPoints.size(); k++)
	{
		if (isEqual(zPoints[k], zPoints[k-1]))
			zIntervals.push_back(zeroThickness);
		else if (zRanges.isType(zPoints[k],RangeType::DEEP))
			zIntervals.push_back(deep);
		else if (zRanges.isType(zPoints[k],RangeType::NEAR))
			zIntervals.push_back(near);
	}

	xMeshData.intervals = xIntervals;
	yMeshData.intervals = yIntervals;
	zMeshData.intervals = zIntervals;

}

#endif /* INPUT_CPP_ */



