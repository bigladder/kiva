/* Copyright (c) 2012-2016 Big Ladder Software. All rights reserved.
* See the LICENSE file for additional terms and conditions. */

#include "Foundation.hpp"

namespace Kiva {

static const double PI = 4.0*atan(1.0);

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

inline bool compareRanges(RangeType first,  RangeType second)
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

void Foundation::createMeshData()
{
  std::size_t nV = polygon.outer().size();

  for (std::size_t v = 0; v < nV; v++)
  {
    double thisX = polygon.outer()[v].get<0>();
    double thisY = polygon.outer()[v].get<1>();
    double nextX, nextY;

    if (v < nV -1)
    {
      nextX = polygon.outer()[v+1].get<0>();
      nextY = polygon.outer()[v+1].get<1>();
    }
    else
    {
      nextX = polygon.outer()[0].get<0>();
      nextY = polygon.outer()[0].get<1>();
    }

    Polygon3 poly;
    poly.outer().push_back(Point3(thisX,thisY,0.0));
    poly.outer().push_back(Point3(thisX,thisY,buildingHeight));
    poly.outer().push_back(Point3(nextX,nextY,buildingHeight));
    poly.outer().push_back(Point3(nextX,nextY,0.0));
    buildingSurfaces.push_back(poly);
  }

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
  near.maxGrowthCoeff = mesh.maxNearGrowthCoeff;
  near.minCellDim = mesh.minCellDim;
  near.growthDir = Interval::CENTERED;

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
  midInterior.growthDir = Interval::CENTERED;

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
  double zSlab = zMax - foundationDepth;

  double zNearDeep = std::min(zSlab, zGrade);  // This will change depending on configuration

  // Set misc. "X/Y" dimensions (relative to foundation outline --
  // currently the interior of the wall)

  double xyWallExterior;
  if (hasWall)
  {
    if (hasExteriorVerticalInsulation)
      xyWallExterior = wall.totalWidth()
          + exteriorVerticalInsulation.layer.thickness;

    else
      xyWallExterior = wall.totalWidth();
  }
  else
    xyWallExterior = 0.0;

  double xyWallInterior;
  if (hasInteriorVerticalInsulation)
    xyWallInterior = -interiorVerticalInsulation.layer.thickness;
  else
    xyWallInterior = 0.0;

  double xySlabPerimeter = xyWallInterior;

  double xyNearInt = xyWallInterior;  // This will change depending on configuration
  double xyNearExt = xyWallExterior;  // This will change depending on configuration

  double xyIntHIns = 0.0;
  double zIntHIns = 0.0;
  if (hasInteriorHorizontalInsulation)
  {
    xyIntHIns = -interiorHorizontalInsulation.width;
    zIntHIns = zMax - interiorHorizontalInsulation.depth - interiorHorizontalInsulation.layer.thickness;

    if(zIntHIns < zNearDeep)
      zNearDeep = zIntHIns;
    if(xyIntHIns < xyNearInt)
      xyNearInt = xyIntHIns;
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

    if (zIntVIns < zNearDeep)
      zNearDeep = zIntVIns;

    if (isGreaterThan(zIntVIns,zSlab))
      xySlabPerimeter = 0.0;
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
    xyPerimeterSurface = xySlabPerimeter;

  double zWall = 0.0;
  if (hasWall)
  {
    zWall = zSlabBottom - wall.footerDepth;
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
    xyExtHIns = wall.totalWidth() + exteriorHorizontalInsulation.width;
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

  double xMin, xMax, yMin, yMax;

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


  // Set 3D foundation areas (for calculation of total heat transfer rates)
  surfaceAreas[Surface::ST_SLAB_CORE] = boost::geometry::area(offset(polygon, xyPerimeterSurface));
  surfaceAreas[Surface::ST_SLAB_PERIM] = boost::geometry::area(offset(polygon, xyPerimeterSurface)) - surfaceAreas[Surface::ST_SLAB_CORE];

  if (hasInteriorVerticalInsulation && isGreaterThan(zIntVIns, zSlab)) {
    surfaceAreas[Surface::ST_WALL_INT] = boost::geometry::perimeter(offset(polygon, xyWallInterior))*(zMax - zIntVIns);
    surfaceAreas[Surface::ST_WALL_INT] += boost::geometry::area(polygon) - boost::geometry::area(offset(polygon, xyWallInterior));
    surfaceAreas[Surface::ST_WALL_INT] += boost::geometry::perimeter(polygon)*(zIntVIns - zSlab);
  }
  else {
    surfaceAreas[Surface::ST_WALL_INT] = boost::geometry::perimeter(offset(polygon, xyWallInterior))*foundationDepth;
  }

  hasSurface[Surface::ST_SLAB_CORE] = true;
  hasSurface[Surface::ST_SLAB_PERIM] = hasPerimeterSurface;
  hasSurface[Surface::ST_WALL_INT] = (foundationDepth > 0);

  if (numberOfDimensions == 2 )
  {

    // TODO: 2D
    double area = boost::geometry::area(polygon);  // [m2] Area of foundation
    double perimeter = boost::geometry::perimeter(polygon);  // [m] Perimeter of foundation

    linearAreaMultiplier = 1.0;

    double ap = area/perimeter;

    if (reductionStrategy == RS_AP)
    {
      twoParameters = false;
      if (coordinateSystem == CS_CYLINDRICAL)
      {
        reductionLength2 = 2.0*ap;
      }
      else if (coordinateSystem == CS_CARTESIAN)
      {
        reductionLength2 = ap;
      }
    }

    if (reductionStrategy == RS_RR)
    {
      twoParameters = false;
      double rrA = (perimeter - sqrt(perimeter*perimeter - 4*PI*area))/PI;
      double rrB = (perimeter - PI*rrA)*0.5;
      reductionLength2 = (rrA)*0.5;
      linearAreaMultiplier = rrB;
    }

    xMin = 0.0;
    xMax = reductionLength2 + farFieldWidth;

    yMin = 0.0;
    yMax = 1.0;

    double xRef2 = reductionLength2;
    double xRef1 = reductionLength1;

    // Symmetry Surface
    {
      Surface surface;
      surface.type = Surface::ST_SYMMETRY;
      surface.xMin = xMin;
      surface.xMax = xMin;
      surface.yMin = 0.0;
      surface.yMax = 1.0;
      surface.setSquarePolygon();
      surface.zMin = zMin;
      if (twoParameters)
        surface.zMax = zGrade;
      else
        surface.zMax = zSlab;
      surface.boundaryConditionType = Surface::ZERO_FLUX;
      surface.orientation = Surface::X_NEG;
      surfaces.push_back(surface);
    }

    if(isGreaterThan(foundationDepth,0.0))
    {
      // Interior Wall Surface

      if (hasInteriorVerticalInsulation && isGreaterThan(zIntVIns,zSlab))
      {
        {
          Surface surface;
          surface.type = Surface::ST_WALL_INT;
          surface.xMin = xRef2 + xyWallInterior;
          surface.xMax = xRef2 + xyWallInterior;
          surface.yMin = 0.0;
          surface.yMax = 1.0;
          surface.setSquarePolygon();
          surface.zMin = zIntVIns;
          surface.zMax = zMax;
          surface.boundaryConditionType = Surface::INTERIOR_FLUX;
          surface.orientation = Surface::X_NEG;
          surface.emissivity = wall.interiorEmissivity;
          surfaces.push_back(surface);
        }
        {
          Surface surface;
          surface.type = Surface::ST_WALL_INT;
          surface.xMin = xRef2 + xyWallInterior;
          surface.xMax = xRef2;
          surface.yMin = 0.0;
          surface.yMax = 1.0;
          surface.setSquarePolygon();
          surface.zMin = zIntVIns;
          surface.zMax = zIntVIns;
          surface.boundaryConditionType = Surface::INTERIOR_FLUX;
          surface.orientation = Surface::Z_NEG;
          surface.emissivity = wall.interiorEmissivity;
          surfaces.push_back(surface);
        }
        {
          Surface surface;
          surface.type = Surface::ST_WALL_INT;
          surface.xMin = xRef2;
          surface.xMax = xRef2;
          surface.yMin = 0.0;
          surface.yMax = 1.0;
          surface.setSquarePolygon();
          surface.zMin = zSlab;
          surface.zMax = zIntVIns;
          surface.boundaryConditionType = Surface::INTERIOR_FLUX;
          surface.orientation = Surface::X_NEG;
          surface.emissivity = wall.interiorEmissivity;
          surfaces.push_back(surface);
        }
        if (twoParameters)
        {
          {
            Surface surface;
            surface.type = Surface::ST_WALL_INT;
            surface.xMin = xRef1 - xyWallInterior;
            surface.xMax = xRef1 - xyWallInterior;
            surface.yMin = 0.0;
            surface.yMax = 1.0;
            surface.setSquarePolygon();
            surface.zMin = zIntVIns;
            surface.zMax = zMax;
            surface.boundaryConditionType = Surface::INTERIOR_FLUX;
            surface.orientation = Surface::X_POS;
            surface.emissivity = wall.interiorEmissivity;
            surfaces.push_back(surface);
          }
          {
            Surface surface;
            surface.type = Surface::ST_WALL_INT;
            surface.xMin = xRef1;
            surface.xMax = xRef1 - xyWallInterior;
            surface.yMin = 0.0;
            surface.yMax = 1.0;
            surface.setSquarePolygon();
            surface.zMin = zIntVIns;
            surface.zMax = zIntVIns;
            surface.boundaryConditionType = Surface::INTERIOR_FLUX;
            surface.orientation = Surface::Z_NEG;
            surface.emissivity = wall.interiorEmissivity;
            surfaces.push_back(surface);
          }
          {
            Surface surface;
            surface.type = Surface::ST_WALL_INT;
            surface.xMin = xRef1;
            surface.xMax = xRef1;
            surface.yMin = 0.0;
            surface.yMax = 1.0;
            surface.setSquarePolygon();
            surface.zMin = zSlab;
            surface.zMax = zIntVIns;
            surface.boundaryConditionType = Surface::INTERIOR_FLUX;
            surface.orientation = Surface::X_POS;
            surface.emissivity = wall.interiorEmissivity;
            surfaces.push_back(surface);
          }
        }

      }
      else
      {
        Surface surface;
        surface.type = Surface::ST_WALL_INT;
        surface.xMin = xRef2 + xyWallInterior;
        surface.xMax = xRef2 + xyWallInterior;
        surface.yMin = 0.0;
        surface.yMax = 1.0;
        surface.setSquarePolygon();
        surface.zMin = zSlab;
        surface.zMax = zMax;
        surface.boundaryConditionType = Surface::INTERIOR_FLUX;
        surface.orientation = Surface::X_NEG;
        surface.emissivity = wall.interiorEmissivity;
        surfaces.push_back(surface);

        if (twoParameters)
        {
          Surface surface;
          surface.type = Surface::ST_WALL_INT;
          surface.xMin = xRef1 - xyWallInterior;
          surface.xMax = xRef1 - xyWallInterior;
          surface.yMin = 0.0;
          surface.yMax = 1.0;
          surface.setSquarePolygon();
          surface.zMin = zSlab;
          surface.zMax = zMax;
          surface.boundaryConditionType = Surface::INTERIOR_FLUX;
          surface.orientation = Surface::X_POS;
          surface.emissivity = wall.interiorEmissivity;
          surfaces.push_back(surface);
        }
      }

      // Interior Air Symmetry Temperature
      if (!twoParameters)
      {
        Surface surface;
        surface.type = Surface::ST_SYMMETRY_AIR;
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

    if(isGreaterThan(zMax,0.0))
    {
      // Exterior Wall Surface
      {
        Surface surface;
        surface.type = Surface::ST_WALL_EXT;
        surface.xMin = xRef2 + xyWallExterior;
        surface.xMax = xRef2 + xyWallExterior;
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

      if(twoParameters)
      {
        Surface surface;
        surface.type = Surface::ST_WALL_EXT;
        surface.xMin = xRef1 - xyWallExterior;
        surface.xMax = xRef1 - xyWallExterior;
        surface.yMin = 0.0;
        surface.yMax = 1.0;
        surface.setSquarePolygon();
        surface.zMin = zGrade;
        surface.zMax = zMax;
        surface.boundaryConditionType = Surface::EXTERIOR_FLUX;
        surface.orientation = Surface::X_NEG;
        surface.emissivity = wall.exteriorEmissivity;
        surface.absorptivity = wall.exteriorAbsorptivity;
        surfaces.push_back(surface);
      }

      // Exterior Air Right Surface
      {
        Surface surface;
        surface.type = Surface::ST_FAR_FIELD_AIR;
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
      surface.type = Surface::ST_FAR_FIELD;
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
      surface.type = Surface::ST_DEEP_GROUND;
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
      surface.type = Surface::ST_DEEP_GROUND;
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
      surface.type = Surface::ST_SLAB_CORE;
      if (!twoParameters)
      {
        surface.xMin = xMin;
        surface.xMax = xRef2 + xyPerimeterSurface;
      }
      else
      {
        surface.xMin = xRef1 - xyPerimeterSurface;
        surface.xMax = xRef2 + xyPerimeterSurface;
      }
      surface.yMin = 0.0;
      surface.yMax = 1.0;
      surface.setSquarePolygon();
      surface.zMin = zSlab;
      surface.zMax = zSlab;
      surface.boundaryConditionType = Surface::INTERIOR_FLUX;
      surface.orientation = Surface::Z_POS;
      surface.emissivity = slab.emissivity;
      surfaces.push_back(surface);
    }
    if (hasPerimeterSurface)
    {
      {
        Surface surface;
        surface.type = Surface::ST_SLAB_PERIM;
        surface.xMin = xRef2 + xyPerimeterSurface;
        surface.xMax = xRef2 + xySlabPerimeter;
        surface.yMin = 0.0;
        surface.yMax = 1.0;
        surface.setSquarePolygon();
        surface.zMin = zSlab;
        surface.zMax = zSlab;
        surface.boundaryConditionType = Surface::INTERIOR_FLUX;
        surface.orientation = Surface::Z_POS;
        surface.emissivity = slab.emissivity;
        surfaces.push_back(surface);
      }
      if (twoParameters)
      {
        Surface surface;
        surface.type = Surface::ST_SLAB_PERIM;
        surface.xMin = xRef1 - xySlabPerimeter;
        surface.xMax = xRef1 - xyPerimeterSurface;
        surface.yMin = 0.0;
        surface.yMax = 1.0;
        surface.setSquarePolygon();
        surface.zMin = zSlab;
        surface.zMax = zSlab;
        surface.boundaryConditionType = Surface::INTERIOR_FLUX;
        surface.orientation = Surface::Z_POS;
        surface.emissivity = slab.emissivity;
        surfaces.push_back(surface);
      }
    }

    // Grade
    {
      Surface surface;
      surface.type = Surface::ST_GRADE;
      surface.xMin = xRef2 + xyWallExterior;
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
    if (twoParameters)
    {
      Surface surface;
      surface.type = Surface::ST_GRADE;
      surface.xMin = xMin;
      surface.xMax = xRef1 - xyWallExterior;
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
    if(foundationDepth > 0.0)
    {
      // Interior Air Top Surface
      Surface surface;
      surface.type = Surface::ST_TOP_AIR_INT;
      if (!twoParameters)
      {
        surface.xMin = xMin;
        surface.xMax = xRef2 + xyWallInterior;
      }
      else
      {
        surface.xMin = xRef1 - xyWallInterior;
        surface.xMax = xRef2 + xyWallInterior;
      }
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
      {
        // Exterior Air Top Surface
        Surface surface;
        surface.type = Surface::ST_TOP_AIR_EXT;
        surface.xMin = xRef2 + xyWallExterior;
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
      if (twoParameters)
      {
        // Exterior Air Top Surface
        Surface surface;
        surface.type = Surface::ST_TOP_AIR_EXT;
        surface.xMin = xMin;
        surface.xMax = xRef1 - xyWallExterior;
        surface.yMin = 0.0;
        surface.yMax = 1.0;
        surface.setSquarePolygon();
        surface.zMin = zMax;
        surface.zMax = zMax;
        surface.boundaryConditionType = Surface::EXTERIOR_TEMPERATURE;
        surface.orientation = Surface::Z_POS;
        surfaces.push_back(surface);
      }

    }

    if (wallTopBoundary == WTB_LINEAR_DT)
    {
      // Wall Top
      if(hasWall)
      {
        double position = 0.0;
        double& Tin = wallTopInteriorTemperature;
        double& Tout = wallTopExteriorTemperature;
        double Tdiff = (Tin - Tout);
        std::size_t N = xyWallExterior/mesh.minCellDim;
        double temperature = Tin - (1.0/N)/2*Tdiff;

        for (std::size_t n = 1; n <= N; n++)
        {
          Surface surface;
          surface.type = Surface::ST_WALL_TOP;
          surface.xMin = xRef2 + position;
          surface.xMax = xRef2 + position + xyWallExterior/N;
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
          temperature -= (1.0/N)*Tdiff;

        }

        if (twoParameters)
        {
          double& Tin = wallTopInteriorTemperature;
          double& Tout = wallTopExteriorTemperature;
          double Tdiff = (Tin - Tout);
          std::size_t N = xyWallExterior/mesh.minCellDim;
          double temperature = Tin - (1.0/N)/2*Tdiff;

          for (std::size_t n = 1; n <= N; n++)
          {
            Surface surface;
            surface.type = Surface::ST_WALL_TOP;
            surface.xMin = xRef1 - position - xyWallExterior/N;
            surface.xMax = xRef1 - position;
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
            temperature -= (1.0/N)*Tdiff;

          }

        }
      }
    }
    else
    {
      // Wall Top
      if(hasWall)
      {
        {
          Surface surface;
          surface.type = Surface::ST_WALL_TOP;
          surface.xMin = xRef2 + xyWallInterior;
          surface.xMax = xRef2 + xyWallExterior;
          surface.yMin = 0.0;
          surface.yMax = 1.0;
          surface.setSquarePolygon();
          surface.zMin = zMax;
          surface.zMax = zMax;
          surface.boundaryConditionType = Surface::ZERO_FLUX;
          surface.orientation = Surface::Z_POS;
          surfaces.push_back(surface);
        }
        if (twoParameters)
        {
          Surface surface;
          surface.type = Surface::ST_WALL_TOP;
          surface.xMin = xRef1 - xyWallExterior;
          surface.xMax = xRef1 - xyWallInterior;
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
    }

    // Interior Horizontal Insulation
    if (hasInteriorHorizontalInsulation)
    {
      {
        Block block;
        block.material = interiorHorizontalInsulation.layer.material;
        block.blockType = Block::SOLID;
        block.xMin = xRef2 + xyIntHIns;
        block.xMax = xRef2;
        block.yMin = 0.0;
        block.yMax = 1.0;
        block.setSquarePolygon();
        block.zMin = zIntHIns;
        block.zMax = zIntHIns + interiorHorizontalInsulation.layer.thickness;
        blocks.push_back(block);
      }
      if (twoParameters)
      {
        Block block;
        block.material = interiorHorizontalInsulation.layer.material;
        block.blockType = Block::SOLID;
        block.xMin = xRef1;
        block.xMax = xRef1 - xyIntHIns;
        block.yMin = 0.0;
        block.yMax = 1.0;
        block.setSquarePolygon();
        block.zMin = zIntHIns;
        block.zMax = zIntHIns + interiorHorizontalInsulation.layer.thickness;
        blocks.push_back(block);
      }
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
        if (!twoParameters)
        {
          block.xMin = xMin;
          block.xMax = xRef2;
        }
        else
        {
          block.xMin = xRef1;
          block.xMax = xRef2;
        }
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
      {
        Block block;
        block.material = interiorVerticalInsulation.layer.material;
        block.blockType = Block::SOLID;
        block.xMin = xRef2 + xyWallInterior;
        block.xMax = xRef2;
        block.yMin = 0.0;
        block.yMax = 1.0;
        block.setSquarePolygon();
        block.zMin = zIntVIns;
        block.zMax = zMax;
        blocks.push_back(block);
      }

      if (isGreaterThan(zIntVIns,zSlab))
      {
        Block block;
        block.material = air;
        block.blockType = Block::INTERIOR_AIR;
        block.xMin = xRef2 + xyWallInterior;
        block.xMax = xRef2;
        block.yMin = 0.0;
        block.yMax = 1.0;
        block.setSquarePolygon();
        block.zMin = zSlab;
        block.zMax = zIntVIns;
        blocks.push_back(block);
      }

      if (twoParameters)
      {
        {
          Block block;
          block.material = interiorVerticalInsulation.layer.material;
          block.blockType = Block::SOLID;
          block.xMin = xRef1;
          block.xMax = xRef1 - xyWallInterior;
          block.yMin = 0.0;
          block.yMax = 1.0;
          block.setSquarePolygon();
          block.zMin = zIntVIns;
          block.zMax = zMax;
          blocks.push_back(block);
        }

        if (isGreaterThan(zIntVIns,zSlab))
        {
          Block block;
          block.material = air;
          block.blockType = Block::INTERIOR_AIR;
          block.xMin = xRef1;
          block.xMax = xRef1 - xyWallInterior;
          block.yMin = 0.0;
          block.yMax = 1.0;
          block.setSquarePolygon();
          block.zMin = zSlab;
          block.zMax = zIntVIns;
          blocks.push_back(block);
        }
      }
    }

    // Indoor Air
    {
      Block block;
      block.material = air;
      block.blockType = Block::INTERIOR_AIR;
      if (twoParameters)
      {
        block.xMin = xRef1 - xyWallInterior;
        block.xMax = xRef2 + xyWallInterior;
      }
      else
      {
        block.xMin = xMin;
        block.xMax = xRef2 + xyWallInterior;
      }
      block.yMin = 0.0;
      block.yMax = 1.0;
      block.setSquarePolygon();
      block.zMin = zSlab;
      block.zMax = zMax;
      blocks.push_back(block);
    }

    if (hasWall)
    {
      {
        double xPosition = xRef2;

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

      if (twoParameters)
      {
        double xPosition = xRef1;

        // Foundation Wall
        for (int n = wall.layers.size() - 1; n >= 0; n--)
        {
          Block block;
          block.material = wall.layers[n].material;
          block.blockType = Block::SOLID;
          block.xMin = xPosition - wall.layers[n].thickness;
          block.xMax = xPosition;
          block.yMin = 0.0;
          block.yMax = 1.0;
          block.setSquarePolygon();
          block.zMin = zWall;
          block.zMax = zMax;
          xPosition = block.xMin;
          blocks.push_back(block);
        }
      }
    }

    // Exterior Vertical Insulation
    if (hasExteriorVerticalInsulation)
    {
      {
        Block block;
        block.material = exteriorVerticalInsulation.layer.material;
        block.blockType = Block::SOLID;
        block.xMin = xRef2 + wall.totalWidth();
        block.xMax = xRef2 + xyWallExterior;
        block.yMin = 0.0;
        block.yMax = 1.0;
        block.setSquarePolygon();
        block.zMin = zExtVIns;
        block.zMax = zMax;
        blocks.push_back(block);
      }
      if (twoParameters)
      {
        Block block;
        block.material = exteriorVerticalInsulation.layer.material;
        block.blockType = Block::SOLID;
        block.xMin = xRef1 - xyWallExterior;
        block.xMax = xRef1 - wall.totalWidth();
        block.yMin = 0.0;
        block.yMax = 1.0;
        block.setSquarePolygon();
        block.zMin = zExtVIns;
        block.zMax = zMax;
        blocks.push_back(block);
      }
    }

    // Exterior Horizontal Insulation
    if (hasExteriorHorizontalInsulation)
    {
      {
        Block block;
        block.material = exteriorHorizontalInsulation.layer.material;
        block.blockType = Block::SOLID;
        block.xMin = xRef2 + wall.totalWidth();
        block.xMax = xRef2 + xyExtHIns;
        block.yMin = 0.0;
        block.yMax = 1.0;
        block.setSquarePolygon();
        block.zMin = zExtHIns;
        block.zMax = zExtHIns + exteriorHorizontalInsulation.layer.thickness;
        blocks.push_back(block);
      }
      if (twoParameters)
      {
        Block block;
        block.material = exteriorHorizontalInsulation.layer.material;
        block.blockType = Block::SOLID;
        block.xMin = xRef1 - xyExtHIns;
        block.xMax = xRef1 - wall.totalWidth();
        block.yMin = 0.0;
        block.yMax = 1.0;
        block.setSquarePolygon();
        block.zMin = zExtHIns;
        block.zMax = zExtHIns + exteriorHorizontalInsulation.layer.thickness;
        blocks.push_back(block);
      }
    }

    // Exterior Air
    {
      Block block;
      block.material = air;
      block.blockType = Block::EXTERIOR_AIR;
      block.xMin = xRef2 + xyWallExterior;
      block.xMax = xMax;
      block.yMin = 0.0;
      block.yMax = 1.0;
      block.setSquarePolygon();
      block.zMin = zGrade;
      block.zMax = zMax;
      blocks.push_back(block);
    }
    if (twoParameters)
    {
      Block block;
      block.material = air;
      block.blockType = Block::EXTERIOR_AIR;
      block.xMin = xMin;
      block.xMax = xRef1 - xyWallExterior;
      block.yMin = 0.0;
      block.yMax = 1.0;
      block.setSquarePolygon();
      block.zMin = zGrade;
      block.zMax = zMax;
      blocks.push_back(block);
    }

    // Set range types

    if (!twoParameters)
    {
      RangeType xInteriorRange;
      xInteriorRange.range.first = xMin;
      xInteriorRange.range.second = xRef2 + xyNearInt;
      xInteriorRange.type = RangeType::MIN_INTERIOR;
      xRanges.ranges.push_back(xInteriorRange);

      RangeType xNearRange;
      xNearRange.range.first = xRef2 + xyNearInt;
      xNearRange.range.second = xRef2 + xyNearExt;
      xNearRange.type = RangeType::NEAR;
      xRanges.ranges.push_back(xNearRange);

      RangeType xExteriorRange;
      xExteriorRange.range.first = xRef2 + xyNearExt;
      xExteriorRange.range.second = xMax;
      xExteriorRange.type = RangeType::MAX_EXTERIOR;
      xRanges.ranges.push_back(xExteriorRange);
    }
    else
    {
      RangeType xMinExteriorRange;
      xMinExteriorRange.range.first = xMin;
      xMinExteriorRange.range.second = xRef1 - xyNearExt;
      xMinExteriorRange.type = RangeType::MIN_INTERIOR;
      xRanges.ranges.push_back(xMinExteriorRange);

      RangeType xNearRange1;
      xNearRange1.range.first = xRef1 - xyNearExt;
      xNearRange1.range.second = xRef1 - xyNearInt;
      xNearRange1.type = RangeType::NEAR;
      xRanges.ranges.push_back(xNearRange1);

      RangeType xInteriorRange;
      xInteriorRange.range.first = xRef1 - xyNearInt;
      xInteriorRange.range.second = xRef2 + xyNearInt;
      xInteriorRange.type = RangeType::MID_INTERIOR;
      xRanges.ranges.push_back(xInteriorRange);

      RangeType xNearRange2;
      xNearRange2.range.first = xRef2 + xyNearInt;
      xNearRange2.range.second = xRef2 + xyNearExt;
      xNearRange2.type = RangeType::NEAR;
      xRanges.ranges.push_back(xNearRange2);

      RangeType xMaxExteriorRange;
      xMaxExteriorRange.range.first = xRef2 + xyNearExt;
      xMaxExteriorRange.range.second = xMax;
      xMaxExteriorRange.type = RangeType::MAX_EXTERIOR;
      xRanges.ranges.push_back(xMaxExteriorRange);
    }

  }
  else if(numberOfDimensions == 3 && !useSymmetry)
  {
    // TODO 3D
    Box boundingBox;
    boost::geometry::envelope(polygon, boundingBox);

    double xMinBB = boundingBox.min_corner().get<0>();
    double yMinBB = boundingBox.min_corner().get<1>();

    double xMaxBB = boundingBox.max_corner().get<0>();
    double yMaxBB = boundingBox.max_corner().get<1>();

    xMin = xMinBB - farFieldWidth;
    yMin = yMinBB - farFieldWidth;

    xMax = xMaxBB + farFieldWidth;
    yMax = yMaxBB + farFieldWidth;

    if(isGreaterThan(foundationDepth,0.0))
    {
      // Interior Wall Surface
      if (hasInteriorVerticalInsulation && isGreaterThan(zIntVIns,zSlab))
      {
        {
          Polygon poly;
          poly = offset(polygon, xyWallInterior);

          for (std::size_t v = 0; v < nV; v++)
          {
            Surface surface;
            surface.type = Surface::ST_WALL_INT;
            surface.xMin = getXmin(poly,v);
            surface.xMax = getXmax(poly,v);
            surface.yMin = getYmin(poly,v);
            surface.yMax = getYmax(poly,v);
            surface.setSquarePolygon();
            surface.zMin = zIntVIns;
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
        {
          Polygon poly;
          poly = polygon;

          Polygon temp;
          temp = offset(polygon, xyWallInterior);
          Ring ring;
          boost::geometry::convert(temp, ring);
          boost::geometry::reverse(ring);

          poly.inners().push_back(ring);

          Surface surface;
          surface.type = Surface::ST_WALL_INT;
          surface.polygon = poly;
          surface.zMin = zIntVIns;
          surface.zMax = zIntVIns;
          surface.boundaryConditionType = Surface::INTERIOR_FLUX;
          surface.orientation = Surface::Z_NEG;
          surface.emissivity = wall.interiorEmissivity;
          surfaces.push_back(surface);
        }
        {
          Polygon poly;
          poly = polygon;

          for (std::size_t v = 0; v < nV; v++)
          {
            Surface surface;
            surface.type = Surface::ST_WALL_INT;
            surface.xMin = getXmin(poly,v);
            surface.xMax = getXmax(poly,v);
            surface.yMin = getYmin(poly,v);
            surface.yMax = getYmax(poly,v);
            surface.setSquarePolygon();
            surface.zMin = zSlab;
            surface.zMax = zIntVIns;
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
      else
      {
        Polygon poly;
        poly = offset(polygon, xyWallInterior);

        for (std::size_t v = 0; v < nV; v++)
        {
          Surface surface;
          surface.type = Surface::ST_WALL_INT;
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

    if(isGreaterThan(zMax,0.0))
    {
      // Exterior Wall Surface
      {
        Polygon poly;
        poly = offset(polygon, xyWallExterior);

        for (std::size_t v = 0; v < nV; v++)
        {
          Surface surface;
          surface.type = Surface::ST_WALL_EXT;
          surface.xMin = getXmin(poly,v);
          surface.xMax = getXmax(poly,v);
          surface.yMin = getYmin(poly,v);
          surface.yMax = getYmax(poly,v);
          surface.zMin = zGrade;
          surface.zMax = zMax;
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

      // Exterior Air Perimeter Surfaces
      {
        {
          // X Min
          Surface surface;
          surface.type = Surface::ST_FAR_FIELD_AIR;
          surface.xMin =  xMin;
          surface.xMax =  xMin;
          surface.yMin =  yMin;
          surface.yMax =  yMax;
          surface.setSquarePolygon();
          surface.zMin =  zGrade;
          surface.zMax =  zMax;
          surface.boundaryConditionType = Surface::EXTERIOR_TEMPERATURE;
          surface.orientation =  Surface::X_NEG;
          surfaces.push_back(surface);
        }

        {
          // X Max
          Surface surface;
          surface.type = Surface::ST_FAR_FIELD_AIR;
          surface.xMin =  xMax;
          surface.xMax =  xMax;
          surface.yMin =  yMin;
          surface.yMax =  yMax;
          surface.setSquarePolygon();
          surface.zMin =  zGrade;
          surface.zMax =  zMax;
          surface.boundaryConditionType = Surface::EXTERIOR_TEMPERATURE;
          surface.orientation =  Surface::X_POS;
          surfaces.push_back(surface);
        }

        {
          // Y Min
          Surface surface;
          surface.type = Surface::ST_FAR_FIELD_AIR;
          surface.xMin =  xMin;
          surface.xMax =  xMax;
          surface.yMin =  yMin;
          surface.yMax =  yMin;
          surface.setSquarePolygon();
          surface.zMin =  zGrade;
          surface.zMax =  zMax;
          surface.boundaryConditionType = Surface::EXTERIOR_TEMPERATURE;
          surface.orientation =  Surface::Y_NEG;
          surfaces.push_back(surface);
        }

        {
          // Y Max
          Surface surface;
          surface.type = Surface::ST_FAR_FIELD_AIR;
          surface.xMin =  xMin;
          surface.xMax =  xMax;
          surface.yMin =  yMax;
          surface.yMax =  yMax;
          surface.setSquarePolygon();
          surface.zMin =  zGrade;
          surface.zMax =  zMax;
          surface.boundaryConditionType = Surface::EXTERIOR_TEMPERATURE;
          surface.orientation =  Surface::Y_POS;
          surfaces.push_back(surface);
        }

      }

    }

    // Far Field Surfaces
    {
      // X Min
      Surface surface;
      surface.type = Surface::ST_FAR_FIELD;
      surface.xMin =  xMin;
      surface.xMax =  xMin;
      surface.yMin =  yMin;
      surface.yMax =  yMax;
      surface.setSquarePolygon();
      surface.zMin =  zMin;
      surface.zMax =  zGrade;
      surface.boundaryConditionType = Surface::ZERO_FLUX;
      surface.orientation =  Surface::X_NEG;
      surfaces.push_back(surface);
    }

    {
      // X Max
      Surface surface;
      surface.type = Surface::ST_FAR_FIELD;
      surface.xMin =  xMax;
      surface.xMax =  xMax;
      surface.yMin =  yMin;
      surface.yMax =  yMax;
      surface.setSquarePolygon();
      surface.zMin =  zMin;
      surface.zMax =  zGrade;
      surface.boundaryConditionType = Surface::ZERO_FLUX;
      surface.orientation =  Surface::X_POS;
      surfaces.push_back(surface);
    }

    {
      // Y Min
      Surface surface;
      surface.type = Surface::ST_FAR_FIELD;
      surface.xMin =  xMin;
      surface.xMax =  xMax;
      surface.yMin =  yMin;
      surface.yMax =  yMin;
      surface.setSquarePolygon();
      surface.zMin =  zMin;
      surface.zMax =  zGrade;
      surface.boundaryConditionType = Surface::ZERO_FLUX;
      surface.orientation =  Surface::Y_NEG;
      surfaces.push_back(surface);
    }

    {
      // Y Max
      Surface surface;
      surface.type = Surface::ST_FAR_FIELD;
      surface.xMin =  xMin;
      surface.xMax =  xMax;
      surface.yMin =  yMax;
      surface.yMax =  yMax;
      surface.setSquarePolygon();
      surface.zMin =  zMin;
      surface.zMax =  zGrade;
      surface.boundaryConditionType = Surface::ZERO_FLUX;
      surface.orientation =  Surface::Y_POS;
      surfaces.push_back(surface);
    }

    // Deep ground surface
    if (deepGroundBoundary == DGB_CONSTANT_TEMPERATURE ||
      deepGroundBoundary == DGB_AUTO)
    {
      Surface surface;
      surface.type = Surface::ST_DEEP_GROUND;
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
      surface.type = Surface::ST_DEEP_GROUND;
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
      surface.type = Surface::ST_SLAB_CORE;
      surface.polygon = poly;
      surface.zMin = zSlab;
      surface.zMax = zSlab;
      surface.boundaryConditionType = Surface::INTERIOR_FLUX;
      surface.orientation = Surface::Z_POS;
      surface.emissivity = slab.emissivity;
      surfaces.push_back(surface);
    }
    if (hasPerimeterSurface)
    {
      Polygon poly;
      poly = offset(polygon, xySlabPerimeter);

      Polygon temp;
      temp = offset(polygon, xyPerimeterSurface);
      Ring ring;
      boost::geometry::convert(temp, ring);
      boost::geometry::reverse(ring);

      poly.inners().push_back(ring);

      Surface surface;
      surface.type = Surface::ST_SLAB_PERIM;
      surface.polygon = poly;
      surface.zMin = zSlab;
      surface.zMax = zSlab;
      surface.boundaryConditionType = Surface::INTERIOR_FLUX;
      surface.orientation = Surface::Z_POS;
      surface.emissivity = slab.emissivity;
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
      surface.type = Surface::ST_GRADE;
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

    if(foundationDepth > 0.0)
    {
      // Interior Air Top Surface
      Polygon poly;
      poly = offset(polygon, xyWallInterior);

      Surface surface;
      surface.type = Surface::ST_TOP_AIR_INT;
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
      {
        Polygon poly;
        poly = offset(polygon, xyWallExterior);

        Ring ring;
        boost::geometry::convert(poly, ring);
        boost::geometry::reverse(ring);

        Surface surface;
        surface.type = Surface::ST_TOP_AIR_EXT;
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

    }


    if (wallTopBoundary == WTB_LINEAR_DT)
    {
      // Wall Top
      if(hasWall)
      {
        double position = 0.0;
        double& Tin = wallTopInteriorTemperature;
        double& Tout = wallTopExteriorTemperature;
        double Tdiff = (Tin - Tout);
        std::size_t N = xyWallExterior/mesh.minCellDim;
        double temperature = Tin - (1.0/N)/2*Tdiff;

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
          surface.type = Surface::ST_WALL_TOP;
          surface.polygon = poly;
          surface.zMin = zMax;
          surface.zMax = zMax;
          surface.boundaryConditionType = Surface::CONSTANT_TEMPERATURE;
          surface.orientation = Surface::Z_POS;
          surface.temperature = temperature;
          surfaces.push_back(surface);

          position += xyWallExterior/N;
          temperature -= (1.0/N)*Tdiff;
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
        surface.type = Surface::ST_WALL_TOP;
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

      if (isGreaterThan(zIntVIns,zSlab))
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
        block.material = air;
        block.blockType = Block::INTERIOR_AIR;
        block.polygon = poly;
        block.zMin = zSlab;
        block.zMax = zIntVIns;
        blocks.push_back(block);
      }


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
  else if(numberOfDimensions == 3 && useSymmetry)
  {
    // TODO: 3D Symmetric
    isXSymm = isXSymmetric(polygon);
    isYSymm = isYSymmetric(polygon);

    // Change polygon to minimum symmetric unit
    Polygon symmetricPoly;
    symmetricPoly = symmetricUnit(polygon);

    nV = symmetricPoly.outer().size();

    Box boundingBox;
    boost::geometry::envelope(symmetricPoly, boundingBox);

    double xMinBB = boundingBox.min_corner().get<0>();
    double yMinBB = boundingBox.min_corner().get<1>();

    double xMaxBB = boundingBox.max_corner().get<0>();
    double yMaxBB = boundingBox.max_corner().get<1>();

    if (isXSymm)
      xMin = xMinBB;
    else
      xMin = xMinBB - farFieldWidth;

    if (isYSymm)
      yMin = yMinBB;
    else
      yMin = yMinBB - farFieldWidth;

    xMax = xMaxBB + farFieldWidth;
    yMax = yMaxBB + farFieldWidth;

    Box domainBox(Point(xMin,yMin),Point(xMax,yMax));

    if(isGreaterThan(foundationDepth,0.0))
    {
      // Interior Wall Surface
      if (hasInteriorVerticalInsulation && isGreaterThan(zIntVIns,zSlab))
      {
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
              surface.type = Surface::ST_WALL_INT;
              surface.xMin = getXmin(poly[0],v);
              surface.xMax = getXmax(poly[0],v);
              surface.yMin = getYmin(poly[0],v);
              surface.yMax = getYmax(poly[0],v);
              surface.setSquarePolygon();
              surface.zMin = zIntVIns;
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

          Surface surface;
          surface.type = Surface::ST_WALL_INT;
          surface.polygon = poly[0];
          surface.zMin = zIntVIns;
          surface.zMax = zIntVIns;
          surface.boundaryConditionType = Surface::INTERIOR_FLUX;
          surface.orientation = Surface::Z_NEG;
          surface.emissivity = wall.interiorEmissivity;
          surfaces.push_back(surface);
        }
        {
          Polygon tempPoly;
          tempPoly = polygon;
          MultiPolygon poly;
          boost::geometry::intersection(domainBox,tempPoly,poly);

          for (std::size_t v = 0; v < nV; v++)
          {

            if (!((isEqual(getXmin(poly[0],v),xMin) && isEqual(getXmax(poly[0],v),xMin)) ||
                (isEqual(getYmin(poly[0],v),yMin) && isEqual(getYmax(poly[0],v),yMin))))
            {
              Surface surface;
              surface.type = Surface::ST_WALL_INT;
              surface.xMin = getXmin(poly[0],v);
              surface.xMax = getXmax(poly[0],v);
              surface.yMin = getYmin(poly[0],v);
              surface.yMax = getYmax(poly[0],v);
              surface.setSquarePolygon();
              surface.zMin = zSlab;
              surface.zMax = zIntVIns;
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
      else
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
            surface.type = Surface::ST_WALL_INT;
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
            surface.type = Surface::ST_WALL_EXT;
            surface.xMin = getXmin(poly[0],v);
            surface.xMax = getXmax(poly[0],v);
            surface.yMin = getYmin(poly[0],v);
            surface.yMax = getYmax(poly[0],v);
            surface.setSquarePolygon();
            surface.zMin = zGrade;
            surface.zMax = zMax;
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
      surface.type = Surface::ST_FAR_FIELD;
      surface.xMin =  xMin;
      surface.xMax =  xMin;
      surface.yMin =  yMin;
      surface.yMax =  yMax;
      surface.setSquarePolygon();
      surface.zMin =  zMin;
      surface.zMax =  zMax;
      surface.boundaryConditionType = Surface::ZERO_FLUX;
      surface.orientation = Surface::X_NEG;
      surfaces.push_back(surface);
    }

    {
      // X Max
      Surface surface;
      surface.type = Surface::ST_FAR_FIELD;
      surface.xMin =  xMax;
      surface.xMax =  xMax;
      surface.yMin =  yMin;
      surface.yMax =  yMax;
      surface.setSquarePolygon();
      surface.zMin =  zMin;
      surface.zMax =  zMax;
      surface.boundaryConditionType = Surface::ZERO_FLUX;
      surface.orientation = Surface::X_POS;
      surfaces.push_back(surface);
    }

    {
      // Y Min
      Surface surface;
      surface.type = Surface::ST_FAR_FIELD;
      surface.xMin =  xMin;
      surface.xMax =  xMax;
      surface.yMin =  yMin;
      surface.yMax =  yMin;
      surface.setSquarePolygon();
      surface.zMin =  zMin;
      surface.zMax =  zMax;
      surface.boundaryConditionType = Surface::ZERO_FLUX;
      surface.orientation =  Surface::Y_NEG;
      surfaces.push_back(surface);
    }

    {
      // Y Max
      Surface surface;
      surface.type = Surface::ST_FAR_FIELD;
      surface.xMin =  xMin;
      surface.xMax =  xMax;
      surface.yMin =  yMax;
      surface.yMax =  yMax;
      surface.setSquarePolygon();
      surface.zMin =  zMin;
      surface.zMax =  zMax;
      surface.boundaryConditionType = Surface::ZERO_FLUX;
      surface.orientation =  Surface::Y_POS;
      surfaces.push_back(surface);
    }

    // Deep ground surface
    if (deepGroundBoundary == DGB_CONSTANT_TEMPERATURE ||
      deepGroundBoundary == DGB_AUTO)
    {
      Surface surface;
      surface.type = Surface::ST_DEEP_GROUND;
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
      surface.type = Surface::ST_DEEP_GROUND;
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
      surface.type = Surface::ST_SLAB_CORE;
      surface.polygon = poly[0];
      surface.zMin = zSlab;
      surface.zMax = zSlab;
      surface.boundaryConditionType = Surface::INTERIOR_FLUX;
      surface.orientation = Surface::Z_POS;
      surface.emissivity = slab.emissivity;
      surfaces.push_back(surface);
    }
    if (hasPerimeterSurface)
    {
      Polygon tempPoly;
      tempPoly = offset(polygon, xySlabPerimeter);

      Polygon temp;
      temp = offset(polygon, xyPerimeterSurface);
      Ring ring;
      boost::geometry::convert(temp, ring);
      boost::geometry::reverse(ring);

      tempPoly.inners().push_back(ring);

      MultiPolygon poly;
      boost::geometry::intersection(domainBox,tempPoly,poly);

      Surface surface;
      surface.type = Surface::ST_SLAB_PERIM;
      surface.polygon = poly[0];
      surface.zMin = zSlab;
      surface.zMax = zSlab;
      surface.boundaryConditionType = Surface::INTERIOR_FLUX;
      surface.orientation = Surface::Z_POS;
      surface.emissivity = slab.emissivity;
      surfaces.push_back(surface);
    }

    // Grade
    {
      Polygon tempPoly;
      tempPoly = offset(polygon, xyWallExterior);

      MultiPolygon poly;
      boost::geometry::difference(domainBox,tempPoly,poly);

      Surface surface;
      surface.type = Surface::ST_GRADE;
      surface.polygon = poly[0];
      surface.zMin = zGrade;
      surface.zMax = zGrade;
      surface.boundaryConditionType = Surface::EXTERIOR_FLUX;
      surface.orientation = Surface::Z_POS;
      surface.emissivity = soilEmissivity;
      surface.absorptivity = soilAbsorptivity;
      surfaces.push_back(surface);
    }

    if(foundationDepth > 0.0)
    {
      // Interior Air Top Surface
      Polygon tempPoly;
      tempPoly = offset(polygon, xyWallInterior);
      MultiPolygon poly;
      boost::geometry::intersection(domainBox,tempPoly,poly);

      Surface surface;
      surface.type = Surface::ST_TOP_AIR_INT;
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
      surface.type = Surface::ST_TOP_AIR_EXT;
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
        double& Tin = wallTopInteriorTemperature;
        double& Tout = wallTopExteriorTemperature;
        double Tdiff = (Tin - Tout);
        std::size_t N = xyWallExterior/mesh.minCellDim;
        double temperature = Tin - (1.0/N)/2*Tdiff;

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
          surface.type = Surface::ST_WALL_TOP;
          surface.polygon = poly[0];
          surface.zMin = zMax;
          surface.zMax = zMax;
          surface.boundaryConditionType = Surface::CONSTANT_TEMPERATURE;
          surface.orientation = Surface::Z_POS;
          surface.temperature = temperature;
          surfaces.push_back(surface);

          position += xyWallExterior/N;
          temperature -= (1.0/N)*Tdiff;
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
        surface.type = Surface::ST_WALL_TOP;
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

      if (isGreaterThan(zIntVIns,zSlab))
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
        block.material = air;
        block.blockType = Block::INTERIOR_AIR;
        block.polygon = poly[0];
        block.zMin = zSlab;
        block.zMax = zIntVIns;
        blocks.push_back(block);
      }

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
      // Make sure points are within the domain
      double pointX = std::max(std::min(surfaces[s].polygon.outer()[v].get<0>(),xMax),xMin);
      surfaces[s].polygon.outer()[v].set<0>(pointX);
      double pointY = std::max(std::min(surfaces[s].polygon.outer()[v].get<1>(),yMax),yMin);
      surfaces[s].polygon.outer()[v].set<1>(pointY);

      xPoints.push_back(surfaces[s].polygon.outer()[v].get<0>());
      yPoints.push_back(surfaces[s].polygon.outer()[v].get<1>());
    }
    surfaces[s].zMax = std::max(std::min(surfaces[s].zMax,zMax),zMin);
    surfaces[s].zMin = std::max(std::min(surfaces[s].zMin,zMax),zMin);

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
    for (std::size_t v = 0; v < blocks[b].polygon.outer().size(); v++)
    {
      // Make sure points are within the domain
      double pointX = std::max(std::min(blocks[b].polygon.outer()[v].get<0>(),xMax),xMin);
      blocks[b].polygon.outer()[v].set<0>(pointX);
      double pointY = std::max(std::min(blocks[b].polygon.outer()[v].get<1>(),yMax),yMin);
      blocks[b].polygon.outer()[v].set<1>(pointY);

      xPoints.push_back(blocks[b].polygon.outer()[v].get<0>());
      yPoints.push_back(blocks[b].polygon.outer()[v].get<1>());
    }
    blocks[b].zMax = std::max(std::min(blocks[b].zMax,zMax),zMin);
    blocks[b].zMin = std::max(std::min(blocks[b].zMin,zMax),zMin);

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

  // Make sure ranges are within the domain
  for (size_t r = 0; r < xRanges.ranges.size(); r++)
  {
    xRanges.ranges[r].range.first = std::max(std::min(xRanges.ranges[r].range.first,xMax),xMin);
    xRanges.ranges[r].range.second = std::max(std::min(xRanges.ranges[r].range.second,xMax),xMin);
  }
  for (size_t r = 0; r < yRanges.ranges.size(); r++)
  {
    yRanges.ranges[r].range.first = std::max(std::min(yRanges.ranges[r].range.first,yMax),yMin);
    yRanges.ranges[r].range.second = std::max(std::min(yRanges.ranges[r].range.second,yMax),yMin);
  }
  for (size_t r = 0; r < zRanges.ranges.size(); r++)
  {
    zRanges.ranges[r].range.first = std::max(std::min(zRanges.ranges[r].range.first,zMax),zMin);
    zRanges.ranges[r].range.second = std::max(std::min(zRanges.ranges[r].range.second,zMax),zMin);
  }

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

}
