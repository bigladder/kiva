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
#include "Geometry.h"

class SimulationControl
{
public:
  // Simulation Control
  boost::gregorian::date startDate;
  boost::gregorian::date endDate;
  boost::posix_time::time_duration timestep;
  std::string weatherFile;
  boost::posix_time::ptime startTime;

  void setStartTime();
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
  double heightAboveGrade;  // [m] below grade depth
  double height;  // [m] total height
  std::vector <Layer> layers;

  double totalWidth();
  double totalResistance();
};

class Slab
{
public:

  double emissivity;
  std::vector <Layer> layers;

  double totalWidth();
  double totalResistance();
};

class Mesh
{
public:

  double maxNearGrowthCoeff;
  double maxExteriorGrowthCoeff;
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
  bool contourLabels;
  std::string contourColor;
  bool gradients;
  bool axes;
  bool timestamp;
  int size;
  boost::gregorian::date startDate;
  boost::gregorian::date endDate;
  std::pair<double, double> xRange;
  std::pair<double, double> yRange;
  std::pair<double, double> zRange;

  enum PlotType
  {
    P_TEMP,
    P_FLUX
  };

  PlotType plotType;

  enum FluxDir
  {
    D_M,
    D_X,
    D_Y,
    D_Z
  };
  FluxDir fluxDir;

  enum ColorScheme
  {
    C_CMR,
    C_JET,
    C_NONE
  };
  ColorScheme colorScheme;

  enum Format
  {
    F_PNG,
    F_TEX
  };
  Format format;

  enum OutputUnits
  {
    IP,
    SI
  };
  OutputUnits outputUnits;

  double minimumTemperature;
  double maximumTemperature;

  int numberOfContours;

  bool startDateSet;
  bool endDateSet;
  bool xRangeSet;
  bool yRangeSet;
  bool zRangeSet;

};

class OutputVariable
{
private:
  std::vector<std::string> headers;

public:
  int variableID;
  std::string headerText;

  OutputVariable(int varID);

};

class OutputReport: public std::vector<OutputVariable>
{
public:
  boost::posix_time::time_duration minFrequency;
};

class Block
{
public:

  Polygon polygon;
  double xMin, xMax, yMin, yMax, zMin, zMax;
  Material material;

  enum BlockType
  {
    SOLID,
    INTERIOR_AIR,
    EXTERIOR_AIR
  };

  BlockType blockType;

  void setSquarePolygon();

};

class Surface
{
public:

  Polygon polygon;
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

  std::vector<boost::tuple<std::size_t,std::size_t,std::size_t> > indices;

  double area;

  void setSquarePolygon();
};



class RangeType
{
public:

  typedef std::pair<double,double> Range;

  enum Type
  {
    MIN_INTERIOR,
    MID_INTERIOR,
    MIN_EXTERIOR,
    MAX_EXTERIOR,
    DEEP,
    NEAR
  };
  Type type;
  Range range;
};

class Ranges
{
public:
  std::vector<RangeType> ranges;

  bool isType(double position,RangeType::Type type);

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
  double foundationDepth; // [m] below top of wall
  double orientation;  // [radians] from north

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
  double surfaceRoughness;  // dimensionless (roughly corresponsds to millimeters of relief)

  // Local wind speed characteristics
  double vegetationHeight;  // [m]
  double deltaLocal;  // [m]
  double alphaLocal;  // [-]


  // Geometry
  enum CoordinateSystem
  {
    CS_CARTESIAN,
    CS_CYLINDRICAL
  };
  CoordinateSystem coordinateSystem;

  int numberOfDimensions;  // 2 or 3

  bool useSymmetry;

  enum ReductionStrategy
  {
    RS_AP,
    RS_AP_APNEG,
    RS_RR,
    RS_A_P,
    RS_AP_PNEG,
    RS_CUSTOM,
    RS_BOUNDARY
  };
  ReductionStrategy reductionStrategy;

  bool twoParameters;
  double reductionLength1;
  double reductionLength2;

  double linearAreaMultiplier;

  Polygon polygon;
  bool isXSymm, isYSymm;

  double buildingHeight;
  std::vector<Polygon3> buildingSurfaces;

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

  double perimeterSurfaceWidth;
  bool hasPerimeterSurface;

  // Meshing
  Mesh mesh;


  // Simulation Control
  enum NumericalScheme
  {
    NS_ADE,
    NS_EXPLICIT,
    NS_ADI,
    NS_IMPLICIT,
    NS_CRANK_NICOLSON,
    NS_STEADY_STATE
  };

  NumericalScheme numericalScheme;

  double fADI;  // ADI modified f-factor

  std::string solver;
  std::string preconditioner;
  double tolerance;
  int maxIterations;

  double initialTemperature;
  enum InitializationMethod
  {
    IM_KUSUDA,
    IM_CONSTANT_TEMPERATURE,
    IM_STEADY_STATE
  };

  long warmupDays;
  long implicitAccelTimestep;
  long implicitAccelPeriods;

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

  double wallTopTemperatureDifference;
  enum WallTopBoundary
  {
    WTB_ZERO_FLUX,
    WTB_LINEAR_DT
  };

  WallTopBoundary wallTopBoundary;

  // Output Animations
  std::vector<OutputAnimation> outputAnimations;

  // Output Report
  OutputReport outputReport;

  // Derived variables
  MeshData xMeshData;
  MeshData yMeshData;
  MeshData zMeshData;
  std::vector<Block> blocks;
  std::vector<Surface> surfaces;

  void createMeshData();
};


class Input
{
public:
  SimulationControl simulationControl;
  std::vector <Foundation> foundations;
};

#endif /* INPUT_HPP_ */
