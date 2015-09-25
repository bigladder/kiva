/* InputParser.c++ is part of Kiva (Written by Neal Kruis)
 * Copyright (C) 2012-2015 Big Ladder Software <info@bigladdersoftware.com>
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

#include "Input.h"
#include "yaml-cpp/yaml.h"
#include "WeatherData.h"
#include <fstream>

Input inputParser(std::string inputFile)
{

  Input input;
  SimulationControl simulationControl;
  Foundation foundation1;

  YAML::Node yamlInput = YAML::LoadFile(inputFile);

  simulationControl.startDate =
      boost::gregorian::from_string(yamlInput["Simulation Control"]["Start Date"].as<std::string>());
  simulationControl.endDate =
      boost::gregorian::from_string(yamlInput["Simulation Control"]["End Date"].as<std::string>());
  simulationControl.timestep =
      boost::posix_time::minutes(yamlInput["Simulation Control"]["Timestep"].as<long>());

  // Materials
  std::map<std::string, Material> materials;

  for(YAML::const_iterator it=yamlInput["Materials"].begin();it!=yamlInput["Materials"].end();++it)
  {
    Material tempMaterial;
    tempMaterial.conductivity = it->second["Conductivity"].as<double>();
    tempMaterial.density = it->second["Density"].as<double>();
    tempMaterial.specificHeat = it->second["Specific Heat"].as<double>();

    materials.insert(std::pair<std::string,Material>(it->first.as<std::string>(),tempMaterial));
  }

  // Foundation

  // Soil
  foundation1.soil = materials[yamlInput["Foundation"]["Soil"].as<std::string>()];

  // Grade
  if (yamlInput["Foundation"]["Soil Absorptivity"].IsDefined()) {
    foundation1.soilAbsorptivity = yamlInput["Foundation"]["Soil Absorptivity"].as<double>();
  }
  else {
    foundation1.soilAbsorptivity = 0.8;
  }

  if (yamlInput["Foundation"]["Soil Emissivity"].IsDefined()) {
    foundation1.soilEmissivity = yamlInput["Foundation"]["Soil Emissivity"].as<double>();
  }
  else {
    foundation1.soilEmissivity = 0.8;
  }

  if (yamlInput["Foundation"]["Surface Roughness"].IsDefined()) {
    foundation1.surfaceRoughness = yamlInput["Foundation"]["Surface Roughness"].as<double>();
  }
  else {
    foundation1.surfaceRoughness = 30;
  }

  // Local wind speed characteristics
  if (yamlInput["Foundation"]["Vegetation Height"].IsDefined()) {
    foundation1.vegetationHeight = yamlInput["Foundation"]["Vegetation Height"].as<double>();
  }
  else {
    foundation1.vegetationHeight = 0.3;
  }

  if (yamlInput["Foundation"]["Delta Local"].IsDefined()) {
    foundation1.deltaLocal = yamlInput["Foundation"]["Delta Local"].as<double>();
  }
  else {
    foundation1.deltaLocal = 370;
  }

  if (yamlInput["Foundation"]["Alpha Local"].IsDefined()) {
    foundation1.alphaLocal = yamlInput["Foundation"]["Alpha Local"].as<double>();
  }
  else {
    foundation1.alphaLocal = 0.22;
  }

  // Slab
  if  (yamlInput["Foundation"]["Slab"].IsDefined())
  {
    foundation1.hasSlab = true;

    for (size_t i=0;i<yamlInput["Foundation"]["Slab"]["Layers"].size();i++)
    {

      Layer tempLayer;
      tempLayer.thickness = yamlInput["Foundation"]["Slab"]["Layers"][i]["Thickness"].as<double>();
      tempLayer.material = materials[yamlInput["Foundation"]["Slab"]["Layers"][i]["Material"].as<std::string>()];

      foundation1.slab.layers.push_back(tempLayer);

    }

    if (yamlInput["Foundation"]["Slab"]["Emissivity"].IsDefined()) {
      foundation1.slab.emissivity = yamlInput["Foundation"]["Slab"]["Emissivity"].as<double>();
    }
    else {
      foundation1.slab.emissivity = 0.8;
    }
  }
  else
  {
    foundation1.hasSlab = false;
  }

  // Wall
  if  (yamlInput["Foundation"]["Wall"].IsDefined())
  {
    foundation1.hasWall = true;

    for (size_t i=0;i<yamlInput["Foundation"]["Wall"]["Layers"].size();i++)
    {

      Layer tempLayer;
      tempLayer.thickness = yamlInput["Foundation"]["Wall"]["Layers"][i]["Thickness"].as<double>();
      tempLayer.material = materials[yamlInput["Foundation"]["Wall"]["Layers"][i]["Material"].as<std::string>()];

      foundation1.wall.layers.push_back(tempLayer);

    }

    foundation1.wall.heightAboveGrade = yamlInput["Foundation"]["Wall"]["Height Above Grade"].as<double>();
    foundation1.wall.height = yamlInput["Foundation"]["Wall"]["Height"].as<double>();

    if (yamlInput["Foundation"]["Wall"]["Interior Emissivity"].IsDefined()) {
      foundation1.wall.interiorEmissivity = yamlInput["Foundation"]["Wall"]["Interior Emissivity"].as<double>();
    }
    else {
      foundation1.wall.interiorEmissivity = 0.8;
    }

    if (yamlInput["Foundation"]["Wall"]["Exterior Emissivity"].IsDefined()) {
      foundation1.wall.exteriorEmissivity = yamlInput["Foundation"]["Wall"]["Exterior Emissivity"].as<double>();
    }
    else {
      foundation1.wall.exteriorEmissivity = 0.8;
    }

    if (yamlInput["Foundation"]["Wall"]["Exterior Absorptivity"].IsDefined()) {
      foundation1.wall.exteriorAbsorptivity = yamlInput["Foundation"]["Wall"]["Exterior Absorptivity"].as<double>();
    }
    else {
      foundation1.wall.exteriorAbsorptivity = 0.8;
    }

  }
  else
  {
    foundation1.hasWall = false;
  }

  // Interior Horizontal Insulation
  if  (yamlInput["Foundation"]["Interior Horizontal Insulation"].IsDefined())
  {
    foundation1.hasInteriorHorizontalInsulation = true;

    foundation1.interiorHorizontalInsulation.layer.thickness = yamlInput["Foundation"]["Interior Horizontal Insulation"]["Thickness"].as<double>();
    foundation1.interiorHorizontalInsulation.layer.material = materials[yamlInput["Foundation"]["Interior Horizontal Insulation"]["Material"].as<std::string>()];
    foundation1.interiorHorizontalInsulation.depth = yamlInput["Foundation"]["Interior Horizontal Insulation"]["Depth"].as<double>();
    foundation1.interiorHorizontalInsulation.width = yamlInput["Foundation"]["Interior Horizontal Insulation"]["Width"].as<double>();

  }
  else
  {
    foundation1.hasInteriorHorizontalInsulation = false;
  }

  // Interior Vertical Insulation
  if  (yamlInput["Foundation"]["Interior Vertical Insulation"].IsDefined())
  {
    foundation1.hasInteriorVerticalInsulation = true;

    foundation1.interiorVerticalInsulation.layer.thickness = yamlInput["Foundation"]["Interior Vertical Insulation"]["Thickness"].as<double>();
    foundation1.interiorVerticalInsulation.layer.material = materials[yamlInput["Foundation"]["Interior Vertical Insulation"]["Material"].as<std::string>()];
    foundation1.interiorVerticalInsulation.depth = yamlInput["Foundation"]["Interior Vertical Insulation"]["Depth"].as<double>();

  }
  else
  {
    foundation1.hasInteriorVerticalInsulation = false;
  }

  // Exterior Horizontal Insulation
  if  (yamlInput["Foundation"]["Exterior Horizontal Insulation"].IsDefined())
  {
    foundation1.hasExteriorHorizontalInsulation = true;

    foundation1.exteriorHorizontalInsulation.layer.thickness = yamlInput["Foundation"]["Exterior Horizontal Insulation"]["Thickness"].as<double>();
    foundation1.exteriorHorizontalInsulation.layer.material = materials[yamlInput["Foundation"]["Exterior Horizontal Insulation"]["Material"].as<std::string>()];
    foundation1.exteriorHorizontalInsulation.depth = yamlInput["Foundation"]["Exterior Horizontal Insulation"]["Depth"].as<double>();
    foundation1.exteriorHorizontalInsulation.width = yamlInput["Foundation"]["Exterior Horizontal Insulation"]["Width"].as<double>();

  }
  else
  {
    foundation1.hasExteriorHorizontalInsulation = false;
  }

  // Exterior Vertical Insulation
  if  (yamlInput["Foundation"]["Exterior Vertical Insulation"].IsDefined())
  {
    foundation1.hasExteriorVerticalInsulation = true;

    foundation1.exteriorVerticalInsulation.layer.thickness = yamlInput["Foundation"]["Exterior Vertical Insulation"]["Thickness"].as<double>();
    foundation1.exteriorVerticalInsulation.layer.material = materials[yamlInput["Foundation"]["Exterior Vertical Insulation"]["Material"].as<std::string>()];
    foundation1.exteriorVerticalInsulation.depth = yamlInput["Foundation"]["Exterior Vertical Insulation"]["Depth"].as<double>();

  }
  else
  {
    foundation1.hasExteriorVerticalInsulation = false;
  }


  // Site

  foundation1.foundationDepth = yamlInput["Foundation"]["Foundation Depth"].as<double>();

  if  (yamlInput["Foundation"]["Far-Field Width"].IsDefined()) {
    foundation1.farFieldWidth = yamlInput["Foundation"]["Far-Field Width"].as<double>();
  }
  else {
    foundation1.farFieldWidth = 40;
  }

  if  (yamlInput["Foundation"]["Deep-Ground Depth"].IsDefined()) {
    foundation1.deepGroundDepth = yamlInput["Foundation"]["Deep-Ground Depth"].as<double>();
  }
  else {
    foundation1.deepGroundDepth = 40;
  }

  if  (yamlInput["Foundation"]["Orientation"].IsDefined()) {
    foundation1.orientation = yamlInput["Foundation"]["Orientation"].as<double>();
  }
  else {
    foundation1.orientation = 0.0;
  }

  if (yamlInput["Foundation"]["Deep-Ground Boundary Condition"].IsDefined()) {
    if (yamlInput["Foundation"]["Deep-Ground Boundary Condition"].as<std::string>() == "AUTO")
    {
      foundation1.deepGroundBoundary = Foundation::DGB_AUTO;
    }
    else if (yamlInput["Foundation"]["Deep-Ground Boundary Condition"].as<std::string>() == "CONSTANT-TEMP")
    {
      foundation1.deepGroundBoundary = Foundation::DGB_CONSTANT_TEMPERATURE;
      foundation1.deepGroundTemperature = yamlInput["Foundation"]["Deep-Ground Temperature"].as<double>();
    }
    else if (yamlInput["Foundation"]["Deep-Ground Boundary Condition"].as<std::string>() == "ZERO-FLUX")
    {
      foundation1.deepGroundBoundary = Foundation::DGB_ZERO_FLUX;
    }
  }
  else {
    foundation1.deepGroundBoundary = Foundation::DGB_ZERO_FLUX;
  }

  if  (yamlInput["Foundation"]["Indoor Temperature Method"].IsDefined())
  {
    if (yamlInput["Foundation"]["Indoor Temperature Method"].as<std::string>() == "FILE")
    {
      foundation1.indoorTemperatureMethod = Foundation::ITM_FILE;
      foundation1.indoorAirTemperatureFile.fileName = yamlInput["Foundation"]["Indoor Air Temperature File"]["Name"].as<std::string>();
      foundation1.indoorAirTemperatureFile.firstIndex.first = yamlInput["Foundation"]["Indoor Air Temperature File"]["Index"][0].as<int>();
      foundation1.indoorAirTemperatureFile.firstIndex.second = yamlInput["Foundation"]["Indoor Air Temperature File"]["Index"][1].as<int>();
      foundation1.indoorAirTemperatureFile.readData();
    }
    else if (yamlInput["Foundation"]["Indoor Temperature Method"].as<std::string>() == "CONSTANT")
    {
      foundation1.indoorTemperatureMethod = Foundation::ITM_CONSTANT_TEMPERATURE;
      foundation1.indoorAirTemperature = yamlInput["Foundation"]["Indoor Air Temperature"].as<double>();
    }
  }
  else
  {
    foundation1.indoorTemperatureMethod = Foundation::ITM_CONSTANT_TEMPERATURE;
    foundation1.indoorAirTemperature = yamlInput["Foundation"]["Indoor Air Temperature"].as<double>();
  }

  // Geometry
  if (yamlInput["Foundation"]["Coordinate System"].IsDefined()) {
    if (yamlInput["Foundation"]["Coordinate System"].as<std::string>() == "CARTESIAN")
      foundation1.coordinateSystem = Foundation::CS_CARTESIAN;
    else if (yamlInput["Foundation"]["Coordinate System"].as<std::string>() == "CYLINDRICAL")
      foundation1.coordinateSystem = Foundation::CS_CYLINDRICAL;
  }
  else {
    foundation1.coordinateSystem = Foundation::CS_CARTESIAN;
  }

  if (yamlInput["Foundation"]["Two-Dimensional Approximation"].IsDefined()) {
    if (yamlInput["Foundation"]["Two-Dimensional Approximation"].as<std::string>() == "AP")
      foundation1.reductionStrategy = Foundation::RS_AP;
    else if (yamlInput["Foundation"]["Two-Dimensional Approximation"].as<std::string>() == "NEG")
      foundation1.reductionStrategy = Foundation::RS_AP_APNEG;
    else if (yamlInput["Foundation"]["Two-Dimensional Approximation"].as<std::string>() == "PNEG")
      foundation1.reductionStrategy = Foundation::RS_AP_PNEG;
    else if (yamlInput["Foundation"]["Two-Dimensional Approximation"].as<std::string>() == "RR")
      foundation1.reductionStrategy = Foundation::RS_RR;
    else if (yamlInput["Foundation"]["Two-Dimensional Approximation"].as<std::string>() == "A-P")
      foundation1.reductionStrategy = Foundation::RS_A_P;
    else if (yamlInput["Foundation"]["Two-Dimensional Approximation"].as<std::string>() == "BOUNDARY")
      foundation1.reductionStrategy = Foundation::RS_BOUNDARY;
    else if (yamlInput["Foundation"]["Two-Dimensional Approximation"].as<std::string>() == "CUSTOM")
    {
      foundation1.reductionStrategy = Foundation::RS_CUSTOM;
      if (yamlInput["Foundation"]["Length 1"].IsDefined())
      {
        foundation1.twoParameters = true;
        foundation1.reductionLength1 = yamlInput["Foundation"]["Length 1"].as<double>();
      }
      else
      {
        foundation1.twoParameters = false;
      }
      foundation1.reductionLength2 = yamlInput["Foundation"]["Length 2"].as<double>();
    }
  }
  else {
    foundation1.reductionStrategy = Foundation::RS_BOUNDARY;
  }

  if (yamlInput["Foundation"]["Number of Dimensions"].IsDefined())
    foundation1.numberOfDimensions = yamlInput["Foundation"]["Number of Dimensions"].as<int>();
  else
    foundation1.numberOfDimensions = 2;

  if  (yamlInput["Foundation"]["Use Symmetry"].IsDefined())
    foundation1.useSymmetry = yamlInput["Foundation"]["Use Symmetry"].as<bool>();
  else
    foundation1.useSymmetry = true;

  for (size_t i=0;i<yamlInput["Foundation"]["Polygon"].size();i++)
  {
    foundation1.polygon.outer().push_back(Point(
        yamlInput["Foundation"]["Polygon"][i][0].as<double>(),
        yamlInput["Foundation"]["Polygon"][i][1].as<double>()));
  }

  if  (yamlInput["Foundation"]["Building Height"].IsDefined())
  {
    foundation1.buildingHeight = yamlInput["Foundation"]["Building Height"].as<double>();
  }
  else
  {
    foundation1.buildingHeight = 0.0;
  }

  if  (yamlInput["Foundation"]["Perimeter Surface Width"].IsDefined())
  {
    foundation1.hasPerimeterSurface = true;
    foundation1.perimeterSurfaceWidth = yamlInput["Foundation"]["Perimeter Surface Width"].as<double>();
  }
  else
  {
    foundation1.hasPerimeterSurface = false;
  }

  // Meshing
  if  (yamlInput["Foundation"]["Mesh"].IsDefined())
  {
    foundation1.mesh.minCellDim = yamlInput["Foundation"]["Mesh"]["Minimum Cell Dimension"].as<double>();
    foundation1.mesh.maxNearGrowthCoeff = yamlInput["Foundation"]["Mesh"]["Maximum Near-Field Growth Coefficient"].as<double>();
    foundation1.mesh.maxDepthGrowthCoeff = yamlInput["Foundation"]["Mesh"]["Maximum Deep-Field Growth Coefficient"].as<double>();
    foundation1.mesh.maxInteriorGrowthCoeff = yamlInput["Foundation"]["Mesh"]["Maximum Interior-Field Growth Coefficient"].as<double>();
    foundation1.mesh.maxExteriorGrowthCoeff = yamlInput["Foundation"]["Mesh"]["Maximum Far-Field Growth Coefficient"].as<double>();
  }
  else
  {
    foundation1.mesh.minCellDim = 0.02;
    foundation1.mesh.maxNearGrowthCoeff = 1.5;
    foundation1.mesh.maxDepthGrowthCoeff = 1.5;
    foundation1.mesh.maxInteriorGrowthCoeff = 1.5;
    foundation1.mesh.maxExteriorGrowthCoeff = 1.5;
  }

  // Simulation Control
  if  (yamlInput["Foundation"]["Numerical Scheme"].IsDefined())
  {
    if (yamlInput["Foundation"]["Numerical Scheme"].as<std::string>() == "ADE")
      foundation1.numericalScheme = Foundation::NS_ADE;
    else if (yamlInput["Foundation"]["Numerical Scheme"].as<std::string>() == "EXPLICIT")
      foundation1.numericalScheme = Foundation::NS_EXPLICIT;
    else if (yamlInput["Foundation"]["Numerical Scheme"].as<std::string>() == "ADI")
    {
      foundation1.numericalScheme = Foundation::NS_ADI;
      if  (yamlInput["Foundation"]["f-ADI"].IsDefined())
        foundation1.fADI = yamlInput["Foundation"]["f-ADI"].as<double>();
      else
        foundation1.fADI = 0.00001;
    }
    else if (yamlInput["Foundation"]["Numerical Scheme"].as<std::string>() == "IMPLICIT")
      foundation1.numericalScheme = Foundation::NS_IMPLICIT;
    else if (yamlInput["Foundation"]["Numerical Scheme"].as<std::string>() == "CRANK-NICOLSON")
      foundation1.numericalScheme = Foundation::NS_CRANK_NICOLSON;
    else if (yamlInput["Foundation"]["Numerical Scheme"].as<std::string>() == "STEADY-STATE")
      foundation1.numericalScheme = Foundation::NS_STEADY_STATE;
  }
  else
  {
    foundation1.numericalScheme = Foundation::NS_ADI;
    foundation1.fADI = 0.00001;
}

  if  (yamlInput["Foundation"]["Solver"].IsDefined())
  {
    foundation1.solver = yamlInput["Foundation"]["Solver"].as<std::string>();
  }
  else
  {
    foundation1.solver = "bicgstab";
  }

  if  (yamlInput["Foundation"]["Preconditioner"].IsDefined())
  {
    foundation1.preconditioner = yamlInput["Foundation"]["Preconditioner"].as<std::string>();
  }
  else
  {
    foundation1.preconditioner = "ilu";
  }

  if  (yamlInput["Foundation"]["Maximum Iterations"].IsDefined())
  {
    foundation1.maxIterations = yamlInput["Foundation"]["Maximum Iterations"].as<int>();
  }
  else
  {
    foundation1.maxIterations = 100000;
  }

  if  (yamlInput["Foundation"]["Tolerance"].IsDefined())
  {
    foundation1.tolerance = yamlInput["Foundation"]["Tolerance"].as<double>();
  }
  else
  {
    foundation1.tolerance = 1.0e-6;
  }

  if  (yamlInput["Foundation"]["Initialization Method"].IsDefined())
  {
    if (yamlInput["Foundation"]["Initialization Method"].as<std::string>() == "KUSUDA")
      foundation1.initializationMethod = Foundation::IM_KUSUDA;
    else if (yamlInput["Foundation"]["Initialization Method"].as<std::string>() == "STEADY-STATE")
      foundation1.initializationMethod = Foundation::IM_STEADY_STATE;
    else if (yamlInput["Foundation"]["Initialization Method"].as<std::string>() == "CONSTANT")
    {
      foundation1.initializationMethod = Foundation::IM_CONSTANT_TEMPERATURE;
      foundation1.initialTemperature = yamlInput["Foundation"]["Initial Temperature"].as<double>();
    }
  }
  else
  {
    foundation1.initializationMethod = Foundation::IM_STEADY_STATE;
  }

  if  (yamlInput["Foundation"]["Accelerated Initialization Timestep"].IsDefined())
  {
    foundation1.implicitAccelTimestep = yamlInput["Foundation"]["Accelerated Initialization Timestep"].as<long>();
  }
  else
  {
    foundation1.implicitAccelTimestep = 168;
  }

  if  (yamlInput["Foundation"]["Number of Accelearted Initialization Timesteps"].IsDefined())
  {
    foundation1.implicitAccelPeriods = yamlInput["Foundation"]["Number of Accelearted Initialization Timesteps"].as<long>();
  }
  else
  {
    foundation1.implicitAccelPeriods = 12;
  }

  if  (yamlInput["Foundation"]["Number of Warmup Days in Initialization"].IsDefined())
  {
    foundation1.warmupDays = yamlInput["Foundation"]["Number of Warmup Days in Initialization"].as<long>();
  }
  else
  {
    foundation1.warmupDays = 365;
  }


  if  (yamlInput["Foundation"]["Convection Calculation Method"].IsDefined())
  {
    if (yamlInput["Foundation"]["Convection Calculation Method"].as<std::string>() == "AUTO")
      foundation1.convectionCalculationMethod = Foundation::CCM_AUTO;
    else if (yamlInput["Foundation"]["Convection Calculation Method"].as<std::string>() == "CONSTANT")
    {
      foundation1.convectionCalculationMethod = Foundation::CCM_CONSTANT_COEFFICIENT;
      foundation1.interiorConvectiveCoefficient = yamlInput["Foundation"]["Interior Convective Coefficient"].as<double>();
      foundation1.exteriorConvectiveCoefficient = yamlInput["Foundation"]["Exterior Convective Coefficient"].as<double>();
    }
  }
  else
  {
    foundation1.convectionCalculationMethod = Foundation::CCM_AUTO;
  }

  if  (yamlInput["Foundation"]["Outdoor Temperature Method"].IsDefined())
  {
    if (yamlInput["Foundation"]["Outdoor Temperature Method"].as<std::string>() == "WEATHER-FILE")
      foundation1.outdoorTemperatureMethod = Foundation::OTM_WEATHER_FILE;
    else if (yamlInput["Foundation"]["Outdoor Temperature Method"].as<std::string>() == "CONSTANT")
    {
      foundation1.outdoorTemperatureMethod = Foundation::OTM_CONSTANT_TEMPERATURE;
      foundation1.outdoorDryBulbTemperature = yamlInput["Foundation"]["Outdoor Dry-Bulb Temperature"].as<double>();
    }
  }
  else
  {
    foundation1.outdoorTemperatureMethod = Foundation::OTM_WEATHER_FILE;
  }

  if  (yamlInput["Foundation"]["Wall Top Boundary Condition"].IsDefined())
  {
    if (yamlInput["Foundation"]["Wall Top Boundary Condition"].as<std::string>() == "ZERO-FLUX")
      foundation1.wallTopBoundary = Foundation::WTB_ZERO_FLUX;
    else if (yamlInput["Foundation"]["Wall Top Boundary Condition"].as<std::string>() == "LINEAR-DT")
    {
      foundation1.wallTopBoundary = Foundation::WTB_LINEAR_DT;
      foundation1.wallTopTemperatureDifference = yamlInput["Foundation"]["Wall Top Temperature Difference"].as<double>();
    }
  }
  else
  {
    foundation1.wallTopBoundary = Foundation::WTB_ZERO_FLUX;
  }

  // Output

  // CSV Reports
  for(size_t i=0;i<yamlInput["Foundation"]["Output Report"]["Reports"].size();i++)
  {
    OutputVariable temp(yamlInput["Foundation"]["Output Report"]["Reports"][i].as<int>());
    foundation1.outputReport.push_back(temp);
  }

  if  (yamlInput["Foundation"]["Output Report"]["Minimum Reporting Frequency"].IsDefined())
  {
    foundation1.outputReport.minFrequency =
        boost::posix_time::minutes(
            yamlInput["Foundation"]["Output Report"]["Minimum Reporting Frequency"].as<long>());
  }
  else
  {
    foundation1.outputReport.minFrequency = boost::posix_time::minutes(60);
  }

  // Animations/Plots
  for(size_t i=0;i<yamlInput["Foundation"]["Output Snapshots"].size();i++)
  {
    OutputAnimation temp;
    temp.name = yamlInput["Foundation"]["Output Snapshots"][i]["Name"].as<std::string>();
    temp.frequency = boost::posix_time::hours(yamlInput["Foundation"]["Output Snapshots"][i]["Frequency"].as<long>());
    temp.grid = yamlInput["Foundation"]["Output Snapshots"][i]["Grid"].as<bool>();
    temp.gradients = yamlInput["Foundation"]["Output Snapshots"][i]["Gradients"].as<bool>();
    temp.contours = yamlInput["Foundation"]["Output Snapshots"][i]["Contours"].as<bool>();
    temp.contourLabels = yamlInput["Foundation"]["Output Snapshots"][i]["Contour Labels"].as<bool>();
    temp.axes = yamlInput["Foundation"]["Output Snapshots"][i]["Axes"].as<bool>();
    temp.timestamp = yamlInput["Foundation"]["Output Snapshots"][i]["Timestamp"].as<bool>();

    if (yamlInput["Foundation"]["Output Snapshots"][i]["Plot Type"].IsDefined())
    {
      if (yamlInput["Foundation"]["Output Snapshots"][i]["Plot Type"].as<std::string>() == "TEMPERATURE")
        temp.plotType = OutputAnimation::P_TEMP;
      else if (yamlInput["Foundation"]["Output Snapshots"][i]["Plot Type"].as<std::string>() == "HEAT-FLUX")
        temp.plotType = OutputAnimation::P_FLUX;
    }
    else
      temp.plotType = OutputAnimation::P_TEMP;

    if (yamlInput["Foundation"]["Output Snapshots"][i]["Flux Direction"].IsDefined())
    {
      if (yamlInput["Foundation"]["Output Snapshots"][i]["Flux Direction"].as<std::string>() == "MAG")
        temp.fluxDir = OutputAnimation::D_M;
      else if (yamlInput["Foundation"]["Output Snapshots"][i]["Flux Direction"].as<std::string>() == "X")
        temp.fluxDir = OutputAnimation::D_X;
      else if (yamlInput["Foundation"]["Output Snapshots"][i]["Flux Direction"].as<std::string>() == "Y")
        temp.fluxDir = OutputAnimation::D_Y;
      else if (yamlInput["Foundation"]["Output Snapshots"][i]["Flux Direction"].as<std::string>() == "Z")
        temp.fluxDir = OutputAnimation::D_Z;
    }
    else
      temp.fluxDir = OutputAnimation::D_M;

    if (yamlInput["Foundation"]["Output Snapshots"][i]["Color Scheme"].IsDefined())
    {
      if (yamlInput["Foundation"]["Output Snapshots"][i]["Color Scheme"].as<std::string>() == "CMR")
        temp.colorScheme = OutputAnimation::C_CMR;
      else if (yamlInput["Foundation"]["Output Snapshots"][i]["Color Scheme"].as<std::string>() == "JET")
        temp.colorScheme = OutputAnimation::C_JET;
      else if (yamlInput["Foundation"]["Output Snapshots"][i]["Color Scheme"].as<std::string>() == "NONE")
        temp.colorScheme = OutputAnimation::C_NONE;
    }
    else
      temp.colorScheme = OutputAnimation::C_CMR;

    if (yamlInput["Foundation"]["Output Snapshots"][i]["File Format"].IsDefined())
    {
      if (yamlInput["Foundation"]["Output Snapshots"][i]["File Format"].as<std::string>() == "PNG")
        temp.format = OutputAnimation::F_PNG;
      else if (yamlInput["Foundation"]["Output Snapshots"][i]["File Format"].as<std::string>() == "TEX")
        temp.format = OutputAnimation::F_TEX;
    }
    else
      temp.format = OutputAnimation::F_PNG;

    if (yamlInput["Foundation"]["Output Snapshots"][i]["Unit System"].IsDefined())
    {
      if (yamlInput["Foundation"]["Output Snapshots"][i]["Unit System"].as<std::string>() == "IP")
        temp.outputUnits = OutputAnimation::IP;
      else if (yamlInput["Foundation"]["Output Snapshots"][i]["Unit System"].as<std::string>() == "SI")
        temp.outputUnits = OutputAnimation::SI;
    }
    else
      temp.outputUnits = OutputAnimation::SI;

    if (yamlInput["Foundation"]["Output Snapshots"][i]["Output Range"].IsDefined())
    {
      temp.minimumTemperature = yamlInput["Foundation"]["Output Snapshots"][i]["Output Range"][0].as<double>();
      temp.maximumTemperature = yamlInput["Foundation"]["Output Snapshots"][i]["Output Range"][1].as<double>();
    }
    else
    {
      temp.minimumTemperature = -20;
      temp.maximumTemperature = 40;
    }

    if (yamlInput["Foundation"]["Output Snapshots"][i]["Number of Contours"].IsDefined())
    {
      temp.numberOfContours = yamlInput["Foundation"]["Output Snapshots"][i]["Number of Contours"].as<int>();
    }
    else
    {
      temp.numberOfContours = 13;
    }

    if (yamlInput["Foundation"]["Output Snapshots"][i]["Contour Color"].IsDefined())
    {
      temp.contourColor = yamlInput["Foundation"]["Output Snapshots"][i]["Contour Color"].as<std::string>();
    }
    else
    {
      temp.contourColor = "H";
    }

    temp.size = yamlInput["Foundation"]["Output Snapshots"][i]["Size"].as<int>();

    if (yamlInput["Foundation"]["Output Snapshots"][i]["Start Date"].IsDefined())
    {
      temp.startDate = boost::gregorian::from_string(yamlInput["Foundation"]["Output Snapshots"][i]["Start Date"].as<std::string>());
      temp.startDateSet = true;
    }
    else
      temp.startDateSet = false;

    if (yamlInput["Foundation"]["Output Snapshots"][i]["End Date"].IsDefined())
    {
      temp.endDate = boost::gregorian::from_string(yamlInput["Foundation"]["Output Snapshots"][i]["End Date"].as<std::string>());
      temp.endDateSet = true;
    }
    else
      temp.endDateSet = false;

    if (yamlInput["Foundation"]["Output Snapshots"][i]["X Range"].IsDefined())
    {
      temp.xRange.first = yamlInput["Foundation"]["Output Snapshots"][i]["X Range"][0].as<double>();
      temp.xRange.second = yamlInput["Foundation"]["Output Snapshots"][i]["X Range"][1].as<double>();
      temp.xRangeSet = true;
    }
    else
      temp.xRangeSet = false;

    if (yamlInput["Foundation"]["Output Snapshots"][i]["Y Range"].IsDefined())
    {
      temp.yRange.first = yamlInput["Foundation"]["Output Snapshots"][i]["Y Range"][0].as<double>();
      temp.yRange.second = yamlInput["Foundation"]["Output Snapshots"][i]["Y Range"][1].as<double>();
      temp.yRangeSet = true;
    }
    else
      temp.yRangeSet = false;

    if (yamlInput["Foundation"]["Output Snapshots"][i]["Z Range"].IsDefined())
    {

      temp.zRange.first = yamlInput["Foundation"]["Output Snapshots"][i]["Z Range"][0].as<double>();
      temp.zRange.second = yamlInput["Foundation"]["Output Snapshots"][i]["Z Range"][1].as<double>();
      temp.zRangeSet = true;
    }
    else
      temp.zRangeSet = false;

    foundation1.outputAnimations.push_back(temp);
  }


  // Full Input
  input.simulationControl = simulationControl;
  input.foundations.push_back(foundation1);


  return input;

}
