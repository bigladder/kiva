/* Copyright (c) 2012-2022 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#include "InputParser.hpp"
#include "Errors.hpp"

Input inputParser(std::string inputFile) {

  Input input;
  SimulationControl simulationControl;
  Foundation foundation;
  Boundaries boundaries;
  Initialization initialization;
  Output output;

  std::filesystem::path inputPath(inputFile);

  YAML::Node yamlInput = YAML::LoadFile(inputPath.string());

  // SIMULATION CONTROL
  simulationControl.startDate = boost::gregorian::from_string(
      yamlInput["Simulation Control"]["Start Date"].as<std::string>());
  simulationControl.endDate =
      boost::gregorian::from_string(yamlInput["Simulation Control"]["End Date"].as<std::string>());
  simulationControl.timestep =
      boost::posix_time::minutes(yamlInput["Simulation Control"]["Timestep"].as<long>());

  // MATERIALS
  std::map<std::string, Material> materials;

  for (YAML::const_iterator it = yamlInput["Materials"].begin(); it != yamlInput["Materials"].end();
       ++it) {
    Material tempMaterial;
    tempMaterial.conductivity = it->second["Conductivity"].as<double>();
    tempMaterial.density = it->second["Density"].as<double>();
    tempMaterial.specificHeat = it->second["Specific Heat"].as<double>();

    materials.insert(std::pair<std::string, Material>(it->first.as<std::string>(), tempMaterial));
  }

  // FOUNDATION

  // Soil
  foundation.soil = materials[yamlInput["Foundation"]["Soil"].as<std::string>()];

  // Grade
  if (yamlInput["Foundation"]["Soil Absorptivity"].IsDefined()) {
    foundation.grade.absorptivity = yamlInput["Foundation"]["Soil Absorptivity"].as<double>();
  } else {
    foundation.grade.absorptivity = 0.8;
  }

  if (yamlInput["Foundation"]["Soil Emissivity"].IsDefined()) {
    foundation.grade.emissivity = yamlInput["Foundation"]["Soil Emissivity"].as<double>();
  } else {
    foundation.grade.emissivity = 0.8;
  }

  if (yamlInput["Foundation"]["Surface Roughness"].IsDefined()) {
    foundation.grade.roughness = yamlInput["Foundation"]["Surface Roughness"].as<double>();
  } else {
    foundation.grade.roughness = 0.03;
  }

  foundation.foundationDepth = yamlInput["Foundation"]["Foundation Depth"].as<double>();

  // Slab
  if (yamlInput["Foundation"]["Slab"].IsDefined()) {
    foundation.hasSlab = true;

    for (size_t i = 0; i < yamlInput["Foundation"]["Slab"]["Layers"].size(); i++) {

      Layer tempLayer;
      tempLayer.thickness = yamlInput["Foundation"]["Slab"]["Layers"][i]["Thickness"].as<double>();
      tempLayer.material =
          materials[yamlInput["Foundation"]["Slab"]["Layers"][i]["Material"].as<std::string>()];

      foundation.slab.layers.push_back(tempLayer);
    }

    if (yamlInput["Foundation"]["Slab"]["Emissivity"].IsDefined()) {
      foundation.slab.interior.emissivity =
          yamlInput["Foundation"]["Slab"]["Emissivity"].as<double>();
    } else {
      foundation.slab.interior.emissivity = 0.8;
    }
  } else {
    foundation.hasSlab = false;
  }

  // Wall
  if (yamlInput["Foundation"]["Wall"].IsDefined()) {
    foundation.hasWall = true;

    for (size_t i = 0; i < yamlInput["Foundation"]["Wall"]["Layers"].size(); i++) {

      Layer tempLayer;
      tempLayer.thickness = yamlInput["Foundation"]["Wall"]["Layers"][i]["Thickness"].as<double>();
      tempLayer.material =
          materials[yamlInput["Foundation"]["Wall"]["Layers"][i]["Material"].as<std::string>()];

      foundation.wall.layers.push_back(tempLayer);
    }

    if (yamlInput["Foundation"]["Wall"]["Height Above Grade"].IsDefined()) {
      foundation.wall.heightAboveGrade =
          yamlInput["Foundation"]["Wall"]["Height Above Grade"].as<double>();
    } else {
      foundation.wall.heightAboveGrade = 0.2;
    }

    if (yamlInput["Foundation"]["Wall"]["Depth Below Slab"].IsDefined()) {
      foundation.wall.depthBelowSlab =
          yamlInput["Foundation"]["Wall"]["Depth Below Slab"].as<double>();
    } else {
      foundation.wall.depthBelowSlab = 0.0;
    }

    if (yamlInput["Foundation"]["Wall"]["Interior Emissivity"].IsDefined()) {
      foundation.wall.interior.emissivity =
          yamlInput["Foundation"]["Wall"]["Interior Emissivity"].as<double>();
    } else {
      foundation.wall.interior.emissivity = 0.8;
    }

    if (yamlInput["Foundation"]["Wall"]["Exterior Emissivity"].IsDefined()) {
      foundation.wall.exterior.emissivity =
          yamlInput["Foundation"]["Wall"]["Exterior Emissivity"].as<double>();
    } else {
      foundation.wall.exterior.emissivity = 0.8;
    }

    if (yamlInput["Foundation"]["Wall"]["Exterior Absorptivity"].IsDefined()) {
      foundation.wall.exterior.absorptivity =
          yamlInput["Foundation"]["Wall"]["Exterior Absorptivity"].as<double>();
    } else {
      foundation.wall.exterior.absorptivity = 0.8;
    }

  } else {
    foundation.hasWall = false;
  }

  // Misc. Material blocks
  if (yamlInput["Foundation"]["Material Blocks"].IsDefined()) {
    for (auto block : yamlInput["Foundation"]["Material Blocks"]) {
      InputBlock tempBlock;

      if (block["X Position"].IsDefined()) {
        tempBlock.x = block["X Position"].as<double>();
      } else {
        showMessage(MSG_ERR, "'X Position' is required for Material Blocks.");
      }

      if (block["Z Position"].IsDefined()) {
        tempBlock.z = block["Z Position"].as<double>();
      } else {
        showMessage(MSG_ERR, "'Z Position' is required for Material Blocks.");
      }

      if (block["Width"].IsDefined()) {
        tempBlock.width = block["Width"].as<double>();
      } else {
        showMessage(MSG_ERR, "'Width' is required for Material Blocks.");
      }

      if (block["Depth"].IsDefined()) {
        tempBlock.depth = block["Depth"].as<double>();
      } else {
        showMessage(MSG_ERR, "'Depth' is required for Material Blocks.");
      }

      if (block["Material"].IsDefined()) {
        tempBlock.material = materials[block["Material"].as<std::string>()];
      } else {
        showMessage(MSG_ERR, "'Material' is required for Material Blocks.");
      }

      foundation.inputBlocks.push_back(tempBlock);
    }
  }

  // Interior Horizontal Insulation
  if (yamlInput["Foundation"]["Interior Horizontal Insulation"].IsDefined()) {

    InputBlock tempBlock;
    tempBlock.x = 0.0;

    if (yamlInput["Foundation"]["Interior Horizontal Insulation"]["Depth"].IsDefined()) {
      tempBlock.z = foundation.foundationDepth + foundation.slab.totalWidth() +
                    yamlInput["Foundation"]["Interior Horizontal Insulation"]["Depth"].as<double>();
    } else {
      tempBlock.z = foundation.foundationDepth + foundation.slab.totalWidth();
    }

    tempBlock.width =
        -yamlInput["Foundation"]["Interior Horizontal Insulation"]["Width"].as<double>();
    tempBlock.depth =
        yamlInput["Foundation"]["Interior Horizontal Insulation"]["Thickness"].as<double>();
    tempBlock.material =
        materials[yamlInput["Foundation"]["Interior Horizontal Insulation"]["Material"]
                      .as<std::string>()];

    foundation.inputBlocks.push_back(tempBlock);
  }

  // Interior Vertical Insulation
  if (yamlInput["Foundation"]["Interior Vertical Insulation"].IsDefined()) {
    InputBlock tempBlock;
    tempBlock.x = 0.0;
    tempBlock.z = 0.0;

    tempBlock.width =
        -yamlInput["Foundation"]["Interior Vertical Insulation"]["Thickness"].as<double>();
    tempBlock.depth = yamlInput["Foundation"]["Interior Vertical Insulation"]["Depth"].as<double>();
    tempBlock.material =
        materials[yamlInput["Foundation"]["Interior Vertical Insulation"]["Material"]
                      .as<std::string>()];

    foundation.inputBlocks.push_back(tempBlock);
  }

  // Exterior Horizontal Insulation
  if (yamlInput["Foundation"]["Exterior Horizontal Insulation"].IsDefined()) {
    InputBlock tempBlock;
    tempBlock.x = foundation.wall.totalWidth();

    if (yamlInput["Foundation"]["Exterior Horizontal Insulation"]["Depth"].IsDefined()) {
      tempBlock.z = foundation.wall.heightAboveGrade +
                    yamlInput["Foundation"]["Exterior Horizontal Insulation"]["Depth"].as<double>();
    } else {
      tempBlock.z = foundation.wall.heightAboveGrade;
    }

    tempBlock.width =
        yamlInput["Foundation"]["Exterior Horizontal Insulation"]["Width"].as<double>();
    tempBlock.depth =
        yamlInput["Foundation"]["Exterior Horizontal Insulation"]["Thickness"].as<double>();
    tempBlock.material =
        materials[yamlInput["Foundation"]["Exterior Horizontal Insulation"]["Material"]
                      .as<std::string>()];

    foundation.inputBlocks.push_back(tempBlock);
  }

  // Exterior Vertical Insulation
  if (yamlInput["Foundation"]["Exterior Vertical Insulation"].IsDefined()) {
    InputBlock tempBlock;
    tempBlock.x = foundation.wall.totalWidth();
    tempBlock.z = 0.0;

    tempBlock.width =
        yamlInput["Foundation"]["Exterior Vertical Insulation"]["Thickness"].as<double>();
    tempBlock.depth = yamlInput["Foundation"]["Exterior Vertical Insulation"]["Depth"].as<double>();
    tempBlock.material =
        materials[yamlInput["Foundation"]["Exterior Vertical Insulation"]["Material"]
                      .as<std::string>()];

    foundation.inputBlocks.push_back(tempBlock);
  }

  // Footing
  if (yamlInput["Foundation"]["Footing"].IsDefined()) {
    InputBlock tempBlock;

    if (yamlInput["Foundation"]["Footing"]["Width"].IsDefined()) {
      tempBlock.width = yamlInput["Foundation"]["Footing"]["Width"].as<double>();
    } else {
      showMessage(MSG_ERR, "'Width' is required for Footing.");
    }
    if (yamlInput["Foundation"]["Footing"]["Depth"].IsDefined()) {
      tempBlock.depth = yamlInput["Foundation"]["Footing"]["Depth"].as<double>();
    } else {
      showMessage(MSG_ERR, "'Depth' is required for Footing.");
    }

    tempBlock.z = foundation.foundationDepth +
                  (foundation.hasSlab ? foundation.slab.totalWidth() : 0.0) +
                  (foundation.hasWall ? foundation.wall.depthBelowSlab : 0.0);
    double center = (foundation.hasWall ? foundation.wall.totalWidth() / 2.0 : 0.0);
    tempBlock.x = center - tempBlock.width / 2.0; // Centered on wall

    if (yamlInput["Foundation"]["Footing"]["Material"].IsDefined()) {
      tempBlock.material =
          materials[yamlInput["Foundation"]["Footing"]["Material"].as<std::string>()];
    } else {
      showMessage(MSG_ERR, "'Material' is required for Footing.");
    }

    foundation.inputBlocks.push_back(tempBlock);
  }

  // Site
  if (yamlInput["Foundation"]["Orientation"].IsDefined()) {
    foundation.orientation = yamlInput["Foundation"]["Orientation"].as<double>();
  } else {
    foundation.orientation = 0.0;
  }

  // Geometry
  for (size_t i = 0; i < yamlInput["Foundation"]["Polygon"].size(); i++) {
    foundation.polygon.outer().push_back(
        Point(yamlInput["Foundation"]["Polygon"][i][0].as<double>(),
              yamlInput["Foundation"]["Polygon"][i][1].as<double>()));
  }

  if (yamlInput["Foundation"]["Exposed Perimeter"].IsDefined()) {
    foundation.useDetailedExposedPerimeter = true;
    if (yamlInput["Foundation"]["Exposed Perimeter"].size() !=
        yamlInput["Foundation"]["Polygon"].size()) {
      showMessage(MSG_ERR,
                  "'Exposed Perimeter' list must be the same size as 'Polygon' vetex list.");
    } else {
      for (size_t i = 0; i < yamlInput["Foundation"]["Exposed Perimeter"].size(); i++) {
        foundation.isExposedPerimeter.push_back(
            yamlInput["Foundation"]["Exposed Perimeter"][i].as<bool>());
      }
    }
  } else {
    foundation.useDetailedExposedPerimeter = true;
    if (!yamlInput["Foundation"]["Exposed Fraction"].IsDefined()) {
      for (size_t i = 0; i < foundation.polygon.outer().size(); i++) {
        foundation.isExposedPerimeter.push_back(true);
      }
    }
  }

  if (yamlInput["Foundation"]["Exposed Fraction"].IsDefined()) {
    if (yamlInput["Foundation"]["Exposed Perimeter"].IsDefined()) {
      showMessage(MSG_ERR, "Cannot specify both 'Exposed Fraction' and 'Exposed Perimeter'.");
    }
    foundation.useDetailedExposedPerimeter = false;
    foundation.exposedFraction = yamlInput["Foundation"]["Exposed Fraction"].as<double>();
  }

  if (yamlInput["Foundation"]["Building Height"].IsDefined()) {
    foundation.buildingHeight = yamlInput["Foundation"]["Building Height"].as<double>();
  } else {
    foundation.buildingHeight = 0.0;
  }

  if (yamlInput["Foundation"]["Perimeter Surface Width"].IsDefined()) {
    foundation.hasPerimeterSurface = true;
    foundation.perimeterSurfaceWidth =
        yamlInput["Foundation"]["Perimeter Surface Width"].as<double>();
  } else {
    foundation.hasPerimeterSurface = false;
    foundation.perimeterSurfaceWidth = 0.0;
  }

  // NUMERICAL SETTINGS

  // Coordinate System
  if (yamlInput["Numerical Settings"]["Coordinate System"].IsDefined()) {
    if (yamlInput["Numerical Settings"]["Coordinate System"].as<std::string>() == "CARTESIAN")
      foundation.coordinateSystem = Foundation::CS_CARTESIAN;
    else if (yamlInput["Numerical Settings"]["Coordinate System"].as<std::string>() ==
             "CYLINDRICAL")
      foundation.coordinateSystem = Foundation::CS_CYLINDRICAL;
  } else {
    foundation.coordinateSystem = Foundation::CS_CARTESIAN;
  }

  if (yamlInput["Numerical Settings"]["Two-Dimensional Approximation"].IsDefined()) {
    if (yamlInput["Numerical Settings"]["Two-Dimensional Approximation"].as<std::string>() == "AP")
      foundation.reductionStrategy = Foundation::RS_AP;
    else if (yamlInput["Numerical Settings"]["Two-Dimensional Approximation"].as<std::string>() ==
             "RR")
      foundation.reductionStrategy = Foundation::RS_RR;
    else if (yamlInput["Numerical Settings"]["Two-Dimensional Approximation"].as<std::string>() ==
             "BOUNDARY")
      foundation.reductionStrategy = Foundation::RS_BOUNDARY;
    else if (yamlInput["Numerical Settings"]["Two-Dimensional Approximation"].as<std::string>() ==
             "CUSTOM") {
      foundation.reductionStrategy = Foundation::RS_CUSTOM;
      if (yamlInput["Numerical Settings"]["Length 1"].IsDefined()) {
        foundation.twoParameters = true;
        foundation.reductionLength1 = yamlInput["Numerical Settings"]["Length 1"].as<double>();
      } else {
        foundation.twoParameters = false;
      }
      foundation.reductionLength2 = yamlInput["Numerical Settings"]["Length 2"].as<double>();
    }
  } else {
    foundation.reductionStrategy = Foundation::RS_BOUNDARY;
  }

  if (yamlInput["Numerical Settings"]["Number of Dimensions"].IsDefined())
    foundation.numberOfDimensions =
        yamlInput["Numerical Settings"]["Number of Dimensions"].as<int>();
  else
    foundation.numberOfDimensions = 2;

  if (yamlInput["Numerical Settings"]["Use Symmetry"].IsDefined())
    foundation.useSymmetry = yamlInput["Numerical Settings"]["Use Symmetry"].as<bool>();
  else
    foundation.useSymmetry = true;

  // Meshing
  if (yamlInput["Numerical Settings"]["Mesh"].IsDefined()) {

    if (yamlInput["Numerical Settings"]["Mesh"]["Minimum Cell Dimension"].IsDefined()) {
      foundation.mesh.minCellDim =
          yamlInput["Numerical Settings"]["Mesh"]["Minimum Cell Dimension"].as<double>();
    } else {
      foundation.mesh.minCellDim = 0.02;
    }

    if (yamlInput["Numerical Settings"]["Mesh"]["Maximum Near-Field Growth Coefficient"]
            .IsDefined()) {
      foundation.mesh.maxNearGrowthCoeff =
          yamlInput["Numerical Settings"]["Mesh"]["Maximum Near-Field Growth Coefficient"]
              .as<double>();
    } else {
      foundation.mesh.maxNearGrowthCoeff = 1.5;
    }

    if (yamlInput["Numerical Settings"]["Mesh"]["Maximum Deep-Field Growth Coefficient"]
            .IsDefined()) {
      foundation.mesh.maxDepthGrowthCoeff =
          yamlInput["Numerical Settings"]["Mesh"]["Maximum Deep-Field Growth Coefficient"]
              .as<double>();
    } else {
      foundation.mesh.maxDepthGrowthCoeff = 1.5;
    }

    if (yamlInput["Numerical Settings"]["Mesh"]["Maximum Interior-Field Growth Coefficient"]
            .IsDefined()) {
      foundation.mesh.maxInteriorGrowthCoeff =
          yamlInput["Numerical Settings"]["Mesh"]["Maximum Interior-Field Growth Coefficient"]
              .as<double>();
    } else {
      foundation.mesh.maxInteriorGrowthCoeff = 1.5;
    }

    if (yamlInput["Numerical Settings"]["Mesh"]["Maximum Far-Field Growth Coefficient"]
            .IsDefined()) {
      foundation.mesh.maxExteriorGrowthCoeff =
          yamlInput["Numerical Settings"]["Mesh"]["Maximum Far-Field Growth Coefficient"]
              .as<double>();
    } else {
      foundation.mesh.maxExteriorGrowthCoeff = 1.5;
    }
  } else {
    foundation.mesh.minCellDim = 0.02;
    foundation.mesh.maxNearGrowthCoeff = 1.5;
    foundation.mesh.maxDepthGrowthCoeff = 1.5;
    foundation.mesh.maxInteriorGrowthCoeff = 1.5;
    foundation.mesh.maxExteriorGrowthCoeff = 1.5;
  }

  // Simulation Control
  if (yamlInput["Numerical Settings"]["Numerical Scheme"].IsDefined()) {
    if (yamlInput["Numerical Settings"]["Numerical Scheme"].as<std::string>() == "ADE")
      foundation.numericalScheme = Foundation::NS_ADE;
    else if (yamlInput["Numerical Settings"]["Numerical Scheme"].as<std::string>() == "EXPLICIT")
      foundation.numericalScheme = Foundation::NS_EXPLICIT;
    else if (yamlInput["Numerical Settings"]["Numerical Scheme"].as<std::string>() == "ADI") {
      foundation.numericalScheme = Foundation::NS_ADI;
      if (yamlInput["Numerical Settings"]["f-ADI"].IsDefined())
        foundation.fADI = yamlInput["Numerical Settings"]["f-ADI"].as<double>();
      else
        foundation.fADI = 0.00001;
    } else if (yamlInput["Numerical Settings"]["Numerical Scheme"].as<std::string>() == "IMPLICIT")
      foundation.numericalScheme = Foundation::NS_IMPLICIT;
    else if (yamlInput["Numerical Settings"]["Numerical Scheme"].as<std::string>() ==
             "CRANK-NICOLSON")
      foundation.numericalScheme = Foundation::NS_CRANK_NICOLSON;
    else if (yamlInput["Numerical Settings"]["Numerical Scheme"].as<std::string>() ==
             "STEADY-STATE")
      foundation.numericalScheme = Foundation::NS_STEADY_STATE;
  } else {
    foundation.numericalScheme = Foundation::NS_ADI;
    foundation.fADI = 0.00001;
  }

  if (yamlInput["Numerical Settings"]["Maximum Iterations"].IsDefined()) {
    foundation.maxIterations = yamlInput["Numerical Settings"]["Maximum Iterations"].as<int>();
  } else {
    foundation.maxIterations = 100000;
  }

  if (yamlInput["Numerical Settings"]["Tolerance"].IsDefined()) {
    foundation.tolerance = yamlInput["Numerical Settings"]["Tolerance"].as<double>();
  } else {
    foundation.tolerance = 1.0e-6;
  }

  // BOUNDARIES
  if (yamlInput["Boundaries"]["Far-Field Width"].IsDefined()) {
    foundation.farFieldWidth = yamlInput["Boundaries"]["Far-Field Width"].as<double>();
  } else {
    foundation.farFieldWidth = 40;
  }

  if (yamlInput["Boundaries"]["Deep-Ground Depth"].IsDefined()) {
    foundation.deepGroundDepth = yamlInput["Boundaries"]["Deep-Ground Depth"].as<double>();
  } else {
    foundation.deepGroundDepth = 40;
  }

  if (yamlInput["Boundaries"]["Deep-Ground Boundary Condition"].IsDefined()) {
    if (yamlInput["Boundaries"]["Deep-Ground Boundary Condition"].as<std::string>() == "AUTO") {
      foundation.deepGroundBoundary = Foundation::DGB_FIXED_TEMPERATURE;
      boundaries.deepGroundBoundaryType = Boundaries::DGBT_AUTO;
    } else if (yamlInput["Boundaries"]["Deep-Ground Boundary Condition"].as<std::string>() ==
               "CONSTANT-TEMP") {
      foundation.deepGroundBoundary = Foundation::DGB_FIXED_TEMPERATURE;
      boundaries.deepGroundBoundaryType = Boundaries::DGBT_CONSTANT_TEMPERATURE;
      boundaries.deepGroundTemperature =
          yamlInput["Boundaries"]["Deep-Ground Temperature"].as<double>();
    } else if (yamlInput["Boundaries"]["Deep-Ground Boundary Condition"].as<std::string>() ==
               "ZERO-FLUX") {
      foundation.deepGroundBoundary = Foundation::DGB_ZERO_FLUX;
      boundaries.deepGroundBoundaryType = Boundaries::DGBT_ZERO_FLUX;
    }
  } else {
    foundation.deepGroundBoundary = Foundation::DGB_ZERO_FLUX;
    boundaries.deepGroundBoundaryType = Boundaries::DGBT_ZERO_FLUX;
  }

  if (yamlInput["Boundaries"]["Indoor Air Temperature Method"].IsDefined()) {
    if (yamlInput["Boundaries"]["Indoor Air Temperature Method"].as<std::string>() == "FILE") {
      boundaries.indoorTemperatureMethod = Boundaries::ITM_FILE;
      boundaries.indoorAirTemperatureFile.fileName =
          yamlInput["Boundaries"]["Indoor Air Temperature File"]["Name"].as<std::string>();
      boundaries.indoorAirTemperatureFile.firstIndex.first =
          yamlInput["Boundaries"]["Indoor Air Temperature File"]["Index"][0].as<int>();
      boundaries.indoorAirTemperatureFile.firstIndex.second =
          yamlInput["Boundaries"]["Indoor Air Temperature File"]["Index"][1].as<int>();
      boundaries.indoorAirTemperatureFile.searchDir = inputPath.parent_path();
      boundaries.indoorAirTemperatureFile.readData();
    } else if (yamlInput["Boundaries"]["Indoor Air Temperature Method"].as<std::string>() ==
               "CONSTANT") {
      boundaries.indoorTemperatureMethod = Boundaries::ITM_CONSTANT_TEMPERATURE;
      boundaries.indoorAirTemperature =
          yamlInput["Boundaries"]["Indoor Air Temperature"].as<double>();
    }
  } else {
    boundaries.indoorTemperatureMethod = Boundaries::ITM_CONSTANT_TEMPERATURE;
    boundaries.indoorAirTemperature =
        yamlInput["Boundaries"]["Indoor Air Temperature"].as<double>();
  }

  if (yamlInput["Boundaries"]["Local Boundary Layer Thickness"].IsDefined()) {
    boundaries.deltaLocal = yamlInput["Boundaries"]["Local Boundary Layer Thickness"].as<double>();
  } else {
    boundaries.deltaLocal = 370;
  }

  if (yamlInput["Boundaries"]["Local Terrain Exponent"].IsDefined()) {
    boundaries.alphaLocal = yamlInput["Boundaries"]["Local Terrain Exponent"].as<double>();
  } else {
    boundaries.alphaLocal = 0.22;
  }

  if (yamlInput["Boundaries"]["Convection Calculation Method"].IsDefined()) {
    if (yamlInput["Boundaries"]["Convection Calculation Method"].as<std::string>() == "AUTO") {
      boundaries.exteriorConvectionAlgorithm = &getDOE2ConvectionCoeff;
      boundaries.interiorConvectionAlgorithm = &getDOE2ConvectionCoeff;
    } else if (yamlInput["Boundaries"]["Convection Calculation Method"].as<std::string>() ==
               "CONSTANT") {
      double hc_ext = yamlInput["Boundaries"]["Exterior Convective Coefficient"].as<double>();
      double hc_int = yamlInput["Boundaries"]["Interior Convective Coefficient"].as<double>();
      boundaries.exteriorConvectionAlgorithm = KIVA_CONST_CONV(hc_ext);
      boundaries.interiorConvectionAlgorithm = KIVA_CONST_CONV(hc_int);
    }
  } else {
    boundaries.exteriorConvectionAlgorithm = &getDOE2ConvectionCoeff;
    boundaries.interiorConvectionAlgorithm = &getDOE2ConvectionCoeff;
  }

  if (yamlInput["Boundaries"]["Outdoor Air Temperature Method"].IsDefined()) {
    if (yamlInput["Boundaries"]["Outdoor Air Temperature Method"].as<std::string>() ==
        "WEATHER-FILE")
      boundaries.outdoorTemperatureMethod = Boundaries::OTM_WEATHER_FILE;
    else if (yamlInput["Boundaries"]["Outdoor Air Temperature Method"].as<std::string>() ==
             "CONSTANT") {
      boundaries.outdoorTemperatureMethod = Boundaries::OTM_CONSTANT_TEMPERATURE;
      boundaries.outdoorDryBulbTemperature =
          yamlInput["Boundaries"]["Outdoor Air Temperature"].as<double>();
    }
  } else {
    boundaries.outdoorTemperatureMethod = Boundaries::OTM_WEATHER_FILE;
  }

  if (yamlInput["Boundaries"]["Wall Top Boundary Condition"].IsDefined()) {
    if (yamlInput["Boundaries"]["Wall Top Boundary Condition"].as<std::string>() == "ZERO-FLUX")
      foundation.wallTopBoundary = Foundation::WTB_ZERO_FLUX;
    else if (yamlInput["Boundaries"]["Wall Top Boundary Condition"].as<std::string>() ==
             "LINEAR-DT") {
      foundation.wallTopBoundary = Foundation::WTB_LINEAR_DT;
      foundation.wallTopInteriorTemperature =
          yamlInput["Boundaries"]["Wall Top Interior Temperature"].as<double>();
      foundation.wallTopExteriorTemperature =
          yamlInput["Boundaries"]["Wall Top Exterior Temperature"].as<double>();
    }
  } else {
    foundation.wallTopBoundary = Foundation::WTB_ZERO_FLUX;
  }

  // INITIALIZATION
  if (yamlInput["Initialization"]["Initialization Method"].IsDefined()) {
    if (yamlInput["Initialization"]["Initialization Method"].as<std::string>() == "KUSUDA")
      initialization.initializationMethod = Initialization::IM_KUSUDA;
    else if (yamlInput["Initialization"]["Initialization Method"].as<std::string>() ==
             "STEADY-STATE")
      initialization.initializationMethod = Initialization::IM_STEADY_STATE;
    else if (yamlInput["Initialization"]["Initialization Method"].as<std::string>() == "CONSTANT") {
      initialization.initializationMethod = Initialization::IM_CONSTANT_TEMPERATURE;
      initialization.initialTemperature =
          yamlInput["Initialization"]["Initial Temperature"].as<double>();
    }
  } else {
    initialization.initializationMethod = Initialization::IM_STEADY_STATE;
  }

  if (yamlInput["Initialization"]["Accelerated Initialization Timestep"].IsDefined()) {
    initialization.implicitAccelTimestep =
        yamlInput["Initialization"]["Accelerated Initialization Timestep"].as<long>();
  } else {
    initialization.implicitAccelTimestep = 168;
  }

  if (yamlInput["Initialization"]["Number of Accelearted Initialization Timesteps"].IsDefined()) {
    initialization.implicitAccelPeriods =
        yamlInput["Initialization"]["Number of Accelearted Initialization Timesteps"].as<long>();
  } else {
    initialization.implicitAccelPeriods = 12;
  }

  if (yamlInput["Initialization"]["Number of Warmup Days in Initialization"].IsDefined()) {
    initialization.warmupDays =
        yamlInput["Initialization"]["Number of Warmup Days in Initialization"].as<long>();
  } else {
    initialization.warmupDays = 365;
  }

  // OUTPUT

  // CSV Reports
  if (yamlInput["Output"]["Output Report"].IsDefined()) {
    for (size_t i = 0; i < yamlInput["Output"]["Output Report"]["Reports"].size(); i++) {
      OutputVariable temp(yamlInput["Output"]["Output Report"]["Reports"][i].as<int>());
      output.outputReport.push_back(temp);
    }

    if (yamlInput["Output"]["Output Report"]["Minimum Reporting Frequency"].IsDefined()) {
      output.outputReport.minFrequency = boost::posix_time::minutes(
          yamlInput["Output"]["Output Report"]["Minimum Reporting Frequency"].as<long>());
    } else {
      output.outputReport.minFrequency = boost::posix_time::minutes(60);
    }
    output.outputReport.setOutputMap();
  }

  // Animations/Plots
  for (size_t i = 0; i < yamlInput["Output"]["Output Snapshots"].size(); i++) {
    OutputSnapshots temp;
    temp.snapshotSettings.dir =
        yamlInput["Output"]["Output Snapshots"][i]["Directory"].as<std::string>();

    if (yamlInput["Output"]["Output Snapshots"][i]["Size"].IsDefined()) {
      temp.snapshotSettings.size = yamlInput["Output"]["Output Snapshots"][i]["Size"].as<int>();
    } else {
      temp.snapshotSettings.size = 800;
    }

    if (yamlInput["Output"]["Output Snapshots"][i]["Frequency"].IsDefined()) {
      temp.snapshotSettings.frequency =
          yamlInput["Output"]["Output Snapshots"][i]["Frequency"].as<double>() * 60.0 * 60.0;
    } else {
      temp.snapshotSettings.frequency = 36 * 60.0 * 60.0;
    }

    if (yamlInput["Output"]["Output Snapshots"][i]["Mesh"].IsDefined()) {
      temp.snapshotSettings.grid = yamlInput["Output"]["Output Snapshots"][i]["Mesh"].as<bool>();
    } else {
      temp.snapshotSettings.grid = false;
    }

    if (yamlInput["Output"]["Output Snapshots"][i]["Gradients"].IsDefined()) {
      temp.snapshotSettings.gradients =
          yamlInput["Output"]["Output Snapshots"][i]["Gradients"].as<bool>();
    } else {
      temp.snapshotSettings.gradients = false;
    }

    if (yamlInput["Output"]["Output Snapshots"][i]["Contours"].IsDefined()) {
      temp.snapshotSettings.contours =
          yamlInput["Output"]["Output Snapshots"][i]["Contours"].as<bool>();
    } else {
      temp.snapshotSettings.contours = true;
    }

    if (yamlInput["Output"]["Output Snapshots"][i]["Contour Labels"].IsDefined()) {
      temp.snapshotSettings.contourLabels =
          yamlInput["Output"]["Output Snapshots"][i]["Contour Labels"].as<bool>();
    } else {
      temp.snapshotSettings.contourLabels = false;
    }

    if (yamlInput["Output"]["Output Snapshots"][i]["Axes"].IsDefined()) {
      temp.snapshotSettings.axes = yamlInput["Output"]["Output Snapshots"][i]["Axes"].as<bool>();
    } else {
      temp.snapshotSettings.axes = true;
    }

    if (yamlInput["Output"]["Output Snapshots"][i]["Timestamp"].IsDefined()) {
      temp.snapshotSettings.timestamp =
          yamlInput["Output"]["Output Snapshots"][i]["Timestamp"].as<bool>();
    } else {
      temp.snapshotSettings.timestamp = true;
    }

    if (yamlInput["Output"]["Output Snapshots"][i]["Plot Type"].IsDefined()) {
      if (yamlInput["Output"]["Output Snapshots"][i]["Plot Type"].as<std::string>() ==
          "TEMPERATURE")
        temp.snapshotSettings.plotType = SnapshotSettings::P_TEMP;
      else if (yamlInput["Output"]["Output Snapshots"][i]["Plot Type"].as<std::string>() ==
               "HEAT-FLUX")
        temp.snapshotSettings.plotType = SnapshotSettings::P_FLUX;
    } else
      temp.snapshotSettings.plotType = SnapshotSettings::P_TEMP;

    if (yamlInput["Output"]["Output Snapshots"][i]["Flux Direction"].IsDefined()) {
      if (yamlInput["Output"]["Output Snapshots"][i]["Flux Direction"].as<std::string>() == "MAG")
        temp.snapshotSettings.fluxDir = SnapshotSettings::D_M;
      else if (yamlInput["Output"]["Output Snapshots"][i]["Flux Direction"].as<std::string>() ==
               "X")
        temp.snapshotSettings.fluxDir = SnapshotSettings::D_X;
      else if (yamlInput["Output"]["Output Snapshots"][i]["Flux Direction"].as<std::string>() ==
               "Y")
        temp.snapshotSettings.fluxDir = SnapshotSettings::D_Y;
      else if (yamlInput["Output"]["Output Snapshots"][i]["Flux Direction"].as<std::string>() ==
               "Z")
        temp.snapshotSettings.fluxDir = SnapshotSettings::D_Z;
    } else
      temp.snapshotSettings.fluxDir = SnapshotSettings::D_M;

    if (yamlInput["Output"]["Output Snapshots"][i]["Color Scheme"].IsDefined()) {
      if (yamlInput["Output"]["Output Snapshots"][i]["Color Scheme"].as<std::string>() == "CMR")
        temp.snapshotSettings.colorScheme = SnapshotSettings::C_CMR;
      else if (yamlInput["Output"]["Output Snapshots"][i]["Color Scheme"].as<std::string>() ==
               "JET")
        temp.snapshotSettings.colorScheme = SnapshotSettings::C_JET;
      else if (yamlInput["Output"]["Output Snapshots"][i]["Color Scheme"].as<std::string>() ==
               "NONE")
        temp.snapshotSettings.colorScheme = SnapshotSettings::C_NONE;
    } else
      temp.snapshotSettings.colorScheme = SnapshotSettings::C_CMR;

    if (yamlInput["Output"]["Output Snapshots"][i]["File Format"].IsDefined()) {
      if (yamlInput["Output"]["Output Snapshots"][i]["File Format"].as<std::string>() == "PNG")
        temp.snapshotSettings.format = SnapshotSettings::F_PNG;
      else if (yamlInput["Output"]["Output Snapshots"][i]["File Format"].as<std::string>() == "TEX")
        temp.snapshotSettings.format = SnapshotSettings::F_TEX;
    } else
      temp.snapshotSettings.format = SnapshotSettings::F_PNG;

    if (yamlInput["Output"]["Output Snapshots"][i]["Unit System"].IsDefined()) {
      if (yamlInput["Output"]["Output Snapshots"][i]["Unit System"].as<std::string>() == "IP")
        temp.snapshotSettings.outputUnits = SnapshotSettings::IP;
      else if (yamlInput["Output"]["Output Snapshots"][i]["Unit System"].as<std::string>() == "SI")
        temp.snapshotSettings.outputUnits = SnapshotSettings::SI;
    } else
      temp.snapshotSettings.outputUnits = SnapshotSettings::SI;

    if (yamlInput["Output"]["Output Snapshots"][i]["Output Range"].IsDefined()) {
      temp.snapshotSettings.minValue =
          yamlInput["Output"]["Output Snapshots"][i]["Output Range"][0].as<double>();
      temp.snapshotSettings.maxValue =
          yamlInput["Output"]["Output Snapshots"][i]["Output Range"][1].as<double>();
    } else {
      temp.snapshotSettings.minValue = -20;
      temp.snapshotSettings.maxValue = 40;
    }

    if (yamlInput["Output"]["Output Snapshots"][i]["Number of Contours"].IsDefined()) {
      temp.snapshotSettings.numberOfContours =
          yamlInput["Output"]["Output Snapshots"][i]["Number of Contours"].as<int>();
    } else {
      temp.snapshotSettings.numberOfContours = 13;
    }

    if (yamlInput["Output"]["Output Snapshots"][i]["Contour Color"].IsDefined()) {
      temp.snapshotSettings.contourColor =
          yamlInput["Output"]["Output Snapshots"][i]["Contour Color"].as<std::string>();
    } else {
      temp.snapshotSettings.contourColor = "H";
    }

    if (yamlInput["Output"]["Output Snapshots"][i]["Start Date"].IsDefined()) {
      temp.startDate = boost::gregorian::from_string(
          yamlInput["Output"]["Output Snapshots"][i]["Start Date"].as<std::string>());
      temp.startDateSet = true;
    } else
      temp.startDateSet = false;

    if (yamlInput["Output"]["Output Snapshots"][i]["End Date"].IsDefined()) {
      temp.endDate = boost::gregorian::from_string(
          yamlInput["Output"]["Output Snapshots"][i]["End Date"].as<std::string>());
      temp.endDateSet = true;
    } else
      temp.endDateSet = false;

    if (yamlInput["Output"]["Output Snapshots"][i]["X Range"].IsDefined()) {
      temp.snapshotSettings.xRange.first =
          yamlInput["Output"]["Output Snapshots"][i]["X Range"][0].as<double>();
      temp.snapshotSettings.xRange.second =
          yamlInput["Output"]["Output Snapshots"][i]["X Range"][1].as<double>();
      temp.xRangeSet = true;
    } else
      temp.xRangeSet = false;

    if (yamlInput["Output"]["Output Snapshots"][i]["Y Range"].IsDefined()) {
      temp.snapshotSettings.yRange.first =
          yamlInput["Output"]["Output Snapshots"][i]["Y Range"][0].as<double>();
      temp.snapshotSettings.yRange.second =
          yamlInput["Output"]["Output Snapshots"][i]["Y Range"][1].as<double>();
      temp.yRangeSet = true;
    } else
      temp.yRangeSet = false;

    if (yamlInput["Output"]["Output Snapshots"][i]["Z Range"].IsDefined()) {

      temp.snapshotSettings.zRange.first =
          yamlInput["Output"]["Output Snapshots"][i]["Z Range"][0].as<double>();
      temp.snapshotSettings.zRange.second =
          yamlInput["Output"]["Output Snapshots"][i]["Z Range"][1].as<double>();
      temp.zRangeSet = true;
    } else
      temp.zRangeSet = false;

    output.outputSnapshots.push_back(temp);
  }

  // Full Input
  input.simulationControl = simulationControl;
  input.foundation = foundation;
  input.boundaries = boundaries;
  input.initialization = initialization;
  input.output = output;

  return input;
}
