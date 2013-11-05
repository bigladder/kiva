/* InputParser.c++ is part of Kiva (Written by Neal Kruis)
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


	simulationControl.weatherFile =
			yamlInput["Simulation Control"]["weatherFile"].as<std::string>();
	simulationControl.startDate =
			boost::gregorian::from_string(yamlInput["Simulation Control"]["startDate"].as<std::string>());
	simulationControl.endDate =
			boost::gregorian::from_string(yamlInput["Simulation Control"]["endDate"].as<std::string>());
	simulationControl.timestep =
			boost::posix_time::minutes(yamlInput["Simulation Control"]["timeStep"].as<long>());

	// Materials
	std::map<std::string, Material> materials;

	for(YAML::const_iterator it=yamlInput["Materials"].begin();it!=yamlInput["Materials"].end();++it)
	{
		Material tempMaterial;
		tempMaterial.conductivity = it->second["k"].as<double>();
		tempMaterial.density = it->second["rho"].as<double>();
		tempMaterial.specificHeat = it->second["cp"].as<double>();

		materials.insert(std::pair<std::string,Material>(it->first.as<std::string>(),tempMaterial));
	}

	// Foundation

	// Soil
	foundation1.soil = materials[yamlInput["Foundation"]["soil"].as<std::string>()];
    foundation1.soilAbsorptivity = yamlInput["Foundation"]["soilAbsorptivity"].as<double>();
    foundation1.soilEmissivity = yamlInput["Foundation"]["soilEmissivity"].as<double>();

    // Local wind speed characteristics
    foundation1.vegetationHeight = yamlInput["Foundation"]["vegetationHeight"].as<double>();
    foundation1.deltaLocal = yamlInput["Foundation"]["deltaLocal"].as<double>();
    foundation1.alphaLocal = yamlInput["Foundation"]["alphaLocal"].as<double>();


	// Slab
	if  (yamlInput["Foundation"]["slab"].IsDefined())
	{
		foundation1.hasSlab = true;

		for (size_t i=0;i<yamlInput["Foundation"]["slab"]["layers"].size();i++)
		{

			Layer tempLayer;
			tempLayer.thickness = yamlInput["Foundation"]["slab"]["layers"][i]["thickness"].as<double>();
			tempLayer.material = materials[yamlInput["Foundation"]["slab"]["layers"][i]["material"].as<std::string>()];

			foundation1.slab.layers.push_back(tempLayer);

		}

		foundation1.slab.emissivity = yamlInput["Foundation"]["slab"]["emissivity"].as<double>();

	}
	else
	{
		foundation1.hasSlab = false;
	}

	// Wall
	if  (yamlInput["Foundation"]["wall"].IsDefined())
	{
		foundation1.hasWall = true;

		for (size_t i=0;i<yamlInput["Foundation"]["wall"]["layers"].size();i++)
		{

			Layer tempLayer;
			tempLayer.thickness = yamlInput["Foundation"]["wall"]["layers"][i]["thickness"].as<double>();
			tempLayer.material = materials[yamlInput["Foundation"]["wall"]["layers"][i]["material"].as<std::string>()];

			foundation1.wall.layers.push_back(tempLayer);

		}

		foundation1.wall.depth = yamlInput["Foundation"]["wall"]["depth"].as<double>();
		foundation1.wall.height = yamlInput["Foundation"]["wall"]["height"].as<double>();
		foundation1.wall.interiorEmissivity = yamlInput["Foundation"]["wall"]["interiorEmissivity"].as<double>();
		foundation1.wall.exteriorEmissivity = yamlInput["Foundation"]["wall"]["exteriorEmissivity"].as<double>();
		foundation1.wall.exteriorAbsorptivity = yamlInput["Foundation"]["wall"]["exteriorAbsorptivity"].as<double>();
	}
	else
	{
		foundation1.hasWall = false;
	}

	// Interior Horizontal Insulation
	if  (yamlInput["Foundation"]["interiorHorizontalInsulation"].IsDefined())
	{
		foundation1.hasInteriorHorizontalInsulation = true;

		foundation1.interiorHorizontalInsulation.layer.thickness = yamlInput["Foundation"]["interiorHorizontalInsulation"]["thickness"].as<double>();
		foundation1.interiorHorizontalInsulation.layer.material = materials[yamlInput["Foundation"]["interiorHorizontalInsulation"]["material"].as<std::string>()];
		foundation1.interiorHorizontalInsulation.depth = yamlInput["Foundation"]["interiorHorizontalInsulation"]["depth"].as<double>();
		foundation1.interiorHorizontalInsulation.width = yamlInput["Foundation"]["interiorHorizontalInsulation"]["width"].as<double>();

	}
	else
	{
		foundation1.hasInteriorHorizontalInsulation = false;
	}

	// Interior Vertical Insulation
	if  (yamlInput["Foundation"]["interiorVerticalInsulation"].IsDefined())
	{
		foundation1.hasInteriorVerticalInsulation = true;

		foundation1.interiorVerticalInsulation.layer.thickness = yamlInput["Foundation"]["interiorVerticalInsulation"]["thickness"].as<double>();
		foundation1.interiorVerticalInsulation.layer.material = materials[yamlInput["Foundation"]["interiorVerticalInsulation"]["material"].as<std::string>()];
		foundation1.interiorVerticalInsulation.depth = yamlInput["Foundation"]["interiorVerticalInsulation"]["depth"].as<double>();

	}
	else
	{
		foundation1.hasInteriorVerticalInsulation = false;
	}

	// Exterior Horizontal Insulation
	if  (yamlInput["Foundation"]["exteriorHorizontalInsulation"].IsDefined())
	{
		foundation1.hasExteriorHorizontalInsulation = true;

		foundation1.exteriorHorizontalInsulation.layer.thickness = yamlInput["Foundation"]["exteriorHorizontalInsulation"]["thickness"].as<double>();
		foundation1.exteriorHorizontalInsulation.layer.material = materials[yamlInput["Foundation"]["exteriorHorizontalInsulation"]["material"].as<std::string>()];
		foundation1.exteriorHorizontalInsulation.depth = yamlInput["Foundation"]["exteriorHorizontalInsulation"]["depth"].as<double>();
		foundation1.exteriorHorizontalInsulation.width = yamlInput["Foundation"]["exteriorHorizontalInsulation"]["width"].as<double>();

	}
	else
	{
		foundation1.hasExteriorHorizontalInsulation = false;
	}

	// Exterior Vertical Insulation
	if  (yamlInput["Foundation"]["exteriorVerticalInsulation"].IsDefined())
	{
		foundation1.hasExteriorVerticalInsulation = true;

		foundation1.exteriorVerticalInsulation.layer.thickness = yamlInput["Foundation"]["exteriorVerticalInsulation"]["thickness"].as<double>();
		foundation1.exteriorVerticalInsulation.layer.material = materials[yamlInput["Foundation"]["exteriorVerticalInsulation"]["material"].as<std::string>()];
		foundation1.exteriorVerticalInsulation.depth = yamlInput["Foundation"]["exteriorVerticalInsulation"]["depth"].as<double>();

	}
	else
	{
		foundation1.hasExteriorVerticalInsulation = false;
	}


	// Site

	foundation1.excavationDepth = yamlInput["Foundation"]["excavationDepth"].as<double>();
	foundation1.farFieldWidth = yamlInput["Foundation"]["farFieldWidth"].as<double>();
	foundation1.deepGroundDepth = yamlInput["Foundation"]["deepGroundDepth"].as<double>();

	if (yamlInput["Foundation"]["deepGroundBoundary"].as<std::string>() == "AUTO")
	{
		foundation1.deepGroundBoundary = Foundation::DGB_AUTO;
		WeatherData tempWeather(simulationControl.weatherFile);
		foundation1.deepGroundTemperature = tempWeather.dryBulbTemp.getAverage();
	}
	else if (yamlInput["Foundation"]["deepGroundBoundary"].as<std::string>() == "CONSTANT-TEMP")
	{
		foundation1.deepGroundBoundary = Foundation::DGB_CONSTANT_TEMPERATURE;
		foundation1.deepGroundTemperature = yamlInput["Foundation"]["deepGroundTemperature"].as<double>();
	}
	else if (yamlInput["Foundation"]["deepGroundBoundary"].as<std::string>() == "ZERO-FLUX")
	{
		foundation1.deepGroundBoundary = Foundation::DGB_ZERO_FLUX;
	}

	foundation1.indoorAirTemperature = yamlInput["Foundation"]["indoorAirTemperature"].as<double>();


	// Geometry
	if (yamlInput["Foundation"]["coordinateSystem"].as<std::string>() == "2DAXIAL")
		foundation1.coordinateSystem = Foundation::CS_2DAXIAL;
	else if (yamlInput["Foundation"]["coordinateSystem"].as<std::string>() == "2DLINEAR")
		foundation1.coordinateSystem = Foundation::CS_2DLINEAR;
	else if (yamlInput["Foundation"]["coordinateSystem"].as<std::string>() == "3D")
		foundation1.coordinateSystem = Foundation::CS_3D;

	for (size_t i=0;i<yamlInput["Foundation"]["polygon"].size();i++)
	{
		foundation1.polygon.outer().push_back(Point(
				yamlInput["Foundation"]["polygon"][i][0].as<double>(),
				yamlInput["Foundation"]["polygon"][i][1].as<double>()));
	}

	// Meshing
	if  (yamlInput["Foundation"]["mesh"].IsDefined())
	{
		foundation1.mesh.minCellDim = yamlInput["Foundation"]["mesh"]["minCellDim"].as<double>();
		foundation1.mesh.maxDepthGrowthCoeff = yamlInput["Foundation"]["mesh"]["maxDepthGrowthCoeff"].as<double>();
		foundation1.mesh.maxInteriorGrowthCoeff = yamlInput["Foundation"]["mesh"]["maxInteriorGrowthCoeff"].as<double>();
		foundation1.mesh.maxExteriorGrowthCoeff = yamlInput["Foundation"]["mesh"]["maxExteriorGrowthCoeff"].as<double>();
	}
	else
	{
		foundation1.mesh.minCellDim = 0.05;
		foundation1.mesh.maxDepthGrowthCoeff = 1.25;
		foundation1.mesh.maxInteriorGrowthCoeff = 1.25;
		foundation1.mesh.maxExteriorGrowthCoeff = 1.25;
	}

	// Simulation Control
	if  (yamlInput["Foundation"]["numericalScheme"].IsDefined())
	{
		if (yamlInput["Foundation"]["numericalScheme"].as<std::string>() == "ADE")
			foundation1.numericalScheme = Foundation::NS_ADE;
		else if (yamlInput["Foundation"]["numericalScheme"].as<std::string>() == "EXPLICIT")
			foundation1.numericalScheme = Foundation::NS_EXPLICIT;
		else if (yamlInput["Foundation"]["numericalScheme"].as<std::string>() == "ADI")
		{
			foundation1.numericalScheme = Foundation::NS_ADI;
			if  (yamlInput["Foundation"]["fADI"].IsDefined())
				foundation1.fADI = yamlInput["Foundation"]["fADI"].as<double>();
			else
				foundation1.fADI = 0.01;
		}
		else if (yamlInput["Foundation"]["numericalScheme"].as<std::string>() == "IMPLICIT")
			foundation1.numericalScheme = Foundation::NS_IMPLICIT;
		else if (yamlInput["Foundation"]["numericalScheme"].as<std::string>() == "CRANK-NICOLSON")
			foundation1.numericalScheme = Foundation::NS_CRANK_NICOLSON;
		else if (yamlInput["Foundation"]["numericalScheme"].as<std::string>() == "STEADY-STATE")
			foundation1.numericalScheme = Foundation::NS_STEADY_STATE;
	}
	else
	{
		foundation1.numericalScheme = Foundation::NS_ADE;
	}


	if  (yamlInput["Foundation"]["initializationMethod"].IsDefined())
	{
		if (yamlInput["Foundation"]["initializationMethod"].as<std::string>() == "KUSUDA")
			foundation1.initializationMethod = Foundation::IM_KUSUDA;
		else if (yamlInput["Foundation"]["initializationMethod"].as<std::string>() == "IMPLICIT-ACCEL")
		{
			foundation1.initializationMethod = Foundation::IM_IMPLICIT_ACCELERATION;
			foundation1.implicitAccelTimestep = yamlInput["Foundation"]["implicitAccelTimestep"].as<long>();
			foundation1.implicitAccelPeriods = yamlInput["Foundation"]["implicitAccelPeriods"].as<long>();
		}
		else if (yamlInput["Foundation"]["initializationMethod"].as<std::string>() == "STEADY-STATE")
			foundation1.initializationMethod = Foundation::IM_STEADY_STATE;
		else if (yamlInput["Foundation"]["initializationMethod"].as<std::string>() == "CONSTANT")
		{
			foundation1.initializationMethod = Foundation::IM_CONSTANT_TEMPERATURE;
			foundation1.initialTemperature = yamlInput["Foundation"]["initialTemperature"].as<double>();
		}
	}
	else
	{
		foundation1.initializationMethod = Foundation::IM_KUSUDA;
	}

	if  (yamlInput["Foundation"]["convectionCalculationMethod"].IsDefined())
	{
		if (yamlInput["Foundation"]["convectionCalculationMethod"].as<std::string>() == "AUTO")
			foundation1.convectionCalculationMethod = Foundation::CCM_AUTO;
		else if (yamlInput["Foundation"]["convectionCalculationMethod"].as<std::string>() == "CONSTANT")
		{
			foundation1.convectionCalculationMethod = Foundation::CCM_CONSTANT_COEFFICIENT;
			foundation1.interiorConvectiveCoefficient = yamlInput["Foundation"]["interiorConvectiveCoefficient"].as<double>();
			foundation1.exteriorConvectiveCoefficient = yamlInput["Foundation"]["exteriorConvectiveCoefficient"].as<double>();
		}
	}
	else
	{
		foundation1.convectionCalculationMethod = Foundation::CCM_AUTO;
	}

	if  (yamlInput["Foundation"]["outdoorTemperatureMethod"].IsDefined())
	{
		if (yamlInput["Foundation"]["outdoorTemperatureMethod"].as<std::string>() == "WEATHER-FILE")
			foundation1.outdoorTemperatureMethod = Foundation::OTM_WEATHER_FILE;
		else if (yamlInput["Foundation"]["outdoorTemperatureMethod"].as<std::string>() == "CONSTANT")
		{
			foundation1.outdoorTemperatureMethod = Foundation::OTM_CONSTANT_TEMPERATURE;
			foundation1.outdoorDryBulbTemperature = yamlInput["Foundation"]["outdoorDryBulbTemperature"].as<double>();
		}
	}
	else
	{
		foundation1.outdoorTemperatureMethod = Foundation::OTM_WEATHER_FILE;
	}

	if  (yamlInput["Foundation"]["wallTopBoundary"].IsDefined())
	{
		if (yamlInput["Foundation"]["wallTopBoundary"].as<std::string>() == "ZERO-FLUX")
			foundation1.wallTopBoundary = Foundation::WTB_ZERO_FLUX;
		else if (yamlInput["Foundation"]["wallTopBoundary"].as<std::string>() == "LINEAR_DT")
		{
			foundation1.wallTopBoundary = Foundation::WTB_LINEAR_DT;
			foundation1.wallTopTemperatureDifference = yamlInput["Foundation"]["wallTopTemperatureDifference"].as<double>();
		}
	}
	else
	{
		foundation1.wallTopBoundary = Foundation::WTB_ZERO_FLUX;
	}

	// Output

	for(size_t i=0;i<yamlInput["Foundation"]["outputAnimations"].size();i++)
	{
		OutputAnimation temp;
		temp.name = yamlInput["Foundation"]["outputAnimations"][i]["name"].as<std::string>();
		temp.frequency = boost::posix_time::hours(yamlInput["Foundation"]["outputAnimations"][i]["frequency"].as<long>());
		temp.grid = yamlInput["Foundation"]["outputAnimations"][i]["grid"].as<bool>();
		temp.gradients = yamlInput["Foundation"]["outputAnimations"][i]["gradients"].as<bool>();
		temp.contours = yamlInput["Foundation"]["outputAnimations"][i]["contours"].as<bool>();
		temp.size = yamlInput["Foundation"]["outputAnimations"][i]["size"].as<int>();

		if (yamlInput["Foundation"]["outputAnimations"][i]["startDate"].IsDefined())
		{
			temp.startDate = boost::gregorian::from_string(yamlInput["Foundation"]["outputAnimations"][i]["startDate"].as<std::string>());
			temp.startDateSet = true;
		}
		else
			temp.startDateSet = false;

		if (yamlInput["Foundation"]["outputAnimations"][i]["endDate"].IsDefined())
		{
			temp.endDate = boost::gregorian::from_string(yamlInput["Foundation"]["outputAnimations"][i]["endDate"].as<std::string>());
			temp.endDateSet = true;
		}
		else
			temp.endDateSet = false;

		if (yamlInput["Foundation"]["outputAnimations"][i]["xRange"].IsDefined())
		{
			temp.xRange.first = yamlInput["Foundation"]["outputAnimations"][i]["xRange"][0].as<double>();
			temp.xRange.second = yamlInput["Foundation"]["outputAnimations"][i]["xRange"][1].as<double>();
			temp.xRangeSet = true;
		}
		else
			temp.xRangeSet = false;

		if (yamlInput["Foundation"]["outputAnimations"][i]["yRange"].IsDefined())
		{
			temp.yRange.first = yamlInput["Foundation"]["outputAnimations"][i]["yRange"][0].as<double>();
			temp.yRange.second = yamlInput["Foundation"]["outputAnimations"][i]["yRange"][1].as<double>();
			temp.yRangeSet = true;
		}
		else
			temp.yRangeSet = false;

		if (yamlInput["Foundation"]["outputAnimations"][i]["zRange"].IsDefined())
		{

			temp.zRange.first = yamlInput["Foundation"]["outputAnimations"][i]["zRange"][0].as<double>();
			temp.zRange.second = yamlInput["Foundation"]["outputAnimations"][i]["zRange"][1].as<double>();
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


