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
#include <fstream>

using namespace std;

Input inputParser(string inputFile)
{

	Input input;
	SimulationControl simulationControl;
	Foundation foundation1;

	YAML::Node yamlInput = YAML::LoadFile("in.yaml");


	simulationControl.weatherFile =
			yamlInput["Simulation Control"]["weatherFile"].as<string>();
	simulationControl.startDate =
			from_string(yamlInput["Simulation Control"]["startDate"].as<string>());
	simulationControl.endDate =
			from_string(yamlInput["Simulation Control"]["endDate"].as<string>());
	simulationControl.timestep =
			minutes(yamlInput["Simulation Control"]["timeStep"].as<long>());

	// Materials
	map<string, Material> materials;

	for(YAML::const_iterator it=yamlInput["Materials"].begin();it!=yamlInput["Materials"].end();++it)
	{
		Material tempMaterial;
		tempMaterial.conductivity = it->second["k"].as<double>();
		tempMaterial.density = it->second["rho"].as<double>();
		tempMaterial.specificHeat = it->second["cp"].as<double>();

		materials.insert(std::pair<string,Material>(it->first.as<string>(),tempMaterial));
	}

	// Foundation

	// Soil
	foundation1.soil = materials[yamlInput["Foundation"]["soil"].as<string>()];
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
			tempLayer.material = materials[yamlInput["Foundation"]["slab"]["layers"][i]["material"].as<string>()];

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
			tempLayer.material = materials[yamlInput["Foundation"]["slab"]["layers"][i]["material"].as<string>()];

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
		foundation1.interiorHorizontalInsulation.layer.material = materials[yamlInput["Foundation"]["interiorHorizontalInsulation"]["material"].as<string>()];
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
		foundation1.interiorVerticalInsulation.layer.material = materials[yamlInput["Foundation"]["interiorVerticalInsulation"]["material"].as<string>()];
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
		foundation1.exteriorHorizontalInsulation.layer.material = materials[yamlInput["Foundation"]["exteriorHorizontalInsulation"]["material"].as<string>()];
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
		foundation1.exteriorVerticalInsulation.layer.material = materials[yamlInput["Foundation"]["exteriorVerticalInsulation"]["material"].as<string>()];
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

	if (yamlInput["Foundation"]["deepGroundBoundary"].as<string>() == "AUTO")
		foundation1.deepGroundBoundary = Foundation::DGB_AUTO;
	else if (yamlInput["Foundation"]["deepGroundBoundary"].as<string>() == "CONSTANT-TEMP")
	{
		foundation1.deepGroundBoundary = Foundation::DGB_CONSTANT_TEMPERATURE;
		foundation1.deepGroundTemperature = yamlInput["Foundation"]["deepGroundTemperature"].as<double>();
	}


	foundation1.indoorAirTemperature = yamlInput["Foundation"]["indoorAirTemperature"].as<double>();


	// Geometry
	foundation1.area = yamlInput["Foundation"]["area"].as<double>();
	foundation1.perimeter = yamlInput["Foundation"]["perimeter"].as<double>();

	// Meshing
	if  (yamlInput["Foundation"]["mesh"].IsDefined())
	{
		foundation1.mesh.minCellDim = yamlInput["Foundation"]["mesh"]["minCellDim"].as<double>();
		foundation1.mesh.maxDepthGrowthCoeff = yamlInput["Foundation"]["mesh"]["maxDepthGrowthCoeff"].as<double>();
		foundation1.mesh.maxInteriorGrowthCoeff = yamlInput["Foundation"]["mesh"]["maxInteriorGrowthCoeff"].as<double>();
		foundation1.mesh.maxRadialGrowthCoeff = yamlInput["Foundation"]["mesh"]["maxExteriorGrowthCoeff"].as<double>();
	}
	else
	{
		foundation1.mesh.minCellDim = 0.05;
		foundation1.mesh.maxDepthGrowthCoeff = 1.25;
		foundation1.mesh.maxInteriorGrowthCoeff = 1.25;
		foundation1.mesh.maxRadialGrowthCoeff = 1.25;
	}

	// Simulation Control
	if  (yamlInput["Foundation"]["numericalScheme"].IsDefined())
	{
		if (yamlInput["Foundation"]["numericalScheme"].as<string>() == "ADE")
			foundation1.numericalScheme = Foundation::NS_ADE;
		else if (yamlInput["Foundation"]["numericalScheme"].as<string>() == "EXPLICIT")
			foundation1.numericalScheme = Foundation::NS_EXPLICIT;
		else if (yamlInput["Foundation"]["numericalScheme"].as<string>() == "IMPLICIT")
			foundation1.numericalScheme = Foundation::NS_IMPLICIT;
		else if (yamlInput["Foundation"]["numericalScheme"].as<string>() == "CRANK-NICOLSON")
			foundation1.numericalScheme = Foundation::NS_CRANK_NICOLSON;
		else if (yamlInput["Foundation"]["numericalScheme"].as<string>() == "STEADY-STATE")
			foundation1.numericalScheme = Foundation::NS_STEADY_STATE;
	}
	else
	{
		foundation1.numericalScheme = Foundation::NS_ADE;
	}

	if  (yamlInput["Foundation"]["initializationMethod"].IsDefined())
	{
		if (yamlInput["Foundation"]["initializationMethod"].as<string>() == "KUSUDA")
			foundation1.initializationMethod = Foundation::IM_KUSUDA;
		else if (yamlInput["Foundation"]["initializationMethod"].as<string>() == "IMPLICIT-ACCEL")
		{
			foundation1.initializationMethod = Foundation::IM_IMPLICIT_ACCELERATION;
			foundation1.implicitAccelTimestep = yamlInput["Foundation"]["implicitAccelTimestep"].as<long>();
			foundation1.implicitAccelPeriods = yamlInput["Foundation"]["implicitAccelPeriods"].as<long>();
		}
		else if (yamlInput["Foundation"]["initializationMethod"].as<string>() == "STEADY-STATE")
			foundation1.initializationMethod = Foundation::IM_STEADY_STATE;
		else if (yamlInput["Foundation"]["initializationMethod"].as<string>() == "CONSTANT")
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
		if (yamlInput["Foundation"]["convectionCalculationMethod"].as<string>() == "AUTO")
			foundation1.convectionCalculationMethod = Foundation::CCM_AUTO;
		else if (yamlInput["Foundation"]["convectionCalculationMethod"].as<string>() == "CONSTANT")
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
		if (yamlInput["Foundation"]["outdoorTemperatureMethod"].as<string>() == "WEATHER-FILE")
			foundation1.outdoorTemperatureMethod = Foundation::OTM_WEATHER_FILE;
		else if (yamlInput["Foundation"]["outdoorTemperatureMethod"].as<string>() == "CONSTANT")
		{
			foundation1.outdoorTemperatureMethod = Foundation::OTM_CONSTANT_TEMPERATURE;
			foundation1.outdoorDryBulbTemperature = yamlInput["Foundation"]["outdoorDryBulbTemperature"].as<double>();
		}
	}
	else
	{
		foundation1.outdoorTemperatureMethod = Foundation::OTM_WEATHER_FILE;
	}

	// Output

	foundation1.outputAnimation.name = yamlInput["Foundation"]["outputAnimation"]["name"].as<string>();
	foundation1.outputAnimation.frequency = hours(yamlInput["Foundation"]["outputAnimation"]["frequency"].as<long>());
	foundation1.outputAnimation.grid = yamlInput["Foundation"]["outputAnimation"]["grid"].as<bool>();
	foundation1.outputAnimation.gradients = yamlInput["Foundation"]["outputAnimation"]["gradients"].as<bool>();
	foundation1.outputAnimation.contours = yamlInput["Foundation"]["outputAnimation"]["contours"].as<bool>();
	foundation1.outputAnimation.size = yamlInput["Foundation"]["outputAnimation"]["size"].as<int>();

	// Full Input
	input.simulationControl = simulationControl;
	input.foundations.push_back(foundation1);


	return input;

}


