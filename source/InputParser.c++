/*
 * inputParser.cpp
 *
 *  Created on: Oct 25, 2012
 *      Author: nkruis
 */

#include "Input.h"

Input inputParser()
{

	Input input;

	// Simulation Control

	SimulationControl simulationControl;
	simulationControl.weatherFile = "/Applications/EnergyPlus-7-2-0/WeatherData/USA_IL_Chicago-OHare.Intl.AP.725300_TMY3.epw";
	date dayStart(2013,Jan,1);
	date dayEnd(2014,Dec,31);
	simulationControl.startDate = dayStart;
	simulationControl.endDate = dayEnd;
	simulationControl.timestep = hours(1) + minutes(0);

	// Materials

	// Soils (from wolfram alpha)
	Material wetSandySoil;  // 3.349e-7  Extreme high!
	wetSandySoil.k = 0.586;  // [W/m-K]
	wetSandySoil.rho = 1750.0;  // [kg/m3]
	wetSandySoil.cp = 1000.0;  // [J/kg-K]

	Material wetLoamSoil;  // 2.488e-7
	wetLoamSoil.k = 0.418;
	wetLoamSoil.rho = 1600.0;
	wetLoamSoil.cp = 1050.0;

	Material wetClaySoil;  // 2.617e-7
	wetClaySoil.k = 1.51;
	wetClaySoil.rho = 1500.0;
	wetClaySoil.cp = 2930.0;

	Material drySandySoil;  // 2.013e-7  Extreme low!
	drySandySoil.k = 0.264;
	drySandySoil.rho = 1650.0;
	drySandySoil.cp = 795.0;

	Material dryLoamSoil;  // 2.399e-7
	dryLoamSoil.k = 0.251;
	dryLoamSoil.rho = 1250.0;
	dryLoamSoil.cp = 837.0;

	// Other
	Material concrete;
	concrete.k = 1.98;
	concrete.rho = 1900.0;
	concrete.cp = 665.0;

	Material xps;
	xps.k = 0.03;
	xps.rho = 28.0;
	xps.cp = 1450;

	Material carpet;
	carpet.k = 0.08;
	carpet.rho = 28.0;
	carpet.cp = 1450;

	Material s140a;
	s140a.k = 1.9;
	s140a.rho = 1490;
	s140a.cp = 1800;

	Material krartiConcrete;
	krartiConcrete.k = 1.5;
	krartiConcrete.rho = 1490;
	krartiConcrete.cp = 1800;

	Material krartiSoil;
	krartiSoil.k = 1.0;
	krartiSoil.rho = 1490;
	krartiSoil.cp = 1800;

	// Layers
	Layer concrete05;
	concrete05.d = 0.05;  // [m]
	concrete05.material = concrete;

	Layer concrete10;
	concrete10.d = 0.1;  // [m]
	concrete10.material = concrete;

	Layer concrete6_5;
	concrete6_5.d = 0.1651;  // [m]
	concrete6_5.material = concrete;

	Layer s140slab;
	s140slab.d = 0.1;  // [m]
	s140slab.material = s140a;

	Layer xps05;
	xps05.d = 0.05;
	xps05.material = xps;

	Layer xps10;
	xps10.d = 0.10;
	xps10.material = xps;

	Layer xps20;
	xps20.d = 0.20;
	xps20.material = xps;

	Layer xps025;
	xps025.d = 0.025;
	xps025.material = xps;

	Layer xps125;
	xps125.d = 0.125;
	xps125.material = xps;

	Layer xps25;
	xps25.d = 0.25;
	xps25.material = xps;

	Layer xps15;
	xps15.d = 0.15;
	xps15.material = xps;
	Layer xps263;

	xps263.d = 0.263;
	xps263.material = xps;

	Layer xps239;
	xps239.d = 0.239;
	xps239.material = xps;

	Layer carpet02;
	carpet02.d = 0.02;
	carpet02.material = carpet;

	Layer s140wall;
	s140wall.d = 0.24;
	s140wall.material = s140a;

	Layer krartiPile;
	krartiPile.d = 4.0*0.3048;
	krartiPile.material = krartiConcrete;

	Layer null;
	null.d = 0.0;
	null.material = xps;


	// Slab
	Slab slab;
	slab.layers.push_back(s140slab);
	//slab.layers.push_back(concrete10);
	//slab.layers.push_back(carpet02);

	// Wall
	Wall wall;
	wall.depth = 0.61;  // [m]
	wall.height = 0.61;  // [m]
	wall.layers.push_back(s140wall);
	//wall.layers.push_back(xps05);

	// Interior Horizontal Insulation
	HorizontalInsulation ihi;
	ihi.depth = 0.225;  // [m]
	ihi.width = 1.22;  // [m]
	ihi.layer = xps125;

	// Interior Vertical Insulation
	VerticalInsulation ivi;
	ivi.depth = 0.1;  // [m]
	ivi.layer = xps10;

	// Exterior Horizontal Insulation
	HorizontalInsulation ehi;
	ehi.depth = 0.1;  // [m]
	ehi.width = 0.6;  // [m]
	ehi.layer = null;

	// Exterior Vertical Insulation
	VerticalInsulation evi;
	evi.depth = 0.61;  // [m]
	evi.layer = xps10;


	// Foundations
	Foundation foundation1;

	foundation1.soil = drySandySoil;
	foundation1.slab = slab;
	foundation1.hasSlab = true;
	foundation1.wall = wall;
	foundation1.hasWall = true;
	foundation1.interiorHorizontalInsulation = ihi;
	foundation1.hasInteriorHorizontalInsulation = false;
	foundation1.interiorVerticalInsulation = ivi;
	foundation1.hasInteriorVerticalInsulation = false;
	foundation1.exteriorHorizontalInsulation = ehi;
	foundation1.hasExteriorHorizontalInsulation = false;
	foundation1.exteriorVerticalInsulation = evi;
	foundation1.hasExteriorVerticalInsulation = false;

	// Site
	//foundation1.farFieldWidth = 5.0;  // [m]
	foundation1.farFieldWidth = 20.0 + foundation1.wall.layers[0].d/2;  // [m]
	foundation1.deepGroundDepth = 30.0;  // [m]
	foundation1.deepGroundTemperature = 14.0 + 273.15;  // [K]
	foundation1.deepGroundBoundary = Foundation::DGB_AUTO;
	foundation1.excavationDepth = 0.0;  // [m]

	double tempInF = 70.0;  // [F]
	foundation1.indoorAirTemperature = (tempInF - 32.0)*5.0/9.0 + 273.15;  // [K]

	// Geometry
	double length = 10 - foundation1.wall.layers[0].d; // [m]
	double width = length; // [m]
	foundation1.area = length*width;  // [m2]
	foundation1.perimeter = 2*length + 2*width;  // [m]


	// Meshing
	Mesh uniformMesh;
	uniformMesh.minCellDim = 0.05; // [m]
	uniformMesh.maxDepthGrowthCoeff = 1.05;
	uniformMesh.maxInteriorGrowthCoeff = 1.05;
	uniformMesh.maxRadialGrowthCoeff = 1.05;

	foundation1.mesh = uniformMesh;

	// Simulation Control
	foundation1.numericalScheme = Foundation::NS_ADE;
	foundation1.initialTemperature = 10.0 + 273.15;  // [K]
	foundation1.initializationMethod = Foundation::IM_KUSUDA;
	foundation1.interiorConvectiveCoefficient = 1000000;  // [W/m2-K]
	foundation1.exteriorConvectiveCoefficient = 1000000;  // [W/m2-K]
	foundation1.convectionCalculationMethod = Foundation::CCM_AUTO;
	foundation1.outdoorDryBulbTemperature = 8 + 273.15;  // [K]
	foundation1.outdoorTemperatureMethod = Foundation::OTM_WEATHER_FILE;

	// Output
	OutputAnimation output1;
	output1.name = "Test";
	output1.frequency = hours(6) + minutes(0);
	output1.grid = true;
	foundation1.outputAnimation = output1;

	// Full Input
	input.simulationControl = simulationControl;
	input.foundations.push_back(foundation1);


	return input;

}


