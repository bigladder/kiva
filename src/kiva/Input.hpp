/* Copyright (c) 2012-2019 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#ifndef INPUT_HPP_
#define INPUT_HPP_

#include <fstream>

#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include <boost/lexical_cast.hpp>

#include "Algorithms.hpp"
#include "Foundation.hpp"
#include "Geometry.hpp"
#include "GroundOutput.hpp"
#include "GroundPlot.hpp"
#include "Mesher.hpp"
#include "WeatherData.hpp"

using namespace Kiva;

class SimulationControl {
public:
  // Simulation Control
  boost::gregorian::date startDate;
  boost::gregorian::date endDate;
  boost::posix_time::time_duration timestep;
  std::string weatherFile;
  boost::posix_time::ptime startTime;

  void setStartTime();
};

class Initialization {
public:
  double initialTemperature;
  enum InitializationMethod { IM_KUSUDA, IM_CONSTANT_TEMPERATURE, IM_STEADY_STATE };

  long warmupDays;
  long implicitAccelTimestep;
  long implicitAccelPeriods;

  InitializationMethod initializationMethod;
};

class OutputSnapshots {
public:
  boost::gregorian::date startDate;
  boost::gregorian::date endDate;
  SnapshotSettings snapshotSettings;
  bool startDateSet;
  bool endDateSet;
  bool xRangeSet;
  bool yRangeSet;
  bool zRangeSet;
};

class OutputVariable {
private:
  static const std::vector<std::string> headers;
  static const std::vector<std::vector<Surface::SurfaceType>> surfaceTypes;
  static const std::vector<GroundOutput::OutputType> outTypes;

public:
  int variableID;
  std::string headerText;
  std::vector<Surface::SurfaceType> surfaces;
  GroundOutput::OutputType outType;

  OutputVariable(int varID);
};

class OutputReport : public std::vector<OutputVariable> {
public:
  void setOutputMap();
  std::vector<Surface::SurfaceType> outputMap;
  boost::posix_time::time_duration minFrequency;
};

class Output {
public:
  OutputReport outputReport;
  std::vector<OutputSnapshots> outputSnapshots;
};

class DataFile {
public:
  std::string fileName;
  std::pair<int, int> firstIndex;
  HourlyData data;
  std::filesystem::path searchDir;

  void readData();
};

class Boundaries {
public:
  double indoorAirTemperature; // [K]
  DataFile indoorAirTemperatureFile;

  enum IndoorTemperatureMethod { ITM_FILE, ITM_CONSTANT_TEMPERATURE };

  IndoorTemperatureMethod indoorTemperatureMethod;

  // Local wind speed characteristics
  double deltaLocal; // [m]
  double alphaLocal; // [-]

  double outdoorDryBulbTemperature;
  enum OutdoorTemperatureMethod { OTM_WEATHER_FILE, OTM_CONSTANT_TEMPERATURE };

  OutdoorTemperatureMethod outdoorTemperatureMethod;

  // Deep Ground options
  double deepGroundTemperature;
  enum DeepGroundBoundaryType { DGBT_AUTO, DGBT_CONSTANT_TEMPERATURE, DGBT_ZERO_FLUX };

  DeepGroundBoundaryType deepGroundBoundaryType;

  // Convection
  ConvectionAlgorithm interiorConvectionAlgorithm;
  ConvectionAlgorithm exteriorConvectionAlgorithm;
};

class Input {
public:
  SimulationControl simulationControl;
  Foundation foundation;
  Boundaries boundaries;
  Initialization initialization;
  Output output;
};

#endif /* INPUT_HPP_ */
