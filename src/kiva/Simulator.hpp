/* Copyright (c) 2012-2022 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#ifndef Simulator_HPP
#define Simulator_HPP

#include <iostream>

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations" // sprintf in lexical_cast
#endif
#include <boost/date_time/posix_time/posix_time.hpp>
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

#include "BoundaryConditions.hpp"
#include "Geometry.hpp"
#include "Ground.hpp"
#include "GroundOutput.hpp"
#include "GroundPlot.hpp"
#include "Input.hpp"
#include "WeatherData.hpp"

using namespace Kiva;

class Simulator {
public:
  // Constructor
  Simulator(WeatherData &weatherData, Input &input, std::string outputFileName);

  virtual ~Simulator();
  void simulate();

  WeatherData &weatherData;
  Input &input;

  double annualAverageDryBulbTemperature;

  double percentComplete;

private:
  Ground ground;
  BoundaryConditions bcs;

  std::vector<GroundPlot> plots;
  std::ofstream outputFile;
  filesys::path outputDir;
  void initializePlots();
  void initializeConditions();

  void printStatus(boost::posix_time::ptime t);

  std::string printOutputHeaders();
  std::string printOutputLine();

  void plot(boost::posix_time::ptime t);

  boost::posix_time::ptime prevStatusUpdate;
  boost::posix_time::ptime prevOutputTime;
  bool initPeriod;

  double getInitialTemperature(boost::posix_time::ptime t, double z);

  void updateBoundaryConditions(boost::posix_time::ptime t);
};

#endif // Simulator_HPP
