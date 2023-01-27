/* Copyright (c) 2012-2022 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#include "Simulator.hpp"
#include "Errors.hpp"

#include <fmt/core.h>

using namespace Kiva;

static const double PI = 4.0 * atan(1.0);

Simulator::Simulator(WeatherData &weatherData, Input &input, std::string outputFileName)
    : weatherData(weatherData), input(input),
      ground(input.foundation, input.output.outputReport.outputMap) {
  // set up output file
  std::filesystem::path outputPath(outputFileName);
  outputDir = outputPath.parent_path();
  outputFile.open(outputFileName.c_str());
  outputFile << "Timestamp" << printOutputHeaders() << std::endl;

  annualAverageDryBulbTemperature = weatherData.dryBulbTemp.getAverage();

  showMessage(MSG_INFO, "Creating Domain...");

  if (input.boundaries.deepGroundBoundaryType == Boundaries::DGBT_AUTO)
    input.boundaries.deepGroundTemperature = annualAverageDryBulbTemperature;

  if (!input.foundation.useDetailedExposedPerimeter || !isConvex(input.foundation.polygon)) {
    if (input.foundation.reductionStrategy == Foundation::RS_BOUNDARY) {
      input.foundation.reductionStrategy = Foundation::RS_AP;
    }
  }

  bcs = input.boundaries.bcs;

  if (input.foundation.reductionStrategy == Foundation::RS_BOUNDARY) {
    ground.calculateBoundaryLayer();
    ground.setNewBoundaryGeometry();
  }

  ground.buildDomain();

  std::stringstream ss;

  ss << "  X Cells: " << ground.nX << "\n"
     << "  Y Cells: " << ground.nY << "\n"
     << "  Z Cells: " << ground.nZ << "\n"
     << "  Total Cells: " << ground.nX * ground.nY * ground.nZ;

  showMessage(MSG_INFO, ss.str());

  // Initial Conditions
  initializeConditions();

  initializePlots();
}

Simulator::~Simulator() { outputFile.close(); }

void Simulator::initializeConditions() {

  prevStatusUpdate = boost::posix_time::second_clock::local_time();

  initPeriod = true;

  if (input.foundation.numericalScheme !=
      Foundation::NS_STEADY_STATE) // Intialization not necessary for steady state calculations
  {
    showMessage(MSG_INFO, "Initializing Temperatures...");

    // Calculate initial time in seconds (simulation start minus warmup and acceleration periods)
    boost::posix_time::time_duration &simulationTimestep = input.simulationControl.timestep;
    boost::posix_time::time_duration accelTimestep =
        boost::posix_time::hours(input.initialization.implicitAccelTimestep);

    boost::posix_time::time_duration accelDuration =
        accelTimestep * input.initialization.implicitAccelPeriods;
    boost::posix_time::time_duration warmupDuration =
        boost::posix_time::hours(input.initialization.warmupDays * 24);

    boost::posix_time::ptime tInit = input.simulationControl.startTime - warmupDuration -
                                     simulationTimestep - accelDuration - accelTimestep;

    // Calculate initial conditions
    if (input.initialization.initializationMethod == Initialization::IM_STEADY_STATE) {
      Foundation::NumericalScheme tempNS = input.foundation.numericalScheme;
      input.foundation.numericalScheme = Foundation::NS_STEADY_STATE;
      updateBoundaryConditions(tInit);
      ground.calculate(bcs);
      printStatus(tInit);
      input.foundation.numericalScheme = tempNS;
    } else {
      for (size_t i = 0; i < ground.nX; ++i) {
        for (size_t j = 0; j < ground.nY; ++j) {
          for (size_t k = 0; k < ground.nZ; ++k) {
            std::size_t index = ground.domain.getIndex(i, j, k);
            ground.TOld[index] = getInitialTemperature(tInit, ground.domain.mesh[2].centers[k]);
          }
        }
      }
    }

    // Calculate implicit acceleration
    if (input.initialization.implicitAccelPeriods > 0) {
      boost::posix_time::ptime tAccelStart = input.simulationControl.startTime - warmupDuration -
                                             simulationTimestep -
                                             accelDuration; // [s] Acceleration start time
      boost::posix_time::ptime tAccelEnd = input.simulationControl.startTime - warmupDuration -
                                           simulationTimestep; // [s] Acceleration end time

      Foundation::NumericalScheme tempNS = input.foundation.numericalScheme;
      input.foundation.numericalScheme = Foundation::NS_IMPLICIT;

      for (boost::posix_time::ptime t = tAccelStart; t <= tAccelEnd; t += accelTimestep) {
        updateBoundaryConditions(t);
        ground.calculate(bcs, static_cast<double>(accelTimestep.total_seconds()));
        printStatus(t);
      }

      input.foundation.numericalScheme = tempNS;
    }

    // Calculate warmup
    if (input.initialization.warmupDays > 0) {

      boost::posix_time::ptime tWarmupStart =
          input.simulationControl.startTime - warmupDuration; // [s] Acceleration start time
      boost::posix_time::ptime tWarmupEnd =
          input.simulationControl.startTime - simulationTimestep; // [s] Simulation end time

      for (boost::posix_time::ptime t = tWarmupStart; t <= tWarmupEnd; t += simulationTimestep) {
        updateBoundaryConditions(t);
        ground.calculate(bcs, static_cast<double>(simulationTimestep.total_seconds()));
        printStatus(t);
      }
    }
  }
  initPeriod = false;
}

void Simulator::initializePlots() {
  for (std::size_t p = 0; p < input.output.outputSnapshots.size(); p++) {
    if (!input.output.outputSnapshots[p].startDateSet)
      input.output.outputSnapshots[p].startDate = input.simulationControl.startDate;

    if (!input.output.outputSnapshots[p].endDateSet)
      input.output.outputSnapshots[p].endDate = input.simulationControl.endDate;

    if (ground.foundation.numberOfDimensions == 3) {
      if (!input.output.outputSnapshots[p].xRangeSet) {
        input.output.outputSnapshots[p].snapshotSettings.xRange.first =
            ground.domain.mesh[0].dividers[0];
        input.output.outputSnapshots[p].snapshotSettings.xRange.second =
            ground.domain.mesh[0].dividers[ground.nX];
      }

      if (!input.output.outputSnapshots[p].yRangeSet) {
        input.output.outputSnapshots[p].snapshotSettings.yRange.first =
            ground.domain.mesh[1].dividers[0];
        input.output.outputSnapshots[p].snapshotSettings.yRange.second =
            ground.domain.mesh[1].dividers[ground.nY];
      }

      if (!input.output.outputSnapshots[p].zRangeSet) {

        if (!input.output.outputSnapshots[p].xRangeSet &&
            !input.output.outputSnapshots[p].yRangeSet) {
          input.output.outputSnapshots[p].snapshotSettings.zRange.first = 0;
          input.output.outputSnapshots[p].snapshotSettings.zRange.second = 0;
        } else {
          input.output.outputSnapshots[p].snapshotSettings.zRange.first =
              ground.domain.mesh[2].dividers[0];
          input.output.outputSnapshots[p].snapshotSettings.zRange.second =
              ground.domain.mesh[2].dividers[ground.nZ];
        }
      }

    } else if (ground.foundation.numberOfDimensions == 2) {
      if (!input.output.outputSnapshots[p].xRangeSet) {
        input.output.outputSnapshots[p].snapshotSettings.xRange.first =
            ground.domain.mesh[0].dividers[0];
        input.output.outputSnapshots[p].snapshotSettings.xRange.second =
            ground.domain.mesh[0].dividers[ground.nX];
      }

      if (!input.output.outputSnapshots[p].yRangeSet) {
        input.output.outputSnapshots[p].snapshotSettings.yRange.first = 0.5;
        input.output.outputSnapshots[p].snapshotSettings.yRange.second = 0.5;
      }

      if (!input.output.outputSnapshots[p].zRangeSet) {
        input.output.outputSnapshots[p].snapshotSettings.zRange.first =
            ground.domain.mesh[2].dividers[0];
        input.output.outputSnapshots[p].snapshotSettings.zRange.second =
            ground.domain.mesh[2].dividers[ground.nZ];
      }
    } else {
      if (!input.output.outputSnapshots[p].xRangeSet) {
        input.output.outputSnapshots[p].snapshotSettings.xRange.first = 0.5;
        input.output.outputSnapshots[p].snapshotSettings.xRange.second = 0.5;
      }

      if (!input.output.outputSnapshots[p].yRangeSet) {
        input.output.outputSnapshots[p].snapshotSettings.yRange.first = 0.5;
        input.output.outputSnapshots[p].snapshotSettings.yRange.second = 0.5;
      }

      if (!input.output.outputSnapshots[p].zRangeSet) {
        input.output.outputSnapshots[p].snapshotSettings.zRange.first =
            ground.domain.mesh[2].dividers[0];
        input.output.outputSnapshots[p].snapshotSettings.zRange.second =
            ground.domain.mesh[2].dividers[ground.nZ];
      }
    }

    input.output.outputSnapshots[p].snapshotSettings.dir =
        (outputDir / input.output.outputSnapshots[p].snapshotSettings.dir).string();

    plots.emplace_back(input.output.outputSnapshots[p].snapshotSettings, ground.domain,
                       input.foundation);

    boost::posix_time::ptime startTime(input.output.outputSnapshots[p].startDate,
                                       boost::posix_time::hours(0));
    boost::posix_time::ptime endTime(input.output.outputSnapshots[p].endDate +
                                     boost::gregorian::days(1));

    plots[p].tStart =
        static_cast<double>((startTime - input.simulationControl.startTime).total_seconds());
    plots[p].nextPlotTime =
        static_cast<double>((startTime - input.simulationControl.startTime).total_seconds());
    plots[p].tEnd =
        static_cast<double>((endTime - input.simulationControl.startTime).total_seconds());
  }
}

void Simulator::simulate() {
  showMessage(MSG_INFO, "Beginning Simulation...");

  boost::posix_time::ptime simStart = input.simulationControl.startTime;
  boost::posix_time::ptime simEnd(input.simulationControl.endDate + boost::gregorian::days(1));
  boost::posix_time::time_duration simDuration = simEnd - simStart;

  prevOutputTime = input.simulationControl.startTime - input.output.outputReport.minFrequency;
  double timestep = static_cast<double>(input.simulationControl.timestep.total_seconds());

  for (boost::posix_time::ptime t = simStart; t < simEnd;
       t = t + input.simulationControl.timestep) {

    percentComplete =
        round(double((t - simStart).total_seconds()) / double(simDuration.total_seconds()) * 1000) /
        10.0;
    updateBoundaryConditions(t);
    ground.calculate(bcs, timestep);
    ground.calculateSurfaceAverages();
    plot(t);
    printStatus(t);

    if (t - prevOutputTime >= input.output.outputReport.minFrequency) {
      outputFile << to_simple_string(t) << printOutputLine() << std::endl;
      prevOutputTime = t;
    }
  }

  showMessage(MSG_INFO,
              "  " + to_simple_string(simEnd - input.simulationControl.timestep) + " (100%)");
}

void Simulator::plot(boost::posix_time::ptime t) {
  for (std::size_t p = 0; p < plots.size(); p++) {
    if (plots[p].makeNewFrame(
            static_cast<double>((t - input.simulationControl.startTime).total_seconds()))) {
      std::string timeStamp = to_simple_string(t);

      plots[p].createFrame(ground, timeStamp.substr(5, timeStamp.size() - 5));
    }
  }
}

void Simulator::printStatus(boost::posix_time::ptime t) {
  boost::posix_time::ptime currentTime = boost::posix_time::second_clock::local_time();

  if (currentTime - prevStatusUpdate > boost::posix_time::milliseconds(500)) {
    if (initPeriod) {
      showMessage(MSG_INFO, "  " + to_simple_string(t));
    } else {
      std::stringstream ss;
      ss << "  " << t << " (" << percentComplete << "%)";
      showMessage(MSG_INFO, ss.str());
    }

    prevStatusUpdate = currentTime;
  }
}

double Simulator::getInitialTemperature(boost::posix_time::ptime t, double z) {
  if (input.initialization.initializationMethod == Initialization::IM_KUSUDA) {
    boost::gregorian::greg_year year = t.date().year();
    boost::gregorian::date dayBegin(year, boost::gregorian::Jan, 1);
    boost::posix_time::ptime tYearStart(dayBegin);
    boost::posix_time::time_duration tSinceYearStart = t - tYearStart;
    double trel = static_cast<double>(tSinceYearStart.total_seconds());
    double tshift = weatherData.hourOfMinimumTemperature * 60.0 * 60.0;
    double seconds_in_day = 60.0 * 60.0 * 24.0;
    double Tamp = (weatherData.maximumAverageMontlyTemperature -
                   weatherData.minimumAverageMontlyTemperature) /
                  2.0;
    double diff = input.foundation.soil.conductivity /
                  (input.foundation.soil.density * input.foundation.soil.specificHeat);

    return annualAverageDryBulbTemperature -
           Tamp * exp(z * sqrt(PI / (365 * seconds_in_day * diff))) *
               cos((2 * PI) / (365 * seconds_in_day) *
                   (trel - tshift + z * 0.5 * sqrt((365 * seconds_in_day) / (PI * diff))));
  } else // Initialization::IM_CONSTANT_TEMPERATURE)
    return input.initialization.initialTemperature;
}

void Simulator::updateBoundaryConditions(boost::posix_time::ptime t) {
  double Tin;
  if (input.boundaries.indoorTemperatureMethod == Boundaries::ITM_FILE)
    Tin = input.boundaries.indoorAirTemperatureFile.data.getValue(t);
  else // Boundaries::ITM_CONSTANT_TEMPERATURE)
    Tin = input.boundaries.indoorAirTemperature;

  bcs.slabConvectiveTemp = bcs.wallConvectiveTemp = bcs.slabRadiantTemp = bcs.wallRadiantTemp = Tin;

  if (input.boundaries.outdoorTemperatureMethod == Boundaries::OTM_WEATHER_FILE)
    bcs.outdoorTemp = weatherData.dryBulbTemp.getValue(t);
  else // Boundaries::OTM_CONSTANT_TEMPERATURE)
    bcs.outdoorTemp = input.boundaries.outdoorDryBulbTemperature;

  double vWS = weatherData.windSpeed.getValue(t);
  const double deltaWS = 270; // [m]
  const double alphaWS = 0.14;
  const double deltaLocal = input.boundaries.deltaLocal; // [m]
  const double alphaLocal = input.boundaries.alphaLocal;
  const double zWS = 10;                                  // [m]
  const double zLocal = input.foundation.grade.roughness; // [m]
  const double vMult = pow(deltaWS / zWS, alphaWS) * pow(zLocal / deltaLocal, alphaLocal);

  bcs.localWindSpeed = vWS * vMult;
  bcs.windDirection = weatherData.windDirection.getValue(t);

  bcs.solarAzimuth = weatherData.azimuth.getValue(t);
  bcs.solarAltitude = weatherData.altitude.getValue(t);
  bcs.directNormalFlux = weatherData.directNormalSolar.getValue(t);
  bcs.diffuseHorizontalFlux = weatherData.diffuseHorizontalSolar.getValue(t);
  bcs.skyEmissivity = weatherData.skyEmissivity.getValue(t);

  bcs.deepGroundTemperature = input.boundaries.deepGroundTemperature;
}

std::string Simulator::printOutputHeaders() {
  std::string outputHeader = "";

  for (size_t o = 0; o < input.output.outputReport.size(); o++) {
    outputHeader += ", " + input.output.outputReport[o].headerText;
  }

  return outputHeader;
}

std::string Simulator::printOutputLine() {
  std::string outputLine = "";

  for (const auto &out : input.output.outputReport) {

    double totalValue = 0.0;
    double totalVA = 0.0;
    double totalArea = 0.0;
    for (auto surface : out.surfaces) {
      if (ground.foundation.hasSurface[surface]) {
        totalValue += ground.getSurfaceAverageValue({surface, out.outType});
        totalVA += ground.getSurfaceAverageValue({surface, out.outType}) *
                   ground.foundation.surfaceAreas[surface];
        totalArea += ground.foundation.surfaceAreas[surface];
      }
    }

    if (totalArea > 0.0) {
      double value = out.outType == GroundOutput::OT_RATE ? totalValue : totalVA / totalArea;
      outputLine += fmt::format(", {}", value);
    } else {
      outputLine += ", NAN";
    }
  }

  return outputLine;
}
