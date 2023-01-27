/* Copyright (c) 2012-2022 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#include <cmath>
#include <iostream>

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations" // sprintf in lexical_cast
#endif
#include <boost/date_time/posix_time/posix_time.hpp>
#if defined(__clang__)
#pragma clang diagnostic pop
#endif
#include <boost/program_options.hpp>

#include "Errors.hpp"
#include "InputParser.hpp"
#include "Simulator.hpp"
#include "Version.hpp"
#include "WeatherData.hpp"

namespace po = boost::program_options;

void errorCallback(const int messageType, const std::string message, void *) {
  if (messageType == Kiva::MSG_INFO) {
    std::cout << message << std::endl;
  } else if (messageType == Kiva::MSG_WARN) {
    std::cout << "WARNING: " << message << std::endl;
  } else /* if (messageType == Kiva::MSG_ERR) */ {
    std::cout << "ERROR: " << message << std::endl;
    exit(EXIT_FAILURE);
  }
}

int main(int argc, char *argv[]) {
  Kiva::setMessageCallback(errorCallback, NULL);

  std::string versionInfo = "kiva ";
  versionInfo.append(Kiva::getVersion());
  std::string copyrightInfo = "Copyright (C) 2012-";
  copyrightInfo.append(Kiva::getYear());
  copyrightInfo.append(" Big Ladder Software LLC\n"
                       "Web: bigladdersoftware.com");
  std::string usageInfo = "Usage: kiva [Input File] [Weather File] [Output File]\n"
                          "   Input format: yaml\n"
                          "   Weather format: epw\n"
                          "   Output format: csv";

  try {

    po::options_description generic("Options");
    generic.add_options()("help,h", "Produce this message")("version,v",
                                                            "Display version information");

    po::options_description hidden("Hidden options");
    hidden.add_options()("input-file", po::value<std::string>(),
                         "input file")("weather-file", po::value<std::string>(), "weather file")(
        "output-file", po::value<std::string>(), "output file");

    po::options_description cmdLine;
    cmdLine.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("input-file", 1).add("weather-file", 1).add("output-file", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(cmdLine).positional(p).run(), vm);
    po::notify(vm);

    if (vm.empty()) {
      std::stringstream ss;
      ss << versionInfo << "\n" << copyrightInfo << "\n\n" << usageInfo << "\n" << generic;

      Kiva::showMessage(MSG_INFO, ss.str());
    } else if (vm.count("help")) {
      std::stringstream ss;
      ss << versionInfo << "\n" << copyrightInfo << "\n\n" << usageInfo << "\n" << generic;

      Kiva::showMessage(MSG_INFO, ss.str());
    } else if (vm.count("version")) {
      Kiva::showMessage(MSG_INFO, versionInfo + "\n" + copyrightInfo);
    } else if (vm.count("input-file") && vm.count("weather-file") && vm.count("output-file")) {
      Kiva::showMessage(MSG_INFO, versionInfo + "\n" + copyrightInfo);

      boost::posix_time::ptime beginCalc = boost::posix_time::microsec_clock::local_time();
      Kiva::showMessage(MSG_INFO, "Starting Program: " + to_simple_string(beginCalc));

      // parse input
      Input input = inputParser(vm["input-file"].as<std::string>());

      // parse weather
      WeatherData weather(vm["weather-file"].as<std::string>());

      input.simulationControl.setStartTime();

      // initialize
      Simulator simulator(weather, input, vm["output-file"].as<std::string>());

      simulator.simulate();

      boost::posix_time::ptime finishCalc = boost::posix_time::microsec_clock::local_time();
      Kiva::showMessage(MSG_INFO, "Finished Program: " + to_simple_string(finishCalc));

      boost::posix_time::time_duration totalCalc = finishCalc - beginCalc;
      Kiva::showMessage(MSG_INFO, "Elapsed Time: " + to_simple_string(totalCalc));

    } else if (!vm.empty()) {
      std::stringstream ss;
      ss << "Incorrect number of arguments\n\n" << usageInfo + "\n" << generic;

      Kiva::showMessage(MSG_ERR, ss.str());
    }

    return 0;

  } catch (std::exception &e) {
    Kiva::showMessage(MSG_ERR, e.what());
  }
}
