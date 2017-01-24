/* Copyright (c) 2012-2016 Big Ladder Software. All rights reserved.
* See the LICENSE file for additional terms and conditions. */

#include <iostream>
#include <fstream>
#include <cmath>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/program_options.hpp>

#include "Version.hpp"
#include "InputParser.hpp"
#include "WeatherData.hpp"
#include "Simulator.hpp"
#include "libkiva_export.h"

namespace po = boost::program_options;

int main(int argc, char *argv[])
{
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

  try
  {

    po::options_description generic("Options");
    generic.add_options()
        ("help,h", "Produce this message")
        ("version,v", "Display version information");

        po::options_description hidden("Hidden options");
        hidden.add_options()
            ("input-file", po::value<std::string>(), "input file")
            ("weather-file", po::value<std::string>(), "weather file")
            ("output-file", po::value<std::string>(), "output file");

        po::options_description cmdLine;
        cmdLine.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("input-file", 1).add("weather-file", 1).add("output-file", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
        options(cmdLine).positional(p).run(), vm);
    po::notify(vm);

    if (vm.empty())
    {
      std::cout << versionInfo << "\n";
      std::cout << copyrightInfo << "\n\n";
      std::cout << usageInfo << "\n";
      std::cout << generic;
    }
    else if (vm.count("help"))
    {
      std::cout << versionInfo << "\n";
      std::cout << copyrightInfo << "\n\n";
      std::cout << usageInfo << "\n";
      std::cout << generic;
    }
    else if (vm.count("version"))
    {
      std::cout << versionInfo << "\n";
      std::cout << copyrightInfo << std::endl;
    }
    else if (vm.count("input-file") && vm.count("weather-file") && vm.count("output-file"))
    {
      std::cout << versionInfo << "\n";
      std::cout << copyrightInfo << "\n" << std::endl;

      boost::posix_time::ptime beginCalc = boost::posix_time::microsec_clock::local_time();
      std::cout << "Starting Program: " << beginCalc << std::endl;

      // parse input
      Input input = inputParser(vm["input-file"].as<std::string>());

      // parse weather
      WeatherData weather(vm["weather-file"].as<std::string>());

      input.simulationControl.setStartTime();

      // initialize
      Simulator simulator(weather,input,vm["output-file"].as<std::string>());

      simulator.simulate();

      boost::posix_time::ptime finishCalc = boost::posix_time::microsec_clock::local_time();
      std::cout << "Finished Program: " << finishCalc << std::endl;

      boost::posix_time::time_duration totalCalc = finishCalc - beginCalc;
      std::cout << "Elapsed Time: " << totalCalc << std::endl;

    }
    else if (!vm.empty())
    {
      std::cerr << "ERROR: Incorrect number of arguments\n\n";

      std::cout << usageInfo << "\n";
      std::cout << generic;

      return 1;
    }

    return 0;

  }
  catch(std::exception& e)
  {
    std::cerr << e.what() << std::endl;
    return 1;
  }

}
