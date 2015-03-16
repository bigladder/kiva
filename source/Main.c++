/* Main.c++ is part of Kiva (Written by Neal Kruis)
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

#include <iostream>
#include <fstream>
#include <cmath>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/program_options.hpp>

#if defined(USE_LIS_SOLVER)
#include "lis.h"
#endif

#include "Version.h"
#include "InputParser.h"
#include "WeatherData.h"
#include "Input.h"
#include "Ground.h"

namespace po = boost::program_options;

int main(int argc, char *argv[])
{
#if defined(USE_LIS_SOLVER)
  lis_initialize(&argc, &argv);
#endif

  std::string versionInfo = "kiva ";
  versionInfo.append(VERSION_NUMBER);
  std::string copyrightInfo = "Copyright (C) 2012-2014 Big Ladder Software\n"
                           "Web: www.bigladdersoftware.com";
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
      Ground ground(weather,input.foundations[0],input.simulationControl,
          vm["output-file"].as<std::string>());

      ground.simulate();

      boost::posix_time::ptime finishCalc = boost::posix_time::microsec_clock::local_time();
      std::cout << "Finished Program: " << finishCalc << std::endl;

      boost::posix_time::time_duration totalCalc = finishCalc - beginCalc;
      std::cout << "Elapsed Time: " << totalCalc << std::endl;

    }
    else if (!vm.empty())
    {
      std::cout << "ERROR: Incorrect number of arguments\n\n";

      std::cout << usageInfo << "\n";
      std::cout << generic;

#if defined(USE_LIS_SOLVER)
      lis_finalize();
#endif
      return 1;
    }

#if defined(USE_LIS_SOLVER)
    lis_finalize();
#endif
    return 0;

  }
  catch(std::exception& e)
  {
    std::cout << e.what() << std::endl;
#if defined(USE_LIS_SOLVER)
    lis_finalize();
#endif
    return 1;
  }

}



