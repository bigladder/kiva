/* Copyright (c) 2012-2016 Big Ladder Software. All rights reserved.
* See the LICENSE file for additional terms and conditions. */

#ifndef INPUTPARSER_HPP_
#define INPUTPARSER_HPP_

#include "Input.hpp"
#include "yaml-cpp/yaml.h"
#include "WeatherData.hpp"
#include <fstream>

Input inputParser(std::string inputFile);


#endif /* INPUTPARSER_HPP_ */
