/* Copyright (c) 2012-2019 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#ifndef INPUTPARSER_HPP_
#define INPUTPARSER_HPP_

#include <fstream>

#include "yaml-cpp/yaml.h"

#include "Input.hpp"
#include "WeatherData.hpp"

Input inputParser(std::string inputFile);

#endif /* INPUTPARSER_HPP_ */
