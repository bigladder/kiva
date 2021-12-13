/* Copyright (c) 2012-2021 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#ifndef INPUT_CPP_
#define INPUT_CPP_

#include "Input.hpp"
#include "Errors.hpp"

using namespace Kiva;

void SimulationControl::setStartTime() {
  boost::posix_time::ptime st(startDate, boost::posix_time::hours(0));
  startTime = st;
}

const std::vector<std::string> OutputVariable::headers = {
    "Slab Core Average Heat Flux [W/m2]",               // 0
    "Slab Core Average Temperature [K]",                // 1
    "Slab Core Average Effective Temperature [C]",      // 2
    "Slab Core Total Heat Transfer Rate [W]",           // 3
    "Slab Perimeter Average Heat Flux [W/m2]",          // 4
    "Slab Perimeter Average Temperature [K]",           // 5
    "Slab Perimeter Average Effective Temperature [C]", // 6
    "Slab Perimeter Total Heat Transfer Rate [W]",      // 7
    "Slab Average Heat Flux [W/m2]",                    // 8
    "Slab Average Temperature [K]",                     // 9
    "Slab Total Heat Transfer Rate [W]",                // 10
    "Wall Average Heat Flux [W/m2]",                    // 11
    "Wall Average Temperature [K]",                     // 12
    "Wall Average Effective Temperature [C]",           // 13
    "Wall Total Heat Transfer Rate [W]",                // 14
    "Foundation Average Heat Flux [W/m2]",              // 15
    "Foundation Average Temperature [K]",               // 16
    "Foundation Total Heat Transfer Rate [W]"           // 17
};

const std::vector<std::vector<Surface::SurfaceType>> OutputVariable::surfaceTypes = {
    {Surface::ST_SLAB_CORE},
    {Surface::ST_SLAB_CORE},
    {Surface::ST_SLAB_CORE},
    {Surface::ST_SLAB_CORE},
    {Surface::ST_SLAB_PERIM},
    {Surface::ST_SLAB_PERIM},
    {Surface::ST_SLAB_PERIM},
    {Surface::ST_SLAB_PERIM},
    {Surface::ST_SLAB_CORE, Surface::ST_SLAB_PERIM},
    {Surface::ST_SLAB_CORE, Surface::ST_SLAB_PERIM},
    {Surface::ST_SLAB_CORE, Surface::ST_SLAB_PERIM},
    {Surface::ST_WALL_INT},
    {Surface::ST_WALL_INT},
    {Surface::ST_WALL_INT},
    {Surface::ST_WALL_INT},
    {Surface::ST_SLAB_CORE, Surface::ST_SLAB_PERIM, Surface::ST_WALL_INT},
    {Surface::ST_SLAB_CORE, Surface::ST_SLAB_PERIM, Surface::ST_WALL_INT},
    {Surface::ST_SLAB_CORE, Surface::ST_SLAB_PERIM, Surface::ST_WALL_INT}};

const std::vector<GroundOutput::OutputType> OutputVariable::outTypes = {
    GroundOutput::OT_FLUX,     GroundOutput::OT_TEMP,     GroundOutput::OT_EFF_TEMP,
    GroundOutput::OT_RATE,     GroundOutput::OT_FLUX,     GroundOutput::OT_TEMP,
    GroundOutput::OT_EFF_TEMP, GroundOutput::OT_RATE,     GroundOutput::OT_FLUX,
    GroundOutput::OT_TEMP,     GroundOutput::OT_RATE,     GroundOutput::OT_FLUX,
    GroundOutput::OT_TEMP,     GroundOutput::OT_EFF_TEMP, GroundOutput::OT_RATE,
    GroundOutput::OT_FLUX,     GroundOutput::OT_TEMP,     GroundOutput::OT_RATE};

OutputVariable::OutputVariable(int varID) {
  variableID = varID;
  headerText = headers[varID];
  surfaces = surfaceTypes[varID];
  outType = outTypes[varID];
}

void OutputReport::setOutputMap() {
  for (auto outVar : *this) {
    for (auto s : outVar.surfaces) {
      if (std::find(outputMap.begin(), outputMap.end(), s) ==
          outputMap.end()) { // If surface isn't in map create it and add to the list of outputs
        outputMap.push_back(s);
      }
    }
  }
}

void DataFile::readData() {
  // Check relative to working directory first
  std::filesystem::path dataFilePath(fileName);

  if (!std::filesystem::exists(dataFilePath)) {
    // Then check relative to input file Directory

    if (!std::filesystem::exists(searchDir / dataFilePath)) {
      // Print an error and exit
      showMessage(MSG_ERR, "Unable to read data file.");
    } else {
      dataFilePath = searchDir / dataFilePath;
    }
  }

  std::ifstream file(dataFilePath.string().c_str());

  if (!file.is_open()) {
    // Print an error and exit
    showMessage(MSG_ERR, "Unable to read data file.");
  }

  // While there's still stuff left to read
  std::string line;
  std::vector<std::string> columns;

  typedef boost::tokenizer<boost::escaped_list_separator<char>> Tokenizer;
  int row = 0;

  while (!safeGetline(file, line).eof()) {
    row += 1;
    Tokenizer tok(line, boost::escaped_list_separator<char>("\\", ",", "\""));

    columns.assign(tok.begin(), tok.end());

    if (row > firstIndex.first) {
      data.push_back(std::stod(columns[firstIndex.second]));
    }
    columns.clear();
  }
  file.close();
}

#endif /* INPUT_CPP_ */
