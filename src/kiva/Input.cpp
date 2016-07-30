/* Copyright (c) 2012-2016 Big Ladder Software. All rights reserved.
* See the LICENSE file for additional terms and conditions. */

#ifndef INPUT_CPP_
#define INPUT_CPP_

#include "Input.hpp"

static const double PI = 4.0*atan(1.0);

void SimulationControl::setStartTime()
{
  boost::posix_time::ptime st(startDate,boost::posix_time::hours(0));
  startTime = st;
}

OutputVariable::OutputVariable(int varID)
{
  variableID = varID;
  headerText = headers[varID];
  surfaces = surfaceTypes[varID];
  outType = outTypes[varID];
}

void OutputReport::setOutputMap()
{
  for (auto outVar : *this) {
    for (auto s : outVar.surfaces) {
      if (!(outputMap.count(s))) { // If surface isn't in map
        outputMap[s].push_back(outVar.outType);
      }
      else if (std::find(outputMap[s].begin(), outputMap[s].end(), outVar.outType) == outputMap[s].end()) {
        outputMap[s].push_back(outVar.outType);
      }
    }
  }
}


void DataFile::readData()
{
  std::ifstream file(fileName.c_str());

  if (!file.is_open())
  {
      // Print an error and exit
      std::cerr << "Unable to read data file" << std::endl;
      exit(1);
  }

  // While there's still stuff left to read
  std::string line;
  std::vector <std::string> columns;

  typedef boost::tokenizer< boost::escaped_list_separator<char> > Tokenizer;
  int row = 0;

  while (!safeGetline(file,line).eof())
  {
    row += 1;
    Tokenizer tok(line, boost::escaped_list_separator<char>("\\",",","\""));

    columns.assign(tok.begin(), tok.end());

    if (row > firstIndex.first)
    {
      data.push_back(double(boost::lexical_cast<double>(columns[firstIndex.second])));
    }
    columns.clear();
  }
  file.close();
}

#endif /* INPUT_CPP_ */
