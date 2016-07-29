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
  headers.resize(18);

  headers[0] = "Slab Core Average Heat Flux [W/m2]";
  headers[1] = "Slab Core Average Temperature [K]";
  headers[2] = "Slab Core Average Effective Temperature [C]";
  headers[3] = "Slab Core Total Heat Transfer Rate [W]";
  headers[4] = "Slab Perimeter Average Heat Flux [W/m2]";
  headers[5] = "Slab Perimeter Average Temperature [K]";
  headers[6] = "Slab Perimeter Average Effective Temperature [C]";
  headers[7] = "Slab Perimeter Total Heat Transfer Rate [W]";
  headers[8] = "Slab Average Heat Flux [W/m2]";
  headers[9] = "Slab Average Temperature [K]";
  headers[10] = "Slab Total Heat Transfer Rate [W]";
  headers[11] = "Wall Average Heat Flux [W/m2]";
  headers[12] = "Wall Average Temperature [K]";
  headers[13] = "Wall Average Effective Temperature [C]";
  headers[14] = "Wall Total Heat Transfer Rate [W]";
  headers[15] = "Foundation Average Heat Flux [W/m2]";
  headers[16] = "Foundation Average Temperature [K]";
  headers[17] = "Foundation Total Heat Transfer Rate [W]";

  variableID = varID;
  headerText = headers[varID];
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
