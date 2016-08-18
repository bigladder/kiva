/* Copyright (c) 2012-2016 Big Ladder Software. All rights reserved.
* See the LICENSE file for additional terms and conditions. */

#ifndef GROUNDPLOT_H_
#define GROUNDPLOT_H_

#include <mgl2/mgl.h>

#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <boost/filesystem/operations.hpp>

#include "Input.hpp"
#include "Domain.hpp"
#include "Functions.hpp"
#include "Geometry.hpp"

using namespace Kiva;

struct Axis
{
  std::size_t nMin;
  std::size_t nMax;
  std::size_t nN;
  Mesher mesh;
};

enum SliceType
{
  XZ_2D,
  XY,
  XZ,
  YZ
};

class GroundPlot
{
private:
  mglData hDat, vDat, cDat, hGrid, vGrid, TGrid;


  double slice;
  double xMin, xMax, yMin, yMax;

  Axis hAxis;
  Axis vAxis;

  int frameNumber;


public:

  //mglGraph gr;
  OutputAnimation outputAnimation;
  std::vector<Block> blocks;
  SliceType sliceType;

  mglData TDat;
  std::size_t iMin, iMax, jMin, jMax, kMin, kMax;

  double distanceUnitConversion;

  boost::posix_time::ptime tStart, tEnd;
  boost::posix_time::ptime nextPlotTime;
  GroundPlot(OutputAnimation &outputAnimation, Domain &domain, std::vector<Block> &blocks);
  void createFrame(std::string timeStamp);
  bool makeNewFrame(boost::posix_time::ptime t);
};


#endif /* GROUNDPLOT_H_ */
