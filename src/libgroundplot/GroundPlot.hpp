/* Copyright (c) 2012-2016 Big Ladder Software. All rights reserved.
* See the LICENSE file for additional terms and conditions. */

#ifndef GROUNDPLOT_H_
#define GROUNDPLOT_H_

#include <mgl2/mgl.h>

#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <boost/filesystem/operations.hpp>

#include "../libkiva/Domain.hpp"
#include "../libkiva/Foundation.hpp"
#include "../libkiva/Functions.hpp"
#include "../libkiva/Geometry.hpp"

namespace Kiva {

class SnapshotSettings
{
public:
  SnapshotSettings();
  std::string dir;
  int size;
  double frequency;
  bool grid;
  bool gradients;
  bool contours;
  bool contourLabels;
  bool axes;
  bool timestamp;

  enum PlotType
  {
    P_TEMP,
    P_FLUX
  };

  PlotType plotType;

  enum FluxDir
  {
    D_M,
    D_X,
    D_Y,
    D_Z
  };
  FluxDir fluxDir;

  enum ColorScheme
  {
    C_CMR,
    C_JET,
    C_NONE
  };
  ColorScheme colorScheme;

  enum Format
  {
    F_PNG,
    F_TEX
  };
  Format format;

  enum OutputUnits
  {
    IP,
    SI
  };
  OutputUnits outputUnits;

  double minValue;
  double maxValue;

  int numberOfContours;
  std::string contourColor;

  std::pair<double, double> xRange;
  std::pair<double, double> yRange;
  std::pair<double, double> zRange;

};

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
  SnapshotSettings snapshotSettings;
  std::vector<Block> blocks;
  std::vector<Surface> surfaces;
  SliceType sliceType;

  mglData TDat;
  std::size_t iMin, iMax, jMin, jMax, kMin, kMax;

  double distanceUnitConversion;

  double tStart, tEnd;
  double nextPlotTime;
  GroundPlot(SnapshotSettings &snapshotSettings, Domain &domain, Foundation &foundation);
  void createFrame(std::string timeStamp="");
  bool makeNewFrame(double t);
};

}
#endif /* GROUNDPLOT_H_ */
