/* Copyright (c) 2012-2022 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#ifndef GROUNDPLOT_H_
#define GROUNDPLOT_H_

#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4068) // Unknown pragma
#pragma warning(disable : 4458) // Declaration hides class member
#elif defined(__GNUC__) && defined(__linux__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#endif
#include <mgl2/mgl.h>
#if defined(_MSC_VER)
#pragma warning(pop)
#elif defined(__GNUC__) && defined(__linux__)
#pragma GCC diagnostic pop
#endif

#include "../libkiva/Functions.hpp"
#include "../libkiva/Geometry.hpp"
#include "../libkiva/Ground.hpp"

namespace Kiva {

class SnapshotSettings {
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

  enum PlotType { P_TEMP, P_FLUX };

  PlotType plotType;

  enum FluxDir { D_M, D_X, D_Y, D_Z };
  FluxDir fluxDir;

  enum ColorScheme { C_CMR, C_JET, C_NONE };
  ColorScheme colorScheme;

  enum Format { F_PNG, F_TEX };
  Format format;

  enum OutputUnits { IP, SI };
  OutputUnits outputUnits;

  double minValue;
  double maxValue;

  int numberOfContours;
  std::string contourColor;

  std::pair<double, double> xRange;
  std::pair<double, double> yRange;
  std::pair<double, double> zRange;
};

struct Axis {
  std::size_t nMin;
  std::size_t nMax;
  std::size_t nN;
  Mesher mesh;
};

enum SliceType { XZ_2D, XY, XZ, YZ, Z_1D };

class GroundPlot {
private:
  mglData hDat, vDat, cDat, hGrid, vGrid, TGrid;

  double slice;
  double xMin, xMax, yMin, yMax;

  Axis hAxis;
  Axis vAxis;

  int frameNumber;

public:
  // mglGraph gr;
  SnapshotSettings snapshotSettings;
  std::vector<Block> blocks;
  std::vector<Surface> surfaces;
  SliceType sliceType;

  mglData TDat;
  std::size_t iMin, iMax, jMin, jMax, kMin, kMax;
  std::size_t nI, nJ, nK;

  double distanceUnitConversion;

  double tStart, tEnd;
  double nextPlotTime;
  GroundPlot(SnapshotSettings &snapshotSettings, Domain &domain, Foundation &foundation);
  std::size_t getPlotIndex(std::size_t i, std::size_t j, std::size_t k);
  void createFrame(Ground &ground, std::string timeStamp = "");
  bool makeNewFrame(double t);
private:
  void setupPlot(Ground &ground);

};

} // namespace Kiva
#endif /* GROUNDPLOT_H_ */
