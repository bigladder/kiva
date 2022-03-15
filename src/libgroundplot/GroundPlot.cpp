/* Copyright (c) 2012-2022 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#ifndef GroundPlot_CPP
#define GroundPlot_CPP

#if __has_include(<filesystem>)
#include <filesystem>
namespace filesys = std::filesystem;
#elif __has_include(<experimental/filesystem>)
#include <experimental/filesystem>
namespace filesys = std::experimental::filesystem;
#else
#error "no filesystem support"
#endif

#include <fmt/core.h>

#include "GroundPlot.hpp"

namespace Kiva {

SnapshotSettings::SnapshotSettings()
    : size(800), frequency(129600.0), grid(false), gradients(false), contours(true),
      contourLabels(false), axes(true), timestamp(true), plotType(SnapshotSettings::P_TEMP),
      fluxDir(SnapshotSettings::D_M), colorScheme(SnapshotSettings::C_CMR),
      format(SnapshotSettings::F_PNG), outputUnits(SnapshotSettings::SI), minValue(-20.0),
      maxValue(40.0), numberOfContours(13), contourColor("H") {}

GroundPlot::GroundPlot(SnapshotSettings &snapshotSettings, Domain &domain, Foundation &foundation)
    : snapshotSettings(snapshotSettings), blocks(foundation.blocks), surfaces(foundation.surfaces) {

  filesys::remove_all(snapshotSettings.dir);
  filesys::create_directory(snapshotSettings.dir);

  frameNumber = 0;

  if (snapshotSettings.outputUnits == SnapshotSettings::IP)
    distanceUnitConversion = 3.28084;
  else
    distanceUnitConversion = 1;

  std::size_t contourLevels = snapshotSettings.numberOfContours;

  if (isEqual(snapshotSettings.xRange.first, snapshotSettings.xRange.second)) {
    iMin = domain.mesh[0].getNearestIndex(snapshotSettings.xRange.first);
    iMax = domain.mesh[0].getNearestIndex(snapshotSettings.xRange.second);
  } else {
    iMin = domain.mesh[0].getPreviousIndex(snapshotSettings.xRange.first);
    iMax = domain.mesh[0].getNextIndex(snapshotSettings.xRange.second);

    // Check for exact match
    if (isEqual(domain.mesh[0].centers[iMin + 1], snapshotSettings.xRange.first))
      iMin += 1;
    if (isEqual(domain.mesh[0].centers[iMax - 1], snapshotSettings.xRange.second))
      iMax -= 1;
  }

  if (isEqual(snapshotSettings.yRange.first, snapshotSettings.yRange.second)) {
    jMin = domain.mesh[1].getNearestIndex(snapshotSettings.yRange.first);
    jMax = domain.mesh[1].getNearestIndex(snapshotSettings.yRange.second);
  } else {
    jMin = domain.mesh[1].getPreviousIndex(snapshotSettings.yRange.first);
    jMax = domain.mesh[1].getNextIndex(snapshotSettings.yRange.second);

    // Check for exact match
    if (isEqual(domain.mesh[1].centers[jMin + 1], snapshotSettings.yRange.first))
      jMin += 1;
    if (isEqual(domain.mesh[1].centers[jMax - 1], snapshotSettings.yRange.second))
      jMax -= 1;
  }

  if (isEqual(snapshotSettings.zRange.first, snapshotSettings.zRange.second)) {
    kMin = domain.mesh[2].getNearestIndex(snapshotSettings.zRange.first);
    kMax = domain.mesh[2].getNearestIndex(snapshotSettings.zRange.second);
  } else {
    kMin = domain.mesh[2].getPreviousIndex(snapshotSettings.zRange.first);
    kMax = domain.mesh[2].getNextIndex(snapshotSettings.zRange.second);

    // Check for exact match
    if (isEqual(domain.mesh[2].centers[kMin + 1], snapshotSettings.zRange.first))
      kMin += 1;
    if (isEqual(domain.mesh[2].centers[kMax - 1], snapshotSettings.zRange.second))
      kMax -= 1;
  }

  std::size_t nI = iMax - iMin + 1;
  std::size_t nJ = jMax - jMin + 1;
  std::size_t nK = kMax - kMin + 1;

  xMin = domain.mesh[0].dividers[0];
  xMax = domain.mesh[0].dividers[domain.dim_lengths[0]];

  yMin = domain.mesh[1].dividers[0];
  yMax = domain.mesh[1].dividers[domain.dim_lengths[1]];

  if (nI == 1 && nJ == 1) {
    sliceType = Z_1D;

    hAxis.nN = nI;
    hAxis.nMin = iMin;
    hAxis.nMax = iMax;
    hAxis.mesh = domain.mesh[0];

    vAxis.nN = nK;
    vAxis.nMin = kMin;
    vAxis.nMax = kMax;
    vAxis.mesh = domain.mesh[2];
  } else if (nI == 1) {
    sliceType = YZ;

    hAxis.nN = nJ;
    hAxis.nMin = jMin;
    hAxis.nMax = jMax;
    hAxis.mesh = domain.mesh[1];

    vAxis.nN = nK;
    vAxis.nMin = kMin;
    vAxis.nMax = kMax;
    vAxis.mesh = domain.mesh[2];

    slice = domain.mesh[0].centers[iMin];
  } else if (nJ == 1) {
    if (domain.mesh[1].centers.size() == 1)
      sliceType = XZ_2D;
    else
      sliceType = XZ;

    hAxis.nN = nI;
    hAxis.nMin = iMin;
    hAxis.nMax = iMax;
    hAxis.mesh = domain.mesh[0];

    vAxis.nN = nK;
    vAxis.nMin = kMin;
    vAxis.nMax = kMax;
    vAxis.mesh = domain.mesh[2];

    slice = domain.mesh[1].centers[jMin];
  } else if (nK == 1) {
    sliceType = XY;

    hAxis.nN = nI;
    hAxis.nMin = iMin;
    hAxis.nMax = iMax;
    hAxis.mesh = domain.mesh[0];

    vAxis.nN = nJ;
    vAxis.nMin = jMin;
    vAxis.nMax = jMax;
    vAxis.mesh = domain.mesh[1];

    slice = domain.mesh[2].centers[kMin];
  }

  if (sliceType == Z_1D) {
    mglData hDatRef(vAxis.nN), hGridRef(vAxis.nN + 1), TDatRef(vAxis.nN), TGridRef(vAxis.nN + 1);
    TDat = TDatRef;
    TGrid = TGridRef;
    hDat = hDatRef;
    hGrid = hGridRef;
  } else {
    mglData hDatRef(hAxis.nN), hGridRef(hAxis.nN + 1), TDatRef(hAxis.nN, vAxis.nN),
        TGridRef(hAxis.nN + 1, vAxis.nN + 1);
    TDat = TDatRef;
    TGrid = TGridRef;
    hDat = hDatRef;
    hGrid = hGridRef;
  }

  mglData vDatRef(vAxis.nN), vGridRef(vAxis.nN + 1), cDatRef(contourLevels);

  vDat = vDatRef;
  vGrid = vGridRef;
  cDat = cDatRef;

  double min = snapshotSettings.minValue;
  double max = snapshotSettings.maxValue;
  double step = (max - min) / (contourLevels - 1);

  for (size_t n = 0; n < contourLevels; n++)
    cDat.a[n] = min + double(n) * step;

  hGrid.a[0] = hAxis.mesh.centers[hAxis.nMin] * distanceUnitConversion;

  for (size_t i = 0; i < hAxis.nN - 1; i++) {
    hDat.a[i] = hAxis.mesh.centers[i + hAxis.nMin] * distanceUnitConversion;
    hGrid.a[i + 1] = hAxis.mesh.dividers[i + hAxis.nMin + 1] * distanceUnitConversion;
  }

  hDat.a[hAxis.nN - 1] = hAxis.mesh.centers[hAxis.nMax] * distanceUnitConversion;
  hGrid.a[hAxis.nN] = hAxis.mesh.centers[hAxis.nMax] * distanceUnitConversion;

  vGrid.a[0] = vAxis.mesh.centers[vAxis.nMin] * distanceUnitConversion;

  for (size_t j = 0; j < vAxis.nN - 1; j++) {
    vDat.a[j] = vAxis.mesh.centers[j + vAxis.nMin] * distanceUnitConversion;
    vGrid.a[j + 1] = vAxis.mesh.dividers[j + vAxis.nMin + 1] * distanceUnitConversion;
  }

  vDat.a[vAxis.nN - 1] = vAxis.mesh.centers[vAxis.nMax] * distanceUnitConversion;
  vGrid.a[vAxis.nN] = vAxis.mesh.centers[vAxis.nMax] * distanceUnitConversion;

  for (size_t j = 0; j <= vAxis.nN; j++) {
    for (size_t i = 0; i <= hAxis.nN; i++) {
      TGrid.a[i + hAxis.nN * j] = 200.0;
    }
  }
}

void GroundPlot::createFrame(std::string timeStamp) {

  std::string distanceUnit;
  std::string temperatureUnit;
  std::string fluxUnit;

  if (snapshotSettings.outputUnits == SnapshotSettings::IP) {
    distanceUnit = "ft";
    temperatureUnit = "\\textdegree F";
    fluxUnit = "W/ft^2";
  } else {
    distanceUnit = "m";
    temperatureUnit = "\\textdegree C";
    fluxUnit = "W/m^2";
  }

  double hMin = hGrid.a[0];
  double hMax = hGrid.a[hAxis.nN];
  double hRange = hMax - hMin;

  double vMin = vGrid.a[0];
  double vMax = vGrid.a[vAxis.nN];
  double vRange = vMax - vMin;

  int nT = cDat.GetNN();
  double Tmin = cDat.a[0];
  double Tmax = cDat.a[nT - 1];
  double Tstep = cDat.a[1] - cDat.a[0];

  // Text properties
  double hText = 0.05;
  double vText = 0.95;

  double vTextSpacing = 0.05;

  mglGraph gr;
  gr.LoadFont("heros");

  // Plot
  gr.Clf(1, 1, 1);
  double aspect = 1.0;
  int height = snapshotSettings.size;
  int width = height * aspect;

  gr.SetSize(width, height);

  gr.SetFontSize(2.0);
  if (sliceType == Z_1D) {
    gr.SetRange('x', cDat);
  } else {
    gr.SetRange('x', hGrid);
    gr.Aspect(hRange, vRange);
  }
  gr.SetRange('y', vGrid);
  gr.SetRange('c', Tmin, Tmax);
  gr.SetRange('z', Tmin, Tmax);
  gr.SetTicks('c', Tstep, nT, Tmin);

  // Timestamp

  if (snapshotSettings.axes) {
    if (snapshotSettings.colorScheme != SnapshotSettings::C_NONE) {
      if (snapshotSettings.plotType == SnapshotSettings::P_TEMP)
        gr.Puts(0.9, 0.056, temperatureUnit.c_str(), ":AL");
      else
        gr.Puts(0.9, 0.056, fluxUnit.c_str(), ":AL");
    }
  }

  if (snapshotSettings.timestamp)
    gr.Puts(hText, vText, timeStamp.c_str(), ":AL");

  switch (sliceType) {
  case XY: {
    std::string sliceString =
        "Z = " + fmt::format("Z = {:0.2f} {}", slice * distanceUnitConversion, distanceUnit);
    if (snapshotSettings.axes)
      gr.Puts(hText, vText - vTextSpacing, sliceString.c_str(), ":AL");
  } break;
  case XZ: {
    std::string sliceString =
        "Y = " + fmt::format("Z = {:0.2f} {}", slice * distanceUnitConversion, distanceUnit);
    if (snapshotSettings.axes)
      gr.Puts(hText, vText - vTextSpacing, sliceString.c_str(), ":AL");
  } break;
  case YZ: {
    std::string sliceString =
        "X = " + fmt::format("Z = {:0.2f} {}", slice * distanceUnitConversion, distanceUnit);
    if (snapshotSettings.axes)
      gr.Puts(hText, vText - vTextSpacing, sliceString.c_str(), ":AL");
  }
  default:
    break;
  }
  gr.SetPlotFactor(1.3);

  if (snapshotSettings.axes) {
    gr.Axis("yU");
    if (sliceType != Z_1D)
      gr.Axis("x");
    if (snapshotSettings.colorScheme == SnapshotSettings::C_CMR) {
      gr.Colorbar("kUrqyw_");
    } else if (snapshotSettings.colorScheme == SnapshotSettings::C_JET) {
      gr.Colorbar("BbcyrR_");
    }
  }

  if (sliceType == Z_1D) {
    if (snapshotSettings.colorScheme == SnapshotSettings::C_CMR) {
      gr.Tens(TDat, vDat, TDat, "kUrqyw2");
    } else if (snapshotSettings.colorScheme == SnapshotSettings::C_JET) {
      gr.Tens(TDat, vDat, TDat, "BbcyrR2");
    } else {
      gr.Tens(TDat, vDat, TDat, "r2");
    }
  } else {
    if (snapshotSettings.colorScheme == SnapshotSettings::C_CMR) {
      gr.Dens(hDat, vDat, TDat, "kUrqyw");
    } else if (snapshotSettings.colorScheme == SnapshotSettings::C_JET) {
      gr.Dens(hDat, vDat, TDat, "BbcyrR");
    }

    gr.Box("k2", false);

    if (snapshotSettings.contours) {
      if (snapshotSettings.contourLabels)
        gr.Cont(cDat, hDat, vDat, TDat, (snapshotSettings.contourColor + "t").c_str());
      else
        gr.Cont(cDat, hDat, vDat, TDat, snapshotSettings.contourColor.c_str());
    }
    if (snapshotSettings.gradients)
      gr.Grad(hDat, vDat, TDat);
    if (snapshotSettings.grid)
      gr.Grid(hGrid, vGrid, TGrid, "W");

    // Draw blocks
    for (auto &b : blocks) {
      if (b.blockType != Block::SOLID) {
        continue;
      }
      switch (sliceType) {
      case XZ_2D: {
        mglPoint bl =
            mglPoint(std::min(std::max(b.xMin * distanceUnitConversion, hMin), hMax),
                     std::min(std::max(b.zMin * distanceUnitConversion, vMin), vMax), 210.0);
        mglPoint br =
            mglPoint(std::max(std::min(b.xMax * distanceUnitConversion, hMax), hMin),
                     std::min(std::max(b.zMin * distanceUnitConversion, vMin), vMax), 210.0);
        mglPoint tr =
            mglPoint(std::max(std::min(b.xMax * distanceUnitConversion, hMax), hMin),
                     std::max(std::min(b.zMax * distanceUnitConversion, vMax), vMin), 210.0);
        mglPoint tl =
            mglPoint(std::min(std::max(b.xMin * distanceUnitConversion, hMin), hMax),
                     std::max(std::min(b.zMax * distanceUnitConversion, vMax), vMin), 210.0);

        gr.Line(bl, br, "k2");
        gr.Line(br, tr, "k2");
        gr.Line(tr, tl, "k2");
        gr.Line(tl, bl, "k2");
      } break;
      case XY: {
        // Find intersection with viewing window
        Polygon view;
        view.outer().push_back(Point(hMin, vMin));
        view.outer().push_back(Point(hMin, vMax));
        view.outer().push_back(Point(hMax, vMax));
        view.outer().push_back(Point(hMax, vMin));

        MultiPolygon intersection;
        boost::geometry::intersection(view, b.polygon, intersection);

        // loop through each polygon in resulting multi_polygon
        for (std::size_t p = 0; p < intersection.size(); p++) {
          // loop through each vertex in each polygon and create a line
          std::size_t nV = intersection[p].outer().size();
          for (std::size_t v = 0; v < nV - 1; v++) {
            gr.Line(mglPoint(intersection[p].outer()[v].get<0>() * distanceUnitConversion,
                             intersection[p].outer()[v].get<1>() * distanceUnitConversion, 210.0),
                    mglPoint(intersection[p].outer()[v + 1].get<0>() * distanceUnitConversion,
                             intersection[p].outer()[v + 1].get<1>() * distanceUnitConversion,
                             210.0),
                    "k2");
          }
          gr.Line(mglPoint(intersection[p].outer()[nV - 1].get<0>() * distanceUnitConversion,
                           intersection[p].outer()[nV - 1].get<1>() * distanceUnitConversion,
                           210.0),
                  mglPoint(intersection[p].outer()[0].get<0>() * distanceUnitConversion,
                           intersection[p].outer()[0].get<1>() * distanceUnitConversion, 210.0),
                  "k2");
        }
      } break;
      case XZ: {
        // Find intersecting point pairs with slicing plane
        Line slicingPlane;
        slicingPlane.push_back(Point(xMin - EPSILON, slice));
        slicingPlane.push_back(Point(xMax + EPSILON, slice));

        MultiPoint intersection;
        boost::geometry::intersection(slicingPlane, b.polygon, intersection);

        // sort points in ascending order
        sort(intersection.begin(), intersection.end(), comparePointsX);

        for (std::size_t p = 0; p < intersection.size(); p++) {
          // Use point pairs and zmin/max to draw rectangles similar to 2D
          // case.
          double p1 = intersection[p].get<0>() * distanceUnitConversion;
          double p2 = intersection[p + 1].get<0>() * distanceUnitConversion;

          mglPoint bl =
              mglPoint(std::min(std::max(p1, hMin), hMax),
                       std::min(std::max(b.zMin * distanceUnitConversion, vMin), vMax), 210.0);
          mglPoint br =
              mglPoint(std::max(std::min(p2, hMax), hMin),
                       std::min(std::max(b.zMin * distanceUnitConversion, vMin), vMax), 210.0);
          mglPoint tr =
              mglPoint(std::max(std::min(p2, hMax), hMin),
                       std::max(std::min(b.zMax * distanceUnitConversion, vMax), vMin), 210.0);
          mglPoint tl =
              mglPoint(std::min(std::max(p1, hMin), hMax),
                       std::max(std::min(b.zMax * distanceUnitConversion, vMax), vMin), 210.0);

          gr.Line(bl, br, "k2");
          gr.Line(br, tr, "k2");
          gr.Line(tr, tl, "k2");
          gr.Line(tl, bl, "k2");

          p += 1; // skip one point, on to the next pair
        }
      } break;
      case YZ: {
        // Find intersecting point pairs with slicing plane
        Line slicingPlane;
        slicingPlane.push_back(Point(slice, yMin - EPSILON));
        slicingPlane.push_back(Point(slice, yMax + EPSILON));

        MultiPoint intersection;
        boost::geometry::intersection(slicingPlane, b.polygon, intersection);

        // sort points in ascending order
        sort(intersection.begin(), intersection.end(), comparePointsY);

        for (std::size_t p = 0; p < intersection.size(); p++) {
          // Use point pairs and zmin/max to draw rectangles similar to 2D
          // case.
          double p1 = intersection[p].get<1>() * distanceUnitConversion;
          double p2 = intersection[p + 1].get<1>() * distanceUnitConversion;

          mglPoint bl =
              mglPoint(std::min(std::max(p1, hMin), hMax),
                       std::min(std::max(b.zMin * distanceUnitConversion, vMin), vMax), 210.0);
          mglPoint br =
              mglPoint(std::max(std::min(p2, hMax), hMin),
                       std::min(std::max(b.zMin * distanceUnitConversion, vMin), vMax), 210.0);
          mglPoint tr =
              mglPoint(std::max(std::min(p2, hMax), hMin),
                       std::max(std::min(b.zMax * distanceUnitConversion, vMax), vMin), 210.0);
          mglPoint tl =
              mglPoint(std::min(std::max(p1, hMin), hMax),
                       std::max(std::min(b.zMax * distanceUnitConversion, vMax), vMin), 210.0);

          gr.Line(bl, br, "k2");
          gr.Line(br, tr, "k2");
          gr.Line(tr, tl, "k2");
          gr.Line(tl, bl, "k2");

          p += 1; // skip one point, on to the next pair
        }
      } break;
      case Z_1D: {

      } break;
      }
    }

    // Draw surfaces
    for (auto &s : surfaces) {
      switch (sliceType) {
      case XZ_2D: {
        mglPoint a =
            mglPoint(std::min(std::max(s.xMin * distanceUnitConversion, hMin), hMax),
                     std::min(std::max(s.zMin * distanceUnitConversion, vMin), vMax), 210.0);
        mglPoint b =
            mglPoint(std::max(std::min(s.xMax * distanceUnitConversion, hMax), hMin),
                     std::min(std::max(s.zMax * distanceUnitConversion, vMin), vMax), 210.0);

        gr.Line(a, b, "k2");
      } break;
      case XY: {
      } break;
      case XZ: {
        Line slicingPlane;
        slicingPlane.push_back(Point(xMin - EPSILON, slice));
        slicingPlane.push_back(Point(xMax + EPSILON, slice));

        MultiPoint intersection;
        boost::geometry::intersection(slicingPlane, s.polygon, intersection);

        // sort points in ascending order
        sort(intersection.begin(), intersection.end(), comparePointsX);

        if (s.orientation == Surface::Z_POS || s.orientation == Surface::Z_NEG) {

          for (std::size_t p = 0; p < intersection.size(); p++) {
            // Use point pairs and zmin/max to draw rectangles similar to 2D
            // case.
            double p1 = intersection[p].get<0>() * distanceUnitConversion;
            double p2 = intersection[p + 1].get<0>() * distanceUnitConversion;

            mglPoint a =
                mglPoint(std::min(std::max(p1, hMin), hMax),
                         std::min(std::max(s.zMin * distanceUnitConversion, vMin), vMax), 210.0);
            mglPoint b =
                mglPoint(std::min(std::max(p2, hMin), hMax),
                         std::max(std::min(s.zMax * distanceUnitConversion, vMax), vMin), 210.0);

            gr.Line(a, b, "k2");

            p += 1; // skip one point, on to the next pair
          }
        } else {
          for (std::size_t p = 0; p < intersection.size(); p++) {
            mglPoint a = mglPoint(
                std::min(std::max(intersection[p].get<0>() * distanceUnitConversion, hMin), hMax),
                std::min(std::max(s.zMin * distanceUnitConversion, vMin), vMax), 210.0);
            mglPoint b = mglPoint(
                std::min(std::max(intersection[p].get<0>() * distanceUnitConversion, hMin), hMax),
                std::max(std::min(s.zMax * distanceUnitConversion, vMax), vMin), 210.0);

            gr.Line(a, b, "k2");
          }
        }
      } break;
      case YZ: {
        Line slicingPlane;
        slicingPlane.push_back(Point(yMin - EPSILON, slice));
        slicingPlane.push_back(Point(yMax + EPSILON, slice));

        MultiPoint intersection;
        boost::geometry::intersection(slicingPlane, s.polygon, intersection);

        // sort points in ascending order
        sort(intersection.begin(), intersection.end(), comparePointsY);

        if (s.orientation == Surface::Z_POS || s.orientation == Surface::Z_NEG) {

          for (std::size_t p = 0; p < intersection.size(); p++) {
            // Use point pairs and zmin/max to draw rectangles similar to 2D
            // case.
            double p1 = intersection[p].get<0>() * distanceUnitConversion;
            double p2 = intersection[p + 1].get<0>() * distanceUnitConversion;

            mglPoint a =
                mglPoint(std::min(std::max(p1, hMin), hMax),
                         std::min(std::max(s.zMin * distanceUnitConversion, vMin), vMax), 210.0);
            mglPoint b =
                mglPoint(std::min(std::max(p2, hMin), hMax),
                         std::max(std::min(s.zMax * distanceUnitConversion, vMax), vMin), 210.0);

            gr.Line(a, b, "k2");

            p += 1; // skip one point, on to the next pair
          }
        } else {
          for (std::size_t p = 0; p < intersection.size(); p++) {
            mglPoint a = mglPoint(
                std::min(std::max(intersection[p].get<0>() * distanceUnitConversion, hMin), hMax),
                std::min(std::max(s.zMin * distanceUnitConversion, vMin), vMax), 210.0);
            mglPoint b = mglPoint(
                std::min(std::max(intersection[p].get<0>() * distanceUnitConversion, hMin), hMax),
                std::max(std::min(s.zMax * distanceUnitConversion, vMax), vMin), 210.0);

            gr.Line(a, b, "k2");
          }
        }
      } break;
      case Z_1D: {

      } break;
      }
    }
  }

  if (snapshotSettings.format == SnapshotSettings::F_PNG)
    gr.WritePNG((fmt::format("{}/{:04d}.png", snapshotSettings.dir, frameNumber)).c_str(), "",
                false);
  else if (snapshotSettings.format == SnapshotSettings::F_TEX)
    gr.WriteTEX((fmt::format("{}/{:04d}.tex", snapshotSettings.dir, frameNumber)).c_str());

  frameNumber += 1;
  nextPlotTime += snapshotSettings.frequency;
}

bool GroundPlot::makeNewFrame(double t) {
  if ((t >= nextPlotTime) && (t >= tStart) && (t <= tEnd))
    return true;
  else
    return false;
}

} // namespace Kiva

#endif
