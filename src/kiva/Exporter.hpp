/* Copyright (c) 2012-2022 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#ifndef Exporter_HPP
#define Exporter_HPP

#include <boost/date_time/posix_time/posix_time.hpp>
#include <nlohmann/json.hpp>

#include "Cell.hpp"
#include "Domain.hpp"
#include "GroundPlot.hpp"

namespace Kiva {

class Exporter {
public:
  Exporter();

  void initialize(const Foundation &foundation, const Domain &domain);

  void addSnapshot(const GroundPlot &groundPlot);
  void addSnapshotResults(const std::size_t &snapshotIndex, const GroundPlot &groundPlot,
                          const boost::posix_time::ptime &timestamp);

  void exportCBOR(const std::filesystem::path &outputDir, const std::filesystem::path &inputPath);
  void exportJSON(const std::filesystem::path &outputDir, const std::filesystem::path &inputPath);

private:
  nlohmann::ordered_json jExport;

  nlohmann::ordered_json createMetadata();
  nlohmann::ordered_json createSurface(const Surface &surface);
  nlohmann::ordered_json createBlock(const Block &block, const Foundation &foundation);
  nlohmann::ordered_json createPolygon(const Polygon &polygon);
  nlohmann::ordered_json createRing(const Ring &ring);
  nlohmann::ordered_json createMesh(const Domain &domain);
  nlohmann::ordered_json createCell(const Cell &cell, const Foundation &foundation);
  nlohmann::ordered_json createSnapshot(const GroundPlot &groundPlot);
  nlohmann::ordered_json createSnapshotResults(const GroundPlot &groundPlot,
                                               const boost::posix_time::ptime &timestamp);

  std::string getSurfaceType(const Surface::SurfaceType &surfaceType);
  std::string getBoundaryConditionType(const Surface::BoundaryConditionType &boundaryConditionType);
  std::string getOrientation(const Surface::Orientation &orientation);
  std::string getBlockType(const Block::BlockType &blockType);
  std::string getCellType(const CellType &cellType);
  std::string getPlotType(const SnapshotSettings::PlotType &plotType);
  std::string getUnits(const SnapshotSettings::PlotType &plotType,
                       const SnapshotSettings::OutputUnits &outputUnits);

  std::string formatTime(const boost::posix_time::ptime &time);

  template <typename T> int indexOf(const std::vector<T> &vector, const T &target);
};

} // namespace Kiva
#endif // Exporter_HPP
