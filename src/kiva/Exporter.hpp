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

  void addSnapshot(const GroundPlot &groundPlot, const Domain &domain,
                   const boost::posix_time::ptime &startTime);
  void addResults(const std::size_t &snapshotIndex);

  void exportCBOR(const std::filesystem::path &outputDir);
  void exportJSON(const std::filesystem::path &outputDir);

private:
  nlohmann::ordered_json jExport;

  nlohmann::ordered_json createMetadata();
  nlohmann::ordered_json createSurface(const Surface &surface);
  nlohmann::ordered_json createBlock(const Block &block, const Foundation &foundation);
  nlohmann::ordered_json createPolygon(const Polygon &polygon);
  nlohmann::ordered_json createRing(const Ring &ring);
  nlohmann::ordered_json createCell(const Cell &cell, const Foundation &foundation);
  nlohmann::ordered_json createSnapshot(const GroundPlot &groundPlot, const Domain &domain,
                                        const boost::posix_time::ptime &startTime);
  nlohmann::ordered_json createMesh(const GroundPlot &groundPlot, const Domain &domain);
  nlohmann::ordered_json createAxis(const std::size_t &min, const std::size_t &max,
                                    const Mesher &mesh);
  nlohmann::ordered_json createTimeInterval(const SnapshotSettings &snapshotSettings,
                                            const boost::posix_time::ptime &startTime);
  nlohmann::ordered_json createTimeSeries(const std::size_t &snapshotIndex);

  std::string getSurfaceType(const Surface::SurfaceType &surfaceType);
  std::string getBlockType(const Block::BlockType &blockType);
  std::string getCellType(const CellType &cellType);
  std::string getBoundaryConditionType(const Surface::BoundaryConditionType &boundaryConditionType);

  std::string formatTime(const boost::posix_time::ptime &time);

  template <typename T> int indexOf(const std::vector<T> &vector, const T &target);
};

} // namespace Kiva
#endif // Exporter_HPP
