/* Copyright (c) 2012-2022 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#ifndef Exporter_HPP
#define Exporter_HPP

#include <boost/date_time/posix_time/posix_time.hpp>
#include <nlohmann/json.hpp>

#include "Cell.hpp"
#include "Domain.hpp"
#include "Subdomain.hpp"

namespace Kiva {

class Exporter {
public:
  Exporter();

  void initialize(const Foundation &foundation, const Domain &domain,
                  const std::filesystem::path &inputPath);

  void addSnapshot(const Subdomain &subdomain);
  void addSnapshotResults(const std::size_t &snapshotIndex, const Subdomain &subdomain,
                          const boost::posix_time::ptime &timestamp);

  void exportCBOR(const std::filesystem::path &outputDir, const std::filesystem::path &inputPath);
  void exportJSON(const std::filesystem::path &outputDir, const std::filesystem::path &inputPath);

private:
  nlohmann::ordered_json jExport;

  nlohmann::ordered_json createMetadata(const std::filesystem::path &inputPath);
  nlohmann::ordered_json createSurface(const Surface &surface);
  nlohmann::ordered_json createBlock(const Block &block, const Foundation &foundation);
  nlohmann::ordered_json createPolygon(const Polygon &polygon);
  nlohmann::ordered_json createRing(const Ring &ring);
  nlohmann::ordered_json createMesh(const Domain &domain);
  nlohmann::ordered_json createCell(const Cell &cell, const Foundation &foundation);
  nlohmann::ordered_json createSnapshot(const Subdomain &subdomain);
  nlohmann::ordered_json createSnapshotResults(const Subdomain &subdomain,
                                               const boost::posix_time::ptime &timestamp);

  std::string getSurfaceType(const Surface::SurfaceType &surfaceType);
  std::string getBoundaryConditionType(const Surface::BoundaryConditionType &boundaryConditionType);
  std::string getOrientation(const Surface::Orientation &orientation);
  std::string getBlockType(const Block::BlockType &blockType);
  std::string getCellType(const CellType &cellType);
  std::string getResultsType(const SubdomainSettings::ResultsType &resultsType);

  std::string formatTime(const boost::posix_time::ptime &time);

  template <typename T> int indexOf(const std::vector<T> &vector, const T &target);
};

} // namespace Kiva
#endif // Exporter_HPP
