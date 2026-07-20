/* Copyright (c) 2012-2022 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#ifndef Exporter_HPP
#define Exporter_HPP

#include <boost/date_time/posix_time/posix_time.hpp>
#include <filesystem>
#include <nlohmann/json.hpp>

#include "Cell.hpp"
#include "Domain.hpp"
#include "Subdomain.hpp"
#include "libkiva_export.h"

namespace Kiva {

struct ExportInstance {
  uint64_t InstanceIndex;
  std::vector<Subdomain> Subdomains;
};

class LIBKIVA_EXPORT Exporter {
public:
  Exporter();

  void addInstance(Ground &ground, const std::vector<SubdomainSettings> &settings);
  void addResults(Ground &ground, const boost::posix_time::ptime &timestamp);

  nlohmann::ordered_json getJson();
  std::vector<uint8_t> getCbor();
  void writeJson(const std::filesystem::path &outputPath);
  void writeCbor(const std::filesystem::path &outputPath);

private:
  uint64_t nextInstanceIndex = 0;
  std::unordered_map<const Ground *, std::unique_ptr<ExportInstance>> exportInstancesMap;

  nlohmann::ordered_json jExport;

  nlohmann::ordered_json createMetadata();
  nlohmann::ordered_json createInstance(const Ground &ground);
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
