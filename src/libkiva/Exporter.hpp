/* Copyright (c) 2012-2022 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#ifndef Exporter_HPP
#define Exporter_HPP

#include <boost/date_time/posix_time/posix_time.hpp>
#include <filesystem>
#include <nlohmann/json.hpp>
#include <unordered_map>

#include "Cell.hpp"
#include "Domain.hpp"
#include "Subdomain.hpp"
#include "libkiva_export.h"

namespace Kiva {

struct SnapshotRecord {
  boost::posix_time::ptime timestamp;
  std::vector<double> values;
};

struct ExportInstance {
  const Ground *ground;
  nlohmann::ordered_json metadata;
  std::vector<Subdomain> Subdomains;
  // One list of records per subdomain, kept as raw values until the export is built
  std::vector<std::vector<SnapshotRecord>> SnapshotRecords;
};

class ExportWriter;

// Instance geometry is serialized when the export is requested, so each Ground added must
// outlive any call to getJson/getCbor/writeJson/writeCbor.
//
// Memory: recorded results are held in memory until the export is written, at 8 bytes per cell
// per record (e.g. a 50,000-cell subdomain recorded hourly for a year is ~3.5 GB). writeJson,
// writeCbor and getCbor stream the output, so they add at most one instance's geometry and one
// record on top of that (getCbor also holds the encoded bytes). getJson builds the whole
// document as a JSON tree, which takes several times the raw results size.
class LIBKIVA_EXPORT Exporter {
public:
  // metadata, if given, is written unchanged as the instance's "metadata" (e.g. to identify
  // instances created automatically by a calling program)
  void addInstance(Ground &ground, const std::vector<SubdomainSettings> &settings,
                   const nlohmann::ordered_json &metadata = nlohmann::ordered_json());
  void addResults(Ground &ground, const boost::posix_time::ptime &timestamp);

  nlohmann::ordered_json getJson();
  std::vector<uint8_t> getCbor();
  void writeJson(const std::filesystem::path &outputPath);
  void writeCbor(const std::filesystem::path &outputPath);

private:
  std::vector<ExportInstance> exportInstances;
  std::unordered_map<const Ground *, std::size_t> exportInstanceIndices;

  void writeExport(ExportWriter &writer);
  void writeInstance(ExportWriter &writer, const ExportInstance &instance);

  nlohmann::ordered_json createMetadata();
  nlohmann::ordered_json createInstance(const ExportInstance &instance);
  // Everything in an instance except its snapshots
  nlohmann::ordered_json createInstanceModel(const ExportInstance &instance);
  nlohmann::ordered_json createSurface(const Surface &surface);
  nlohmann::ordered_json createBlock(const Block &block);
  nlohmann::ordered_json createPolygon(const Polygon &polygon);
  nlohmann::ordered_json createRing(const Ring &ring);
  nlohmann::ordered_json createMesh(const Domain &domain);
  nlohmann::ordered_json createCell(const Cell &cell, const Foundation &foundation);
  nlohmann::ordered_json createSnapshot(const Subdomain &subdomain);
  nlohmann::ordered_json createSnapshotResults(const SnapshotRecord &record);

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
