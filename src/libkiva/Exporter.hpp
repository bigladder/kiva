/* Copyright (c) 2012-2022 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#ifndef Exporter_HPP
#define Exporter_HPP

#include <filesystem>
#include <nlohmann/json.hpp>

#include "Instance.hpp"
#include "libkiva_export.h"

namespace Kiva {

struct ExportItem {
  const Instance *instance;
  // Written unchanged as the instance's "metadata" (e.g. to identify instances created
  // automatically by a calling program); omitted if null
  nlohmann::ordered_json metadata;
};

// Writes each item's geometry and recorded snapshots. Only snapshots whose settings.tag equals tag
// are written; an empty tag writes all snapshots. Items with no matching snapshots are still
// written (geometry only). Geometry is read from each instance's live Ground when called.
//
// Memory: recorded results are held in memory by each Instance until written, at 8 bytes per cell
// per record (e.g. a 50,000-cell subdomain recorded hourly for a year is ~3.5 GB). The write*
// functions and getSnapshotsCbor stream the output, so they add at most one instance's geometry and
// one record on top of that (getSnapshotsCbor also holds the encoded bytes). getSnapshotsJson
// builds the whole document as a JSON tree, which takes several times the raw results size.
LIBKIVA_EXPORT nlohmann::ordered_json getSnapshotsJson(const std::vector<ExportItem> &items,
                                                       const std::string &tag = "");
LIBKIVA_EXPORT std::vector<uint8_t> getSnapshotsCbor(const std::vector<ExportItem> &items,
                                                     const std::string &tag = "");
LIBKIVA_EXPORT void writeSnapshotsJson(const std::filesystem::path &outputPath,
                                       const std::vector<ExportItem> &items,
                                       const std::string &tag = "");
LIBKIVA_EXPORT void writeSnapshotsCbor(const std::filesystem::path &outputPath,
                                       const std::vector<ExportItem> &items,
                                       const std::string &tag = "");

} // namespace Kiva
#endif // Exporter_HPP
