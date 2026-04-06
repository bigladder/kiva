/* Copyright (c) 2012-2022 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#ifndef Exporter_CPP
#define Exporter_CPP

#include "Exporter.hpp"

namespace Kiva {

Exporter::Exporter() {}

void Exporter::initialize(const Foundation &foundation, const Domain &domain,
                          const std::filesystem::path &inputPath) {
  jExport["metadata"] = createMetadata(inputPath);

  for (const Surface &surface : foundation.surfaces) {
    nlohmann::ordered_json jSurface = createSurface(surface);
    jExport["surfaces"].push_back(jSurface);
  }

  for (const Block &block : foundation.blocks) {
    nlohmann::ordered_json jBlock = createBlock(block, foundation);
    jExport["blocks"].push_back(jBlock);
  }

  jExport["mesh"] = createMesh(domain);

  for (const std::shared_ptr<Cell> cellPtr : domain.cell) {
    if (cellPtr) {
      nlohmann::ordered_json jCell = createCell(*cellPtr, foundation);
      jExport["cells"].push_back(jCell);
    }
  }
}

void Exporter::addSnapshot(const Subdomain &subdomain) {
  nlohmann::ordered_json jSnapshot = createSnapshot(subdomain);
  jExport["snapshots"].push_back(jSnapshot);
}

void Exporter::addSnapshotResults(const std::size_t &snapshotIndex, const Subdomain &subdomain,
                                  const boost::posix_time::ptime &timestamp) {
  nlohmann::ordered_json jSnapshotResults = createSnapshotResults(subdomain, timestamp);
  jExport["snapshots"][snapshotIndex]["results"].push_back(jSnapshotResults);
}

void Exporter::exportCBOR(const std::filesystem::path &outputDir,
                          const std::filesystem::path &inputPath) {
  std::ofstream file((outputDir / inputPath.filename()).replace_extension(".cbor"),
                     std::ios::out | std::ios::binary);

  if (file.is_open()) {
    std::vector<uint8_t> cbor = nlohmann::ordered_json::to_cbor(jExport);
    file.write(reinterpret_cast<const char *>(cbor.data()), cbor.size());
    file.close();
  }
}

void Exporter::exportJSON(const std::filesystem::path &outputDir,
                          const std::filesystem::path &inputPath) {
  std::ofstream file((outputDir / inputPath.filename()).replace_extension(".json"), std::ios::out);

  if (file.is_open()) {
    file << std::setw(4) << jExport << std::endl;
    file.close();
  }
}

nlohmann::ordered_json Exporter::createMetadata(const std::filesystem::path &inputPath) {
  nlohmann::ordered_json jMetadata;

  jMetadata["schema_author"] = "Big Ladder Software";
  jMetadata["schema_name"] = "KIVA_EXPORT";
  jMetadata["schema_version"] = "0.1.0";
  jMetadata["description"] = "Kiva model data and snapshot results for " + inputPath.string();
  jMetadata["time_of_creation"] = formatTime(boost::posix_time::microsec_clock::universal_time());

  return jMetadata;
}

nlohmann::ordered_json Exporter::createSurface(const Surface &surface) {
  nlohmann::ordered_json jSurface;

  jSurface["surface_type"] = getSurfaceType(surface.type);
  jSurface["boundary_condition_type"] = getBoundaryConditionType(surface.boundaryConditionType);
  jSurface["orientation"] = getOrientation(surface.orientation);
  jSurface["z_min"] = surface.zMin;
  jSurface["z_max"] = surface.zMax;
  jSurface["polygon"] = createPolygon(surface.polygon);

  return jSurface;
}

nlohmann::ordered_json Exporter::createBlock(const Block &block, const Foundation &foundation) {
  nlohmann::ordered_json jBlock;

  jBlock["block_type"] = getBlockType(block.blockType);
  jBlock["z_min"] = block.zMin;
  jBlock["z_max"] = block.zMax;
  jBlock["polygon"] = createPolygon(block.polygon);

  return jBlock;
}

nlohmann::ordered_json Exporter::createPolygon(const Polygon &polygon) {
  nlohmann::ordered_json jPolygon;

  jPolygon["outer"] = createRing(polygon.outer());

  for (const Ring &ring : polygon.inners()) {
    nlohmann::ordered_json jRing = createRing(ring);
    jPolygon["inners"].push_back(jRing);
  }

  return jPolygon;
}

nlohmann::ordered_json Exporter::createRing(const Ring &ring) {
  nlohmann::ordered_json jRing;

  for (const Point &point : ring) {
    jRing.push_back({point.get<0>(), point.get<1>()});
  }

  return jRing;
}

nlohmann::ordered_json Exporter::createMesh(const Domain &domain) {
  nlohmann::ordered_json jMesh;

  std::vector<std::tuple<std::string, int>> axes = {{"x", 0}, {"y", 1}, {"z", 2}};
  for (const auto &[axis, index] : axes) {
    jMesh[axis] = domain.mesh[index].dividers;
  }

  return jMesh;
}

nlohmann::ordered_json Exporter::createCell(const Cell &cell, const Foundation &foundation) {
  nlohmann::ordered_json jCell;

  jCell["cell_type"] = getCellType(cell.cellType);
  jCell["density"] = cell.density;
  jCell["specific_heat"] = cell.specificHeat;
  jCell["conductivity"] = cell.conductivity;

  if (cell.blockPtr) {
    jCell["block_index"] = indexOf(foundation.blocks, *cell.blockPtr);
  }
  if (cell.surfacePtr) {
    jCell["surface_index"] = indexOf(foundation.surfaces, *cell.surfacePtr);
    jCell["boundary_condition_type"] =
        getBoundaryConditionType(cell.surfacePtr->boundaryConditionType);
  }

  return jCell;
}

nlohmann::ordered_json Exporter::createSnapshot(const Subdomain &subdomain) {
  nlohmann::ordered_json jSnapshot;

  jSnapshot["name"] = (std::filesystem::path(subdomain.settings.name)).filename().string();
  jSnapshot["results_type"] = getResultsType(subdomain.settings.resultsType);

  jSnapshot["x_index_min"] = subdomain.iMin;
  jSnapshot["x_index_max"] = subdomain.iMax;
  jSnapshot["y_index_min"] = subdomain.jMin;
  jSnapshot["y_index_max"] = subdomain.jMax;
  jSnapshot["z_index_min"] = subdomain.kMin;
  jSnapshot["z_index_max"] = subdomain.kMax;

  return jSnapshot;
}

nlohmann::ordered_json Exporter::createSnapshotResults(const Subdomain &subdomain,
                                                       const boost::posix_time::ptime &timestamp) {
  nlohmann::ordered_json jSnapshotResults;

  jSnapshotResults["timestamp"] = formatTime(timestamp);
  jSnapshotResults["values"] = subdomain.results;

  return jSnapshotResults;
}

std::string Exporter::getSurfaceType(const Surface::SurfaceType &surfaceType) {
  switch (surfaceType) {
  case Surface::SurfaceType::ST_SLAB_CORE:
    return "ST_SLAB_CORE";
  case Surface::SurfaceType::ST_SLAB_PERIM:
    return "ST_SLAB_PERIM";
  case Surface::SurfaceType::ST_WALL_INT:
    return "ST_WALL_INT";
  case Surface::SurfaceType::ST_WALL_EXT:
    return "ST_WALL_EXT";
  case Surface::SurfaceType::ST_WALL_TOP:
    return "ST_WALL_TOP";
  case Surface::SurfaceType::ST_GRADE:
    return "ST_GRADE";
  case Surface::SurfaceType::ST_SYMMETRY:
    return "ST_SYMMETRY";
  case Surface::SurfaceType::ST_SYMMETRY_AIR:
    return "ST_SYMMETRY_AIR";
  case Surface::SurfaceType::ST_FAR_FIELD:
    return "ST_FAR_FIELD";
  case Surface::SurfaceType::ST_FAR_FIELD_AIR:
    return "ST_FAR_FIELD_AIR";
  case Surface::SurfaceType::ST_DEEP_GROUND:
    return "ST_DEEP_GROUND";
  case Surface::SurfaceType::ST_TOP_AIR_INT:
    return "ST_TOP_AIR_INT";
  case Surface::SurfaceType::ST_TOP_AIR_EXT:
    return "ST_TOP_AIR_EXT";
  default:
    return "UNKNOWN";
  }
}

std::string
Exporter::getBoundaryConditionType(const Surface::BoundaryConditionType &boundaryConditionType) {
  switch (boundaryConditionType) {
  case Surface::BoundaryConditionType::ZERO_FLUX:
    return "ZERO_FLUX";
  case Surface::BoundaryConditionType::INTERIOR_FLUX:
    return "INTERIOR_FLUX";
  case Surface::BoundaryConditionType::EXTERIOR_FLUX:
    return "EXTERIOR_FLUX";
  case Surface::BoundaryConditionType::CONSTANT_TEMPERATURE:
    return "CONSTANT_TEMPERATURE";
  case Surface::BoundaryConditionType::INTERIOR_TEMPERATURE:
    return "INTERIOR_TEMPERATURE";
  case Surface::BoundaryConditionType::EXTERIOR_TEMPERATURE:
    return "EXTERIOR_TEMPERATURE";
  default:
    return "UNKNOWN";
  }
}

std::string Exporter::getOrientation(const Surface::Orientation &orientation) {
  switch (orientation) {
  case Surface::Orientation::X_POS:
    return "X_POS";
  case Surface::Orientation::X_NEG:
    return "X_NEG";
  case Surface::Orientation::Y_POS:
    return "Y_POS";
  case Surface::Orientation::Y_NEG:
    return "Y_NEG";
  case Surface::Orientation::Z_POS:
    return "Z_POS";
  case Surface::Orientation::Z_NEG:
    return "Z_NEG";
  default:
    return "UNKNOWN";
  }
}

std::string Exporter::getBlockType(const Block::BlockType &blockType) {
  switch (blockType) {
  case Block::BlockType::SOLID:
    return "SOLID";
  case Block::BlockType::INTERIOR_AIR:
    return "INTERIOR_AIR";
  case Block::BlockType::EXTERIOR_AIR:
    return "EXTERIOR_AIR";
  default:
    return "UNKNOWN";
  }
}

std::string Exporter::getCellType(const CellType &cellType) {
  switch (cellType) {
  case CellType::EXTERIOR_AIR:
    return "EXTERIOR_AIR";
  case CellType::INTERIOR_AIR:
    return "INTERIOR_AIR";
  case CellType::NORMAL:
    return "NORMAL";
  case CellType::BOUNDARY:
    return "BOUNDARY";
  case CellType::ZERO_THICKNESS:
    return "ZERO_THICKNESS";
  default:
    return "UNKNOWN";
  }
}

std::string Exporter::getResultsType(const SubdomainSettings::ResultsType &resultsType) {
  switch (resultsType) {
  case SubdomainSettings::ResultsType::TEMPERATURE:
    return "TEMPERATURE";
  case SubdomainSettings::ResultsType::HEAT_FLUX:
    return "HEAT_FLUX";
  default:
    return "UNKNOWN";
  }
}

std::string Exporter::formatTime(const boost::posix_time::ptime &time) {
  std::ostringstream formattedTime;

  formattedTime << std::put_time(&to_tm(time), "%Y-%m-%dT%H:%MZ");

  return formattedTime.str();
}

template <typename T> int Exporter::indexOf(const std::vector<T> &vector, const T &target) {
  int index = -1;

  const auto match = std::find_if(vector.begin(), vector.end(),
                                  [&target](const T &item) { return &item == &target; });

  if (match != vector.end()) {
    index = std::distance(vector.begin(), match);
  }

  return index;
}

} // namespace Kiva

#endif
