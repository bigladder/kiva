/* Copyright (c) 2012-2022 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#ifndef Exporter_CPP
#define Exporter_CPP

#include "Exporter.hpp"

namespace Kiva {

Exporter::Exporter() {}

void Exporter::initialize(const Foundation &foundation, const Domain &domain) {
  jExport["metadata"] = createMetadata();

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

void Exporter::addSnapshot(const GroundPlot &groundPlot) {
  nlohmann::ordered_json jSnapshot = createSnapshot(groundPlot);
  jExport["snapshots"].push_back(jSnapshot);
}

void Exporter::addSnapshotResults(const std::size_t &snapshotIndex, const GroundPlot &groundPlot,
                                  const boost::posix_time::ptime &timestamp) {
  nlohmann::ordered_json jSnapshotResults = createSnapshotResults(groundPlot, timestamp);
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

nlohmann::ordered_json Exporter::createMetadata() {
  nlohmann::ordered_json jMetadata;

  jMetadata["schema_author"] = "Big Ladder Software";
  jMetadata["schema_name"] = "KIVA_EXPORT";
  jMetadata["schema_version"] = "0.1.0";
  jMetadata["author"] = "Big Ladder Software";
  jMetadata["description"] = "Kiva model data and snapshot results";
  jMetadata["time_of_creation"] = formatTime(boost::posix_time::microsec_clock::universal_time());
  jMetadata["version"] = "1.0.0";

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

nlohmann::ordered_json Exporter::createSnapshot(const GroundPlot &groundPlot) {
  nlohmann::ordered_json jSnapshot;

  jSnapshot["directory"] =
      (std::filesystem::path(groundPlot.snapshotSettings.dir)).filename().string();
  jSnapshot["plot_type"] = getPlotType(groundPlot.snapshotSettings.plotType);
  jSnapshot["units"] =
      getUnits(groundPlot.snapshotSettings.plotType, groundPlot.snapshotSettings.outputUnits);
  groundPlot.snapshotSettings.xRange;

  jSnapshot["x_min"] = groundPlot.iMin;
  jSnapshot["x_max"] = groundPlot.iMax;
  jSnapshot["y_min"] = groundPlot.jMin;
  jSnapshot["y_max"] = groundPlot.jMax;
  jSnapshot["z_min"] = groundPlot.kMin;
  jSnapshot["z_max"] = groundPlot.kMax;

  return jSnapshot;
}

nlohmann::ordered_json Exporter::createSnapshotResults(const GroundPlot &groundPlot,
                                                       const boost::posix_time::ptime &timestamp) {
  nlohmann::ordered_json jResults;

  jResults["timestamp"] = formatTime(timestamp);
  jResults["values"] =
      std::vector<double>(groundPlot.TDat.a, groundPlot.TDat.a + groundPlot.TDat.GetNN());

  return jResults;
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

std::string Exporter::getPlotType(const SnapshotSettings::PlotType &plotType) {
  switch (plotType) {
  case SnapshotSettings::PlotType::P_TEMP:
    return "P_TEMP";
  case SnapshotSettings::PlotType::P_FLUX:
    return "P_FLUX";
  default:
    return "UNKNOWN";
  }
}

std::string Exporter::getUnits(const SnapshotSettings::PlotType &plotType,
                               const SnapshotSettings::OutputUnits &outputUnits) {
  switch (plotType) {
  case SnapshotSettings::PlotType::P_TEMP:
    return (outputUnits == SnapshotSettings::OutputUnits::IP) ? "F" : "C";
  case SnapshotSettings::PlotType::P_FLUX:
    return (outputUnits == SnapshotSettings::OutputUnits::IP) ? "W/ft2" : "W/m2";
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
