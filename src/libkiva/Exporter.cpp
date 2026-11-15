/* Copyright (c) 2012-2022 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#ifndef Exporter_CPP
#define Exporter_CPP

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "Exporter.hpp"

namespace Kiva {

namespace {

// Streams an export document one element at a time. Container sizes are given up front
// because CBOR requires definite lengths to match nlohmann's own encoding.
class ExportWriter {
public:
  virtual ~ExportWriter() = default;
  virtual void beginObject(std::size_t size) = 0;
  virtual void beginArray(std::size_t size) = 0;
  virtual void end() = 0;
  virtual void key(const std::string &key) = 0;
  virtual void value(const nlohmann::ordered_json &value) = 0;
};

// CharType is char for streams and uint8_t for writing directly into a byte vector
template <typename CharType> class CborExportWriter : public ExportWriter {
public:
  explicit CborExportWriter(nlohmann::detail::output_adapter<CharType> adapter)
      : adapter(adapter), output(adapter) {}
  void beginObject(std::size_t size) override { writeHeader(5, size); }
  void beginArray(std::size_t size) override { writeHeader(4, size); }
  void end() override {}
  void key(const std::string &key) override {
    nlohmann::ordered_json::to_cbor(nlohmann::ordered_json(key), adapter);
  }
  void value(const nlohmann::ordered_json &value) override {
    nlohmann::ordered_json::to_cbor(value, adapter);
  }

private:
  // Copies of the adapter share one underlying output
  nlohmann::detail::output_adapter<CharType> adapter;
  nlohmann::detail::output_adapter_t<CharType> output;

  void put(uint8_t byte) { output->write_character(static_cast<CharType>(byte)); }

  // Shortest-form CBOR head, as nlohmann writes it
  void writeHeader(uint8_t majorType, std::uint64_t size) {
    uint8_t initialByte = static_cast<uint8_t>(majorType << 5);
    int byteCount;
    if (size < 24) {
      put(static_cast<uint8_t>(initialByte | size));
      return;
    } else if (size <= 0xFF) {
      initialByte |= 24;
      byteCount = 1;
    } else if (size <= 0xFFFF) {
      initialByte |= 25;
      byteCount = 2;
    } else if (size <= 0xFFFFFFFF) {
      initialByte |= 26;
      byteCount = 4;
    } else {
      initialByte |= 27;
      byteCount = 8;
    }
    put(initialByte);
    for (int i = byteCount - 1; i >= 0; i--) {
      put(static_cast<uint8_t>((size >> (8 * i)) & 0xFF));
    }
  }
};

// Produces the same text as nlohmann's dump(4)
class JsonExportWriter : public ExportWriter {
public:
  explicit JsonExportWriter(std::ostream &stream) : stream(stream) {}
  void beginObject(std::size_t size) override { beginContainer(size, '{', '}'); }
  void beginArray(std::size_t size) override { beginContainer(size, '[', ']'); }
  void end() override {
    Container container = containers.back();
    containers.pop_back();
    if (container.closer != '\0') {
      stream << '\n' << indentation() << container.closer;
    }
  }
  void key(const std::string &key) override {
    startElement();
    stream << nlohmann::ordered_json(key).dump() << ": ";
  }
  void value(const nlohmann::ordered_json &value) override {
    startArrayElement();
    std::string text = value.dump(indentStep);
    std::string pad = indentation();
    for (char c : text) {
      stream << c;
      if (c == '\n') {
        stream << pad;
      }
    }
  }

private:
  static constexpr int indentStep = 4;
  struct Container {
    bool isObject;
    bool isFirst;
    char closer; // '\0' for an empty container, which is closed when opened
  };
  std::ostream &stream;
  std::vector<Container> containers;

  std::string indentation() { return std::string(containers.size() * indentStep, ' '); }

  void startElement() {
    Container &container = containers.back();
    if (!container.isFirst) {
      stream << ',';
    }
    container.isFirst = false;
    stream << '\n' << indentation();
  }

  // Object members are started by key()
  void startArrayElement() {
    if (!containers.empty() && !containers.back().isObject) {
      startElement();
    }
  }

  void beginContainer(std::size_t size, char opener, char closer) {
    startArrayElement();
    if (size == 0) {
      stream << opener << closer;
      containers.push_back({opener == '{', true, '\0'});
    } else {
      stream << opener;
      containers.push_back({opener == '{', true, closer});
    }
  }
};

std::string getSurfaceType(const Surface::SurfaceType &surfaceType) {
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

std::string getBoundaryConditionType(const Surface::BoundaryConditionType &boundaryConditionType) {
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

std::string getOrientation(const Surface::Orientation &orientation) {
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

std::string getBlockType(const Block::BlockType &blockType) {
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

std::string getCellType(const CellType &cellType) {
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

std::string getResultsType(const SubdomainSettings::ResultsType &resultsType) {
  switch (resultsType) {
  case SubdomainSettings::ResultsType::TEMPERATURE:
    return "TEMPERATURE";
  case SubdomainSettings::ResultsType::HEAT_FLUX:
    return "HEAT_FLUX";
  default:
    return "UNKNOWN";
  }
}

std::string formatTime(const boost::posix_time::ptime &time) {
  std::ostringstream formattedTime;

  std::tm tm = to_tm(time);
  formattedTime << std::put_time(&tm, "%Y-%m-%dT%H:%MZ");

  return formattedTime.str();
}

template <typename T> int indexOf(const std::vector<T> &vector, const T &target) {
  int index = -1;

  const auto match = std::find_if(vector.begin(), vector.end(),
                                  [&target](const T &item) { return &item == &target; });

  if (match != vector.end()) {
    index = static_cast<int>(std::distance(vector.begin(), match));
  }

  return index;
}

nlohmann::ordered_json createMetadata() {
  nlohmann::ordered_json jMetadata;

  jMetadata["schema_author"] = "Big Ladder Software";
  jMetadata["schema_name"] = "KIVA_EXPORT";
  jMetadata["schema_version"] = "0.1.0";
  jMetadata["description"] = "Kiva model data and snapshot results";
  jMetadata["time_of_creation"] = formatTime(boost::posix_time::microsec_clock::universal_time());

  return jMetadata;
}

nlohmann::ordered_json createRing(const Ring &ring) {
  nlohmann::ordered_json jRing;

  for (const Point &point : ring) {
    jRing.push_back({point.get<0>(), point.get<1>()});
  }

  return jRing;
}

nlohmann::ordered_json createPolygon(const Polygon &polygon) {
  nlohmann::ordered_json jPolygon;

  jPolygon["outer"] = createRing(polygon.outer());

  for (const Ring &ring : polygon.inners()) {
    nlohmann::ordered_json jRing = createRing(ring);
    jPolygon["inners"].push_back(jRing);
  }

  return jPolygon;
}

nlohmann::ordered_json createSurface(const Surface &surface) {
  nlohmann::ordered_json jSurface;

  jSurface["surface_type"] = getSurfaceType(surface.type);
  jSurface["boundary_condition_type"] = getBoundaryConditionType(surface.boundaryConditionType);
  jSurface["orientation"] = getOrientation(surface.orientation);
  jSurface["z_min"] = surface.zMin;
  jSurface["z_max"] = surface.zMax;
  jSurface["polygon"] = createPolygon(surface.polygon);

  return jSurface;
}

nlohmann::ordered_json createBlock(const Block &block) {
  nlohmann::ordered_json jBlock;

  jBlock["block_type"] = getBlockType(block.blockType);
  jBlock["z_min"] = block.zMin;
  jBlock["z_max"] = block.zMax;
  jBlock["polygon"] = createPolygon(block.polygon);

  return jBlock;
}

nlohmann::ordered_json createMesh(const Domain &domain) {
  nlohmann::ordered_json jMesh;

  std::vector<std::tuple<std::string, int>> axes = {{"x", 0}, {"y", 1}, {"z", 2}};
  for (const auto &[axis, index] : axes) {
    jMesh[axis] = domain.mesh[index].dividers;
  }

  return jMesh;
}

nlohmann::ordered_json createCell(const Cell &cell, const Foundation &foundation) {
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

nlohmann::ordered_json createSnapshot(const Subdomain &subdomain) {
  nlohmann::ordered_json jSnapshot;

  jSnapshot["results_type"] = getResultsType(subdomain.settings.resultsType);
  jSnapshot["x_index_min"] = subdomain.iMin;
  jSnapshot["x_index_max"] = subdomain.iMax;
  jSnapshot["y_index_min"] = subdomain.jMin;
  jSnapshot["y_index_max"] = subdomain.jMax;
  jSnapshot["z_index_min"] = subdomain.kMin;
  jSnapshot["z_index_max"] = subdomain.kMax;

  return jSnapshot;
}

nlohmann::ordered_json createSnapshotResults(const SnapshotRecord &record) {
  nlohmann::ordered_json jSnapshotResults;

  jSnapshotResults["timestamp"] = formatTime(record.timestamp);
  jSnapshotResults["values"] = record.values;

  return jSnapshotResults;
}

// Indices of the instance's snapshots written for tag
std::vector<std::size_t> getMatchingSnapshots(const Instance &instance, const std::string &tag) {
  std::vector<std::size_t> indices;
  const std::vector<Subdomain> &snapshots = instance.getSnapshots();
  for (std::size_t i = 0; i < snapshots.size(); i++) {
    if (tag.empty() || snapshots[i].settings.tag == tag) {
      indices.push_back(i);
    }
  }
  return indices;
}

// Everything in an instance except its snapshots
nlohmann::ordered_json createInstanceModel(const ExportItem &item) {
  nlohmann::ordered_json jInstance;
  const Ground &ground = *item.instance->ground;

  if (!item.metadata.is_null()) {
    jInstance["metadata"] = item.metadata;
  }

  for (const Surface &surface : ground.foundation.surfaces) {
    nlohmann::ordered_json jSurface = createSurface(surface);
    jInstance["surfaces"].push_back(jSurface);
  }

  for (const Block &block : ground.foundation.blocks) {
    nlohmann::ordered_json jBlock = createBlock(block);
    jInstance["blocks"].push_back(jBlock);
  }

  jInstance["mesh"] = createMesh(ground.domain);

  for (const auto &cellPtr : ground.domain.cell) {
    if (cellPtr) {
      nlohmann::ordered_json jCell = createCell(*cellPtr, ground.foundation);
      jInstance["cells"].push_back(jCell);
    }
  }

  return jInstance;
}

nlohmann::ordered_json createInstance(const ExportItem &item, const std::string &tag) {
  nlohmann::ordered_json jInstance = createInstanceModel(item);
  const Instance &instance = *item.instance;

  for (std::size_t i : getMatchingSnapshots(instance, tag)) {
    nlohmann::ordered_json jSnapshot = createSnapshot(instance.getSnapshots()[i]);
    // A snapshot has no "results" key until its first results are recorded
    for (const SnapshotRecord &record : instance.getSnapshotRecords()[i]) {
      jSnapshot["results"].push_back(createSnapshotResults(record));
    }
    jInstance["snapshots"].push_back(jSnapshot);
  }

  return jInstance;
}

void writeInstance(ExportWriter &writer, const ExportItem &item, const std::string &tag) {
  nlohmann::ordered_json jInstance = createInstanceModel(item);
  const Instance &instance = *item.instance;
  std::vector<std::size_t> snapshotIndices = getMatchingSnapshots(instance, tag);
  bool hasSnapshots = !snapshotIndices.empty();

  writer.beginObject(jInstance.size() + (hasSnapshots ? 1 : 0));
  for (auto it = jInstance.begin(); it != jInstance.end(); ++it) {
    writer.key(it.key());
    writer.value(it.value());
  }

  if (hasSnapshots) {
    writer.key("snapshots");
    writer.beginArray(snapshotIndices.size());
    for (std::size_t i : snapshotIndices) {
      nlohmann::ordered_json jSnapshot = createSnapshot(instance.getSnapshots()[i]);
      const std::vector<SnapshotRecord> &records = instance.getSnapshotRecords()[i];

      writer.beginObject(jSnapshot.size() + (records.empty() ? 0 : 1));
      for (auto it = jSnapshot.begin(); it != jSnapshot.end(); ++it) {
        writer.key(it.key());
        writer.value(it.value());
      }
      if (!records.empty()) {
        writer.key("results");
        writer.beginArray(records.size());
        for (const SnapshotRecord &record : records) {
          writer.value(createSnapshotResults(record));
        }
        writer.end();
      }
      writer.end();
    }
    writer.end();
  }
  writer.end();
}

// Must produce the same document as getSnapshotsJson
void writeExport(ExportWriter &writer, const std::vector<ExportItem> &items,
                 const std::string &tag) {
  writer.beginObject(items.empty() ? 1 : 2);
  writer.key("metadata");
  writer.value(createMetadata());
  if (!items.empty()) {
    writer.key("instances");
    writer.beginArray(items.size());
    for (const ExportItem &item : items) {
      writeInstance(writer, item, tag);
    }
    writer.end();
  }
  writer.end();
}

} // namespace

nlohmann::ordered_json getSnapshotsJson(const std::vector<ExportItem> &items,
                                        const std::string &tag) {
  nlohmann::ordered_json jExport;

  jExport["metadata"] = createMetadata();
  for (const ExportItem &item : items) {
    jExport["instances"].push_back(createInstance(item, tag));
  }

  return jExport;
}

std::vector<uint8_t> getSnapshotsCbor(const std::vector<ExportItem> &items,
                                      const std::string &tag) {
  std::vector<uint8_t> bytes;
  CborExportWriter<uint8_t> writer(bytes);
  writeExport(writer, items, tag);
  return bytes;
}

void writeSnapshotsJson(const std::filesystem::path &outputPath,
                        const std::vector<ExportItem> &items, const std::string &tag) {
  std::ofstream file(outputPath, std::ios::out);

  if (file.is_open()) {
    JsonExportWriter writer(file);
    writeExport(writer, items, tag);
    file << std::endl;
    file.close();
  } else {
    showMessage(MSG_ERR, "Could not open file for JSON export: " + outputPath.string());
  }
}

void writeSnapshotsCbor(const std::filesystem::path &outputPath,
                        const std::vector<ExportItem> &items, const std::string &tag) {
  std::ofstream file(outputPath, std::ios::out | std::ios::binary);

  if (file.is_open()) {
    CborExportWriter<char> writer(file);
    writeExport(writer, items, tag);
    file.close();
  } else {
    showMessage(MSG_ERR, "Could not open file for CBOR export: " + outputPath.string());
  }
}

} // namespace Kiva

#endif
