/* Copyright (c) 2012-2022 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#include "fixtures/exporter-fixture.hpp"

#include "Exporter.hpp"

using namespace Kiva;

TEST_F(ExporterFixture, export_basic) {
  Instance &instance = instances[0];
  for (const SubdomainSettings &snapshotSettings : settings) {
    instance.addSnapshot(snapshotSettings);
  }
  simulate(instance);

  std::vector<ExportItem> items = {{&instance, {}}};
  std::vector<uint8_t> cbor = getSnapshotsCbor(items);
  TestCbor(cbor, 1);

  nlohmann::ordered_json json = getSnapshotsJson(items);
  std::string expectedTimestamp =
      boost::posix_time::to_iso_extended_string(startTime + boost::gregorian::days(2))
          .substr(0, 16) +
      "Z";
  EXPECT_EQ(json["instances"][0]["snapshots"][0]["results"][0]["timestamp"], expectedTimestamp);
}

TEST_F(ExporterFixture, subdomain_values) {
  boost::posix_time::ptime t = startTime + boost::gregorian::days(2);

  Subdomain temperature(settings[0], *ground);
  temperature.updateResults(*ground, t);
  std::size_t index = ground->domain.getIndex(temperature.iMin, temperature.jMin, temperature.kMin);
  EXPECT_DOUBLE_EQ(temperature.results[0], ground->TNew[index]);

  Subdomain flux(settings[1], *ground);
  flux.updateResults(*ground, t);
  index = ground->domain.getIndex(flux.iMin, flux.jMin, flux.kMin);
  std::array<double, 3> Qflux = ground->calculateHeatFlux(index);
  EXPECT_DOUBLE_EQ(flux.results[0],
                   sqrt(Qflux[0] * Qflux[0] + Qflux[1] * Qflux[1] + Qflux[2] * Qflux[2]));
}

TEST_F(ExporterFixture, subdomain_schedule_catches_up) {
  Subdomain subdomain(settings[0], *ground);
  boost::posix_time::ptime start = settings[0].startTime;

  // Results starting late shouldn't fire every step until the schedule catches up
  boost::posix_time::ptime t = start + boost::posix_time::hours(30);
  ASSERT_TRUE(subdomain.isNextResultsInterval(t));
  subdomain.updateResults(*ground, t);
  EXPECT_EQ(subdomain.nextResultsInterval, start + boost::posix_time::hours(48));
  EXPECT_FALSE(subdomain.isNextResultsInterval(t + timestep));
}

TEST_F(ExporterFixture, subdomain_inverted_range) {
  EXPECT_DEATH(settings[0].setRange(SubdomainSettings::X, 5.0, 1.0),
               "Subdomain range minimum is greater than its maximum.");
}

TEST_F(ExporterFixture, add_snapshot_before_create) {
  Instance instance;
  EXPECT_DEATH(instance.addSnapshot(settings[0]), "Snapshot instance domain has not been built.");
}

TEST_F(ExporterFixture, streamed_export_matches_json) {
  for (const SubdomainSettings &snapshotSettings : settings) {
    instances[0].addSnapshot(snapshotSettings);
  }
  simulate(instances[0]);

  nlohmann::ordered_json metadata = {{"floor", "Floor1"}, {"depth", 0.5}};
  std::vector<ExportItem> items = {{&instances[0], metadata},
                                   {&instances[1], {}}}; // no metadata or snapshots

  nlohmann::ordered_json expected = withoutCreationTime(getSnapshotsJson(items));
  EXPECT_EQ(expected["instances"][0]["metadata"], metadata);
  EXPECT_FALSE(expected["instances"][1].contains("metadata"));
  EXPECT_FALSE(expected["instances"][1].contains("snapshots"));

  EXPECT_EQ(withoutCreationTime(nlohmann::ordered_json::from_cbor(getSnapshotsCbor(items))),
            expected);

  std::filesystem::path directory = std::filesystem::temp_directory_path();
  std::filesystem::path cborPath = directory / "kiva_streamed_export.cbor";
  std::filesystem::path jsonPath = directory / "kiva_streamed_export.json";
  writeSnapshotsCbor(cborPath, items);
  writeSnapshotsJson(jsonPath, items);

  std::ifstream cborFile(cborPath, std::ios::binary);
  EXPECT_EQ(withoutCreationTime(nlohmann::ordered_json::from_cbor(cborFile)), expected);

  // Text should match nlohmann's own pretty printing, not just parse equal
  std::ifstream jsonFile(jsonPath);
  std::string jsonText((std::istreambuf_iterator<char>(jsonFile)),
                       std::istreambuf_iterator<char>());
  std::string createdAt = nlohmann::ordered_json::parse(jsonText)["metadata"]["time_of_creation"];
  expected["metadata"]["time_of_creation"] = createdAt;
  EXPECT_EQ(jsonText, expected.dump(4) + "\n");

  cborFile.close();
  jsonFile.close();
  std::filesystem::remove(cborPath);
  std::filesystem::remove(jsonPath);
}

TEST_F(ExporterFixture, streamed_export_empty) {
  nlohmann::ordered_json expected = withoutCreationTime(getSnapshotsJson({}));
  EXPECT_EQ(withoutCreationTime(nlohmann::ordered_json::from_cbor(getSnapshotsCbor({}))), expected);
}

TEST_F(ExporterFixture, tag_filtering) {
  settings[0].tag = "temperature.cbor";
  settings[1].tag = "flux.json";
  for (const SubdomainSettings &snapshotSettings : settings) {
    instances[0].addSnapshot(snapshotSettings);
  }
  simulate(instances[0]);

  std::vector<ExportItem> items = {{&instances[0], {}}};
  std::filesystem::path directory = std::filesystem::temp_directory_path();
  std::filesystem::path cborPath = directory / settings[0].tag;
  std::filesystem::path jsonPath = directory / settings[1].tag;
  writeSnapshotsCbor(cborPath, items, settings[0].tag);
  writeSnapshotsJson(jsonPath, items, settings[1].tag);

  std::ifstream cborFile(cborPath, std::ios::binary);
  nlohmann::ordered_json temperature = nlohmann::ordered_json::from_cbor(cborFile);
  ASSERT_EQ(temperature["instances"][0]["snapshots"].size(), 1u);
  EXPECT_EQ(temperature["instances"][0]["snapshots"][0]["results_type"], "TEMPERATURE");
  EXPECT_EQ(temperature["instances"][0]["snapshots"][0]["results"].size(), 3u);

  std::ifstream jsonFile(jsonPath);
  nlohmann::ordered_json flux = nlohmann::ordered_json::parse(jsonFile);
  ASSERT_EQ(flux["instances"][0]["snapshots"].size(), 1u);
  EXPECT_EQ(flux["instances"][0]["snapshots"][0]["results_type"], "HEAT_FLUX");
  EXPECT_EQ(flux["instances"][0]["snapshots"][0]["results"].size(), 2u);

  // An empty tag writes every snapshot; an unmatched tag writes geometry only
  EXPECT_EQ(getSnapshotsJson(items)["instances"][0]["snapshots"].size(), 2u);
  nlohmann::ordered_json unmatched = getSnapshotsJson(items, "other");
  EXPECT_FALSE(unmatched["instances"][0].contains("snapshots"));
  EXPECT_GT(unmatched["instances"][0]["cells"].size(), 0u);

  cborFile.close();
  jsonFile.close();
  std::filesystem::remove(cborPath);
  std::filesystem::remove(jsonPath);
}

TEST_F(ExporterFixture, multiple_instances) {
  std::vector<ExportItem> items;
  for (std::size_t i = 0; i < instances.size(); i++) {
    for (const SubdomainSettings &snapshotSettings : settings) {
      instances[i].addSnapshot(snapshotSettings);
    }
    simulate(instances[i]);
    items.push_back({&instances[i], {{"index", i}}});
  }

  std::vector<uint8_t> cbor = getSnapshotsCbor(items);
  TestCbor(cbor, 2);

  nlohmann::ordered_json json = getSnapshotsJson(items);
  EXPECT_EQ(json["instances"][1]["metadata"]["index"], 1);
}
