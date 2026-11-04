/* Copyright (c) 2012-2022 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#include "fixtures/exporter-fixture.hpp"

#include "Exporter.hpp"

using namespace Kiva;

TEST_F(ExporterFixture, export_basic) {
  Exporter exporter;
  exporter.addInstance(*ground, settings);

  for (boost::posix_time::ptime t = startTime; t < endTime; t += timestep) {
    bcs.outdoorTemp = 273 + dbt[(t - startTime).hours() % 24];

    ground->calculate(bcs, timestepSeconds);
    ground->calculateSurfaceAverages();

    exporter.addResults(*ground, t);
  }

  std::vector<uint8_t> cbor = exporter.getCbor();
  TestCbor(cbor, 1);

  nlohmann::ordered_json json = exporter.getJson();
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

TEST_F(ExporterFixture, streamed_export_matches_json) {
  Exporter exporter;
  nlohmann::ordered_json metadata = {{"floor", "Floor1"}, {"depth", 0.5}};
  exporter.addInstance(*instances[0].ground, settings, metadata);
  exporter.addInstance(*instances[1].ground, {}); // no metadata or snapshots

  for (boost::posix_time::ptime t = startTime; t < endTime; t += timestep) {
    bcs.outdoorTemp = 273 + dbt[(t - startTime).hours() % 24];
    instances[0].ground->calculate(bcs, timestepSeconds);
    instances[0].ground->calculateSurfaceAverages();
    exporter.addResults(*instances[0].ground, t);
  }

  nlohmann::ordered_json expected = withoutCreationTime(exporter.getJson());
  EXPECT_EQ(expected["instances"][0]["metadata"], metadata);
  EXPECT_FALSE(expected["instances"][1].contains("metadata"));
  EXPECT_FALSE(expected["instances"][1].contains("snapshots"));

  EXPECT_EQ(withoutCreationTime(nlohmann::ordered_json::from_cbor(exporter.getCbor())), expected);

  std::filesystem::path directory = std::filesystem::temp_directory_path();
  std::filesystem::path cborPath = directory / "kiva_streamed_export.cbor";
  std::filesystem::path jsonPath = directory / "kiva_streamed_export.json";
  exporter.writeCbor(cborPath);
  exporter.writeJson(jsonPath);

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
  Exporter exporter;
  nlohmann::ordered_json expected = withoutCreationTime(exporter.getJson());
  EXPECT_EQ(withoutCreationTime(nlohmann::ordered_json::from_cbor(exporter.getCbor())), expected);
}

TEST_F(ExporterFixture, export_aggregator) {
  Aggregator aggregator;
  for (std::size_t i = 0; i < instances.size(); i++) {
    aggregator.add_instance(instances[i].ground.get(), 1, &settings, {{"index", i}});
  }

  for (Instance &instance : instances) {
    for (boost::posix_time::ptime t = startTime; t < endTime; t += timestep) {
      bcs.outdoorTemp = 273 + dbt[(t - startTime).hours() % 24];

      instance.ground->calculate(bcs, timestepSeconds);
      instance.ground->calculateSurfaceAverages();

      aggregator.add_export_results(instance.ground.get(), t);
    }
  }

  std::vector<uint8_t> cbor = aggregator.get_export_cbor();
  TestCbor(cbor, 2);

  nlohmann::ordered_json json = aggregator.get_export_json();
  EXPECT_EQ(json["instances"][1]["metadata"]["index"], 1);
}
