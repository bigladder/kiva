/* Copyright (c) 2012-2022 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#ifndef EXPORTER_FIXTURE_HPP_
#define EXPORTER_FIXTURE_HPP_

#include "Exporter.hpp"
#include "aggregator-fixture.hpp"

using namespace Kiva;

class ExporterFixture : public AggregatorFixture {
public:
  void SetUp() {
    AggregatorFixture::SetUp();
    ground = instances[0].ground;

    currentYear = boost::gregorian::day_clock::local_day().year();
    startTime = boost::posix_time::ptime(boost::gregorian::date(currentYear, 01, 01));
    endTime = boost::posix_time::ptime(boost::gregorian::date(currentYear, 02, 01));
    timestep = boost::posix_time::hours(1);
    timestepSeconds = static_cast<double>(timestep.total_seconds());

    SubdomainSettings tempSettings(
        SubdomainSettings::ResultsType::TEMPERATURE, startTime + boost::gregorian::days(2),
        startTime + boost::gregorian::days(4), boost::posix_time::hours(24));
    tempSettings.setRange(SubdomainSettings::RangeType::X, 0, 10);
    tempSettings.setRange(SubdomainSettings::RangeType::Y, 0, 0);
    tempSettings.setRange(SubdomainSettings::RangeType::Z, -10, 0);
    settings.push_back(tempSettings);

    SubdomainSettings fluxSettings(
        SubdomainSettings::ResultsType::HEAT_FLUX, startTime + boost::gregorian::days(2),
        startTime + boost::gregorian::days(4), boost::posix_time::hours(36));
    fluxSettings.setRange(SubdomainSettings::RangeType::X, 0, 10);
    fluxSettings.setRange(SubdomainSettings::RangeType::Y, 0, 0);
    fluxSettings.setRange(SubdomainSettings::RangeType::Z, -10, 0);
    settings.push_back(fluxSettings);
  }

  void TestCbor(const std::vector<uint8_t> &cbor, std::size_t expectedInstances) {
    nlohmann::json json;
    EXPECT_NO_THROW(json = nlohmann::json::from_cbor(cbor, true, true));

    nlohmann::json jInstances = json["instances"];
    EXPECT_EQ(jInstances.size(), expectedInstances);

    for (nlohmann::json &jInstance : jInstances) {
      EXPECT_GT(jInstance["surfaces"].size(), 0u);
      EXPECT_GT(jInstance["blocks"].size(), 0u);
      EXPECT_GT(jInstance["cells"].size(), 0u);

      nlohmann::json mesh = jInstance["mesh"];
      EXPECT_EQ(mesh.size(), 3u);
      EXPECT_GT(mesh["x"].size(), 0u);
      EXPECT_GT(mesh["y"].size(), 0u);
      EXPECT_GT(mesh["z"].size(), 0u);

      nlohmann::json snapshots = jInstance["snapshots"];
      EXPECT_EQ(snapshots.size(), 2u);

      for (nlohmann::json &snapshot : snapshots) {
        bool isTemperature = snapshot["results_type"] == "TEMPERATURE";
        EXPECT_EQ(snapshot["results"].size(), isTemperature ? 3u : 2u);

        std::size_t cellCount = (snapshot["x_index_max"].get<std::size_t>() -
                                 snapshot["x_index_min"].get<std::size_t>() + 1) *
                                (snapshot["y_index_max"].get<std::size_t>() -
                                 snapshot["y_index_min"].get<std::size_t>() + 1) *
                                (snapshot["z_index_max"].get<std::size_t>() -
                                 snapshot["z_index_min"].get<std::size_t>() + 1);
        for (nlohmann::json &result : snapshot["results"]) {
          EXPECT_EQ(result["values"].size(), cellCount);
        }
      }
    }
  }

  static nlohmann::ordered_json withoutCreationTime(nlohmann::ordered_json json) {
    json["metadata"].erase("time_of_creation");
    return json;
  }

  unsigned short currentYear;
  boost::posix_time::ptime startTime;
  boost::posix_time::ptime endTime;
  boost::posix_time::time_duration timestep;
  double timestepSeconds;

  std::vector<SubdomainSettings> settings;
};

#endif /* EXPORTER_FIXTURE_HPP_ */
