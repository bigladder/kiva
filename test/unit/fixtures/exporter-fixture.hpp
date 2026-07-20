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
    timestepSeconds = timestep.total_seconds();

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

  void TestCbor(const std::vector<uint8_t> &cbor, int expectedInstances) {
    nlohmann::json json;
    EXPECT_NO_THROW(json = nlohmann::json::from_cbor(cbor, true, true));

    nlohmann::json instances = json["instances"];
    EXPECT_EQ(instances.size(), expectedInstances);

    for (nlohmann::json &instance : instances) {
      EXPECT_GT(instance["surfaces"].size(), 0);
      EXPECT_GT(instance["blocks"].size(), 0);
      EXPECT_GT(instance["cells"].size(), 0);

      nlohmann::json mesh = instance["mesh"];
      EXPECT_EQ(mesh.size(), 3);
      EXPECT_GT(mesh["x"].size(), 0);
      EXPECT_GT(mesh["y"].size(), 0);
      EXPECT_GT(mesh["z"].size(), 0);

      nlohmann::json snapshots = instance["snapshots"];
      EXPECT_EQ(snapshots.size(), 2);

      for (nlohmann::json &snapshot : snapshots) {
        EXPECT_EQ(snapshot["results"].size(), (snapshot["results_type"] == "TEMPERATURE") ? 3 : 2);
      }
    }
  }

  int currentYear;
  boost::posix_time::ptime startTime;
  boost::posix_time::ptime endTime;
  boost::posix_time::time_duration timestep;
  int64_t timestepSeconds;

  std::vector<SubdomainSettings> settings;
};

#endif /* EXPORTER_FIXTURE_HPP_ */
