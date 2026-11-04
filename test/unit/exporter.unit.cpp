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

TEST_F(ExporterFixture, export_aggregator) {
  Aggregator aggregator;
  for (Instance &instance : instances) {
    aggregator.add_instance(instance.ground.get(), 1, &settings);
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
}
