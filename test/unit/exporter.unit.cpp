/* Copyright (c) 2012-2022 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#include "fixtures/exporter-fixture.hpp"

#include "Exporter.hpp"

using namespace Kiva;

TEST_F(ExporterFixture, export_basic) {
  Exporter exporter;
  exporter.addInstance(*ground, settings);

  for (boost::posix_time::ptime t = startTime; t < endTime; t += timestep) {
    bcs.outdoorTemp = 273 + dbt[timestepSeconds % 24];

    ground->calculate(bcs, timestepSeconds);
    ground->calculateSurfaceAverages();

    exporter.addResults(*ground, t);
  }

  std::vector<uint8_t> cbor = exporter.getCbor();
  TestCbor(cbor, 1);

  // exporter.writeJson("C:/Kiva/export_basic.json");
  // exporter.writeCbor("C:/Kiva/export_basic.cbor");
}

TEST_F(ExporterFixture, export_aggregator) {
  Aggregator aggregator;
  for (Instance &instance : instances) {
    aggregator.add_instance(instance.ground.get(), 1, &settings);
  }

  for (Instance &instance : instances) {
    for (boost::posix_time::ptime t = startTime; t < endTime; t += timestep) {
      bcs.outdoorTemp = 273 + dbt[timestepSeconds % 24];

      instance.ground->calculate(bcs, timestepSeconds);
      instance.ground->calculateSurfaceAverages();

      aggregator.add_export_results(instance.ground.get(), t);
    }
  }

  std::vector<uint8_t> cbor = aggregator.get_export_cbor();
  TestCbor(cbor, 2);

  // aggregator.write_export_json("C:/Kiva/export_aggregator.json");
  // aggregator.write_export_cbor("C:/Kiva/export_aggregator.cbor");
}
