/* Copyright (c) 2012-2022 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#ifndef Instance_HPP
#define Instance_HPP

#include <memory>

#include "Ground.hpp"
#include "Subdomain.hpp"

namespace Kiva {

class LIBKIVA_EXPORT Instance {
public:
  Instance() = default;
  Instance(Foundation fnd);
  void create();
  void calculate(double ts = 0.0);
  void calculate_surface_averages();

  // Requires create() first, since the subdomain is resolved against the built domain
  void addSnapshot(const SubdomainSettings &settings);
  // Records results for every snapshot whose next results interval has been reached
  void recordSnapshots(const boost::posix_time::ptime &timestamp);
  const std::vector<Subdomain> &getSnapshots() const;
  // One list of records per snapshot, parallel to getSnapshots()
  const std::vector<std::vector<SnapshotRecord>> &getSnapshotRecords() const;

  std::shared_ptr<Ground> ground;
  std::shared_ptr<Foundation> foundation;
  std::shared_ptr<BoundaryConditions> bcs;

private:
  std::vector<Subdomain> snapshots;
  std::vector<std::vector<SnapshotRecord>> snapshotRecords;
};

} // namespace Kiva

#endif
