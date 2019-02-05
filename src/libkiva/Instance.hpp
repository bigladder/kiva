/* Copyright (c) 2012-2018 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#ifndef Instance_HPP
#define Instance_HPP

#include <memory>

#include "Ground.hpp"

namespace Kiva {

class Instance {
public:
  Instance() = default;
  Instance(Foundation fnd);
  void create();
  std::shared_ptr<Ground> ground;
  std::shared_ptr<Foundation> foundation;

};

} // namespace Kiva

#endif
