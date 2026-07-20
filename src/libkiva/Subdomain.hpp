/* Copyright (c) 2012-2022 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#ifndef Subdomain_HPP
#define Subdomain_HPP

#include <boost/date_time/posix_time/posix_time.hpp>

#include "Ground.hpp"
#include "libkiva_export.h"

namespace Kiva {

class LIBKIVA_EXPORT SubdomainSettings {
public:
  enum ResultsType { TEMPERATURE, HEAT_FLUX };
  ResultsType resultsType;

  boost::posix_time::ptime startTime;
  boost::posix_time::ptime endTime;
  boost::posix_time::time_duration frequency;

  std::pair<double, double> xRange;
  std::pair<double, double> yRange;
  std::pair<double, double> zRange;

  bool xRangeSet = false;
  bool yRangeSet = false;
  bool zRangeSet = false;

  enum RangeType { X, Y, Z };

  SubdomainSettings(ResultsType resultsType, boost::posix_time::ptime startTime,
                    boost::posix_time::ptime endTime, boost::posix_time::time_duration frequency);
  void setRange(RangeType rangeType, double min, double max);
};

class LIBKIVA_EXPORT Subdomain {
public:
  SubdomainSettings settings;

  std::size_t iMin, iMax, iN;
  std::size_t jMin, jMax, jN;
  std::size_t kMin, kMax, kN;
  std::vector<double> results;

  boost::posix_time::ptime nextResultsInterval;

  Subdomain(SubdomainSettings &settings, Ground &ground);
  bool isNextResultsInterval(const boost::posix_time::ptime &timestamp);
  std::size_t getResultsIndex(std::size_t i, std::size_t j, std::size_t k);
  void updateResults(Ground &ground);
};

} // namespace Kiva
#endif /* Subdomain_HPP */
