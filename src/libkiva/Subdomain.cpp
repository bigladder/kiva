/* Copyright (c) 2012-2022 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#ifndef Subdomain_CPP
#define Subdomain_CPP

#include "Subdomain.hpp"

namespace Kiva {

SubdomainSettings::SubdomainSettings(ResultsType resultsType, boost::posix_time::ptime startTime,
                                     boost::posix_time::ptime endTime,
                                     boost::posix_time::time_duration frequency)
    : resultsType(resultsType), startTime(startTime), endTime(endTime), frequency(frequency) {}

void SubdomainSettings::setRange(RangeType rangeType, double min, double max) {
  std::pair range(min, max);
  switch (rangeType) {
  case X:
    xRange = range;
    xRangeSet = true;
    break;
  case Y:
    yRange = range;
    yRangeSet = true;
    break;
  case Z:
    zRange = range;
    zRangeSet = true;
    break;
  }
}

Subdomain::Subdomain(SubdomainSettings &settings, Ground &ground) : settings(settings) {
  if (ground.foundation.numberOfDimensions == 3) {
    if (!settings.xRangeSet) {
      settings.xRange.first = ground.domain.mesh[0].dividers[0];
      settings.xRange.second = ground.domain.mesh[0].dividers[ground.nX];
    }
    if (!settings.yRangeSet) {
      settings.yRange.first = ground.domain.mesh[1].dividers[0];
      settings.yRange.second = ground.domain.mesh[1].dividers[ground.nY];
    }
    if (!settings.zRangeSet) {

      if (!settings.xRangeSet && !settings.yRangeSet) {
        settings.zRange.first = 0;
        settings.zRange.second = 0;
      } else {
        settings.zRange.first = ground.domain.mesh[2].dividers[0];
        settings.zRange.second = ground.domain.mesh[2].dividers[ground.nZ];
      }
    }
  } else if (ground.foundation.numberOfDimensions == 2) {
    if (!settings.xRangeSet) {
      settings.xRange.first = ground.domain.mesh[0].dividers[0];
      settings.xRange.second = ground.domain.mesh[0].dividers[ground.nX];
    }
    if (!settings.yRangeSet) {
      settings.yRange.first = 0.5;
      settings.yRange.second = 0.5;
    }
    if (!settings.zRangeSet) {
      settings.zRange.first = ground.domain.mesh[2].dividers[0];
      settings.zRange.second = ground.domain.mesh[2].dividers[ground.nZ];
    }
  } else {
    if (!settings.xRangeSet) {
      settings.xRange.first = 0.5;
      settings.xRange.second = 0.5;
    }
    if (!settings.yRangeSet) {
      settings.yRange.first = 0.5;
      settings.yRange.second = 0.5;
    }
    if (!settings.zRangeSet) {
      settings.zRange.first = ground.domain.mesh[2].dividers[0];
      settings.zRange.second = ground.domain.mesh[2].dividers[ground.nZ];
    }
  }

  if (isEqual(settings.xRange.first, settings.xRange.second)) {
    iMin = ground.domain.mesh[0].getNearestIndex(settings.xRange.first);
    iMax = ground.domain.mesh[0].getNearestIndex(settings.xRange.second);
  } else {
    iMin = ground.domain.mesh[0].getPreviousIndex(settings.xRange.first);
    iMax = ground.domain.mesh[0].getNextIndex(settings.xRange.second);

    // Check for exact match
    if (isEqual(ground.domain.mesh[0].centers[iMin + 1], settings.xRange.first)) {
      iMin += 1;
    }
    if (isEqual(ground.domain.mesh[0].centers[iMax - 1], settings.xRange.second)) {
      iMax -= 1;
    }
  }
  if (isEqual(settings.yRange.first, settings.yRange.second)) {
    jMin = ground.domain.mesh[1].getNearestIndex(settings.yRange.first);
    jMax = ground.domain.mesh[1].getNearestIndex(settings.yRange.second);
  } else {
    jMin = ground.domain.mesh[1].getPreviousIndex(settings.yRange.first);
    jMax = ground.domain.mesh[1].getNextIndex(settings.yRange.second);

    // Check for exact match
    if (isEqual(ground.domain.mesh[1].centers[jMin + 1], settings.yRange.first)) {
      jMin += 1;
    }
    if (isEqual(ground.domain.mesh[1].centers[jMax - 1], settings.yRange.second)) {
      jMax -= 1;
    }
  }
  if (isEqual(settings.zRange.first, settings.zRange.second)) {
    kMin = ground.domain.mesh[2].getNearestIndex(settings.zRange.first);
    kMax = ground.domain.mesh[2].getNearestIndex(settings.zRange.second);
  } else {
    kMin = ground.domain.mesh[2].getPreviousIndex(settings.zRange.first);
    kMax = ground.domain.mesh[2].getNextIndex(settings.zRange.second);

    // Check for exact match
    if (isEqual(ground.domain.mesh[2].centers[kMin + 1], settings.zRange.first)) {
      kMin += 1;
    }
    if (isEqual(ground.domain.mesh[2].centers[kMax - 1], settings.zRange.second)) {
      kMax -= 1;
    }
  }

  iN = iMax - iMin + 1;
  jN = jMax - jMin + 1;
  kN = kMax - kMin + 1;
  results = std::vector<double>(iN * jN * kN, 0);

  nextResultsInterval = settings.startTime;
}

bool Subdomain::isNextResultsInterval(boost::posix_time::ptime &timestamp) {
  return (timestamp >= nextResultsInterval) && (timestamp >= settings.startTime) &&
         (timestamp <= settings.endTime);
}

std::size_t Subdomain::getResultsIndex(std::size_t i, std::size_t j, std::size_t k) {
  return (i - iMin) + iN * (j - jMin) + iN * jN * (k - kMin);
}

void Subdomain::updateResults(Ground &ground) {
  for (std::size_t k = kMin; k <= kMax; k++) {
    for (std::size_t j = jMin; j <= jMax; j++) {
      for (std::size_t i = iMin; i <= iMax; i++) {
        std::size_t results_index = getResultsIndex(i, j, k);
        std::size_t domain_index = ground.domain.getIndex(i, j, k);
        if (settings.resultsType == SubdomainSettings::ResultsType::TEMPERATURE) {
          results[results_index] = ground.TNew[domain_index];
        } else {
          std::array<double, 3> Qflux = ground.calculateHeatFlux(domain_index);
          double Qx = Qflux[0];
          double Qy = Qflux[1];
          double Qz = Qflux[2];
          results[results_index] = sqrt(Qx * Qx + Qy * Qy + Qz * Qz);
        }
      }
    }
  }

  nextResultsInterval += settings.frequency;
}

} // namespace Kiva

#endif
