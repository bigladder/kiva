/* Copyright (c) 2012-2022 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#ifndef Aggregator_HPP
#define Aggregator_HPP

#include "Exporter.hpp"
#include "Ground.hpp"

namespace Kiva {

class LIBKIVA_EXPORT Aggregator {
public:
  Aggregator();
  Aggregator(Surface::SurfaceType st);
  void add_instance(Surface::SurfaceType st, Ground *grnd, double weight,
                    std::vector<SubdomainSettings> *settings = nullptr);
  void add_instance(Ground *grnd, double weight,
                    std::vector<SubdomainSettings> *settings = nullptr);
  std::size_t size();
  void calc_weighted_results();
  std::pair<Ground *, double> get_instance(std::size_t index);

  void add_export_results(Ground *grnd, const boost::posix_time::ptime &timestamp);
  nlohmann::ordered_json get_export_json();
  std::vector<uint8_t> get_export_cbor();
  void write_export_json(const std::filesystem::path &outputPath);
  void write_export_cbor(const std::filesystem::path &outputPath);

  struct Results {
    double hconv, hrad, qtot, qconv, qrad, Tconv, Tavg, Trad;
    void reset() { hconv = hrad = qtot = qconv = qrad = Tconv = Tavg = Trad = 0.0; }
  };

  Results results;

private:
  void validate();
  std::vector<std::pair<Ground *, double>> instances;
  Surface::SurfaceType surface_type;
  bool surface_type_set, validated;
  Exporter exporter;
};

} // namespace Kiva

#endif
