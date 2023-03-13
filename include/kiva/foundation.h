/* Copyright (c) 2012-2022 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#pragma once

#include <memory>
#include <vector>

namespace Kiva {

class FoundationPrivate;

struct Material {
  Material();
  Material(double conductivity, double density, double specific_heat)
      : conductivity(conductivity), density(density), specific_heat(specific_heat) {}
  double conductivity{0.};  // [W/m-K] conductivity
  double density{0.};       // [kg/m3] density
  double specific_heat{0.}; // [J/kg-K] specific heat
};

struct Layer {
  Material material;
  double thickness; // [m] thickness
};

struct SurfaceProperties {
  SurfaceProperties() : emissivity(0.8), absorptivity(0.8), roughness(0.00208) {}
  SurfaceProperties(double emissivity, double absorptivity, double roughness)
      : emissivity(emissivity), absorptivity(absorptivity), roughness(roughness) {}
  double emissivity, absorptivity, roughness;
};

struct Construction {
  SurfaceProperties interior;
  std::vector<Layer> layers;
  double total_width() const {
    double width{0.};
    for (auto layer : layers) {
      width += layer.thickness;
    }
    return width;
  };
  double total_resistance() const {
    double resistance{0.};
    for (auto layer : layers) {
      resistance += (layer.thickness / layer.material.conductivity);
    }
    return resistance;
  };
};

struct Wall : Construction {
  SurfaceProperties exterior;
  double height_above_grade{0.}; // [m]
  double depth_below_slab{0.};   // [m]
};

struct MaterialBlock {
  Material material;

  enum class XReference { symmetry, wall_interior, wall_center, wall_exterior, far_field };
  enum class ZReference { wall_top, grade, slab_top, slab_bottom, wall_bottom, deep_ground };

  struct Point {
    XReference x_reference{XReference::wall_interior};
    ZReference z_reference{ZReference::wall_top};
    double x{0.}; // relative to x_reference
    double z{0.}; // relative to z_reference
  };

  Point point_1;
  Point point_2;
};

struct Foundation {
  // Describes the construction of the foundation
  std::vector<std::array<double, 2>> polygon;
  std::vector<bool> is_exposed_perimeter;
  Wall wall;
  Construction slab;
  std::vector<MaterialBlock> material_blocks;
  double perimeter_surface_width{0.};

private:
  bool has_wall{true};
  bool has_slab{true};
  bool use_detailed_exposed_perimeter{true};
  bool has_perimeter_surface{false};
};

struct MeshSettings {
  double minimum_cell_dimension{0.02}; // [m]
  double maximum_near_growth_coefficient{
      1.5}; // Centered growth coefficient within the "near field" domain
  double maximum_depth_growth_coefficient{1.5};    // Growth coefficient between "near field" domain
                                                   // and deep ground boundary
  double maximum_interior_growth_coefficient{1.5}; // Growth coefficient between "near field" domain
                                                   // and the interior of the foundation
  double maximum_exterior_growth_coefficient{1.5}; // Growth coefficient between "near field" domain
                                                   // and the far field boundary
};

struct KivaSettings {
  // Physical property settings
  Material soil{Material(1.73, 1842, 419)};
  double deep_ground_depth{40.}; // [m]
  enum class DeepGroundBoundaryType { fixed_temperature, zero_flux };
  DeepGroundBoundaryType deep_ground_boundary_type{DeepGroundBoundaryType::zero_flux};
  SurfaceProperties grade;

  // Boundary condition assumptions
  double far_field_width{40.}; // [m] distance from outside of wall to the edge of the domain
  enum class WallTopBoundaryType { zero_flux, linear_temperature_difference };
  WallTopBoundaryType wall_top_boundary_type{
      WallTopBoundaryType::zero_flux}; // only changed for BESTEST cases

  // Approximations
  enum class CoordinateSystem { cartesian, cylindrical };
  CoordinateSystem coordinate_system{CoordinateSystem::cartesian};
  short unsigned int number_of_dimensions{2};
  bool use_symmetry{true};
  enum class ReductionStrategy { area_perimeter, rounded_rectangle, custom, boundary };
  ReductionStrategy reduction_strategy{ReductionStrategy::boundary};

  // Numerical settings
  MeshSettings mesh_settings;
  enum class NumericalScheme {
    alternating_direction_explicit,
    explicit_forward_difference,
    alternating_direction_implicit,
    implicit,
    crank_nicolson,
    steady_state
  };
  NumericalScheme numerical_scheme{NumericalScheme::alternating_direction_implicit};
  double f_adi{0.00001}; // Alternating Direction Implicit modified f-factor
  double solver_tolerance{1.0e-6};
  double solver_maximum_iterations{100000};
};

} // namespace Kiva