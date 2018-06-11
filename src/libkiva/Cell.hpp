/* Copyright (c) 2012-2017 Big Ladder Software LLC. All rights reserved.
* See the LICENSE file for additional terms and conditions. */

#ifndef Cell_HPP
#define Cell_HPP

#include "Foundation.hpp"
#include "Mesher.hpp"
#include "Functions.hpp"

#include <fstream>
#include <memory>
#include <numeric>

#include <Eigen/SparseCore>

namespace Kiva {


enum CellType
{
  EXTERIOR_AIR,  // 0
  INTERIOR_AIR,  // 1
  NORMAL,  // 2
  BOUNDARY,  // 3
  ZERO_THICKNESS  // 4
};

class Cell
{
public:

  Cell(const std::size_t &index, const CellType cellType,
       const std::size_t &i, const std::size_t &j, const std::size_t &k,
       const Foundation &foundation, Surface *surfacePtr, Block *blockPtr,
       Mesher *meshX, Mesher *meshY, Mesher *meshZ);
  std::size_t i, j, k, index;

  // inherent properties
  double density;
  double specificHeat;
  double conductivity;

  double volume;
  double area, r;
  double heatGain;

  // derived properties
  double cxp_c;
  double cxm_c;
  double cxp, cxm, cyp, cym, czp, czm;
  double dxp, dxm, dyp, dym, dzp, dzm;
  double kxp, kxm, kyp, kym, kzp, kzm;
  double *told_ptr;
  // organizational properties

  CellType cellType;

  Block* blockPtr;

  Surface* surfacePtr;
  Mesher *meshXptr, *meshYptr, *meshZptr;
  Cell* i_up_Ptr, *i_down_Ptr, *j_up_Ptr, *j_down_Ptr, *k_up_Ptr, *k_down_Ptr;

  void setNeighbors(std::vector<Cell> &cell_v, std::size_t stepsize_i,
                    std::size_t stepsize_j, std::size_t stepsize_k);
  void setDistances(double &dxp_in, double &dxm_in, double &dyp_in, double &dym_in,
                    double &dzp_in, double &dzm_in);
  void setConductivities();
  void getKXP(), getKXM(), getKYP(), getKYM(), getKZP(), getKZM();
  void setWhatever(int ndims, bool cylindrical);
  void setZeroThicknessCellProperties(std::vector<Cell*> pointSet);

  void calcCellADEUp(double timestep, const Foundation &foundation,
                     std::vector<double> &U, std::vector<double> &UOld);
  void calcCellADEDown(double timestep, const Foundation &foundation,
                       std::vector<double> &V, std::vector<double> &VOld);
  double calcCellExplicit(double timestep, const Foundation &foundation);
  void calcCellMatrix(Foundation::NumericalScheme scheme, double timestep, const Foundation &foundation,
                std::vector<Eigen::Triplet<double>> &tripletList, Eigen::VectorXd &b,
                std::vector<double> &a1, std::vector<double> &a2, std::vector<double> &a3, std::vector<double> &b_);
  void calcCellSteadyState(const Foundation &foundation,
                std::vector<Eigen::Triplet<double>> &tripletList, Eigen::VectorXd &b,
                std::vector<double> &a1, std::vector<double> &a2, std::vector<double> &a3,
                std::vector<double> &b_);
  void calcCellADI(int dim, const Foundation &foundation, double timestsep,
                   double &Am, double &A, double &Ap, double &bVal);

  };

//  class ExteriorAirCell : public Cell {
//
//  };
//
//  class InteriorAirCell : public Cell {
//
//  };
//
//  class BoundaryCell : public Cell {
//
//  };
//
//  class ZeroThicknessCell : public Cell {
//
//  };
//
}

#endif
