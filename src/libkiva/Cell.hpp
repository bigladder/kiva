/* Copyright (c) 2012-2017 Big Ladder Software LLC. All rights reserved.
* See the LICENSE file for additional terms and conditions. */

#ifndef Cell_HPP
#define Cell_HPP

#include "Foundation.hpp"
#include "Mesher.hpp"
#include "Functions.hpp"
#include "BoundaryConditions.hpp"
#include "Algorithms.hpp"

#include <fstream>
#include <memory>
#include <numeric>

#include <Eigen/SparseCore>

namespace Kiva {


  enum CellType {
    EXTERIOR_AIR,  // 0
    INTERIOR_AIR,  // 1
    NORMAL,  // 2
    BOUNDARY,  // 3
    ZERO_THICKNESS  // 4
  };

  class Cell {
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

    Block *blockPtr;

    Surface *surfacePtr;
    Mesher *meshXptr, *meshYptr, *meshZptr;
    std::shared_ptr<Cell> i_up_Ptr, i_down_Ptr, j_up_Ptr, j_down_Ptr, k_up_Ptr, k_down_Ptr;

    void setNeighbors(std::vector<std::shared_ptr<Cell> > &cell_v,
                      std::size_t stepsize_i, std::size_t stepsize_j, std::size_t stepsize_k,
                      std::size_t nX, std::size_t nY, std::size_t nZ);

    void setDistances(double &dxp_in, double &dxm_in, double &dyp_in, double &dym_in,
                      double &dzp_in, double &dzm_in);

    void setConductivities();

    void getKXP(), getKXM(), getKYP(), getKYM(), getKZP(), getKZM();

    void setWhatever(int ndims, bool cylindrical);

    void setZeroThicknessCellProperties(std::vector< std::shared_ptr<Cell> > pointSet);

    void calcCellADEUp(double timestep, const Foundation &foundation,
                       std::vector<double> &U, std::vector<double> &UOld);

    void calcCellADEDown(double timestep, const Foundation &foundation,
                         std::vector<double> &V, std::vector<double> &VOld);

    double calcCellExplicit(double timestep, const Foundation &foundation);

    virtual void calcCellMatrix(Foundation::NumericalScheme scheme, double timestep, const Foundation &foundation,
                        const BoundaryConditions &bcs,
                        std::vector<Eigen::Triplet<double>> &tripletList, Eigen::VectorXd &b,
                        std::vector<double> &a1, std::vector<double> &a2, std::vector<double> &a3,
                        std::vector<double> &b_);

    void calcCellSteadyState(const Foundation &foundation,
                             std::vector<Eigen::Triplet<double>> &tripletList, Eigen::VectorXd &b,
                             std::vector<double> &a1, std::vector<double> &a2, std::vector<double> &a3,
                             std::vector<double> &b_);

    virtual void calcCellADI(int dim, const Foundation &foundation, double timestep,
                             const BoundaryConditions &bcs,
                             double &Am, double &A, double &Ap, double &bVal);

    void Assemble(const Foundation &foundation);

    void setAmatValue(const int i, const int j, const double val, int ndims,
                      std::vector<Eigen::Triplet<double>> &tripletList,
                      std::vector<double> &a1, std::vector<double> &a2, std::vector<double> &a3);
    void setbValue(const int i, const double val, int ndims,
                   Eigen::VectorXd &b, std::vector<double> &b_);
  };


  class ExteriorAirCell : public Cell {
  public:

    ExteriorAirCell(const std::size_t &index, const CellType cellType,
         const std::size_t &i, const std::size_t &j, const std::size_t &k,
         const Foundation &foundation, Surface *surfacePtr, Block *blockPtr,
         Mesher *meshX, Mesher *meshY, Mesher *meshZ);

    void calcCellMatrix(Foundation::NumericalScheme scheme, double timestep, const Foundation &foundation,
                        const BoundaryConditions &bcs,
                        std::vector<Eigen::Triplet<double>> &tripletList, Eigen::VectorXd &b,
                        std::vector<double> &a1, std::vector<double> &a2, std::vector<double> &a3,
                        std::vector<double> &b_) override;
    void calcCellADI(int dim, const Foundation &foundation, double timestep,
                   const BoundaryConditions &bcs,
                   double &Am, double &A, double &Ap, double &bVal) override;
  };


  class InteriorAirCell : public Cell {
  public:

    InteriorAirCell(const std::size_t &index, const CellType cellType,
                    const std::size_t &i, const std::size_t &j, const std::size_t &k,
                    const Foundation &foundation, Surface *surfacePtr, Block *blockPtr,
                    Mesher *meshX, Mesher *meshY, Mesher *meshZ);

    void calcCellMatrix(Foundation::NumericalScheme scheme, double timestep, const Foundation &foundation,
                        const BoundaryConditions &bcs,
                        std::vector<Eigen::Triplet<double>> &tripletList, Eigen::VectorXd &b,
                        std::vector<double> &a1, std::vector<double> &a2, std::vector<double> &a3,
                        std::vector<double> &b_) override;
    void calcCellADI(int dim, const Foundation &foundation, double timestep,
                     const BoundaryConditions &bcs,
                     double &Am, double &A, double &Ap, double &bVal) override;
  };

  class BoundaryCell : public Cell {
  public:

    BoundaryCell(const std::size_t &index, const CellType cellType,
                    const std::size_t &i, const std::size_t &j, const std::size_t &k,
                    const Foundation &foundation, Surface *surfacePtr, Block *blockPtr,
                    Mesher *meshX, Mesher *meshY, Mesher *meshZ);

    void calcCellMatrix(Kiva::Foundation::NumericalScheme scheme, double timestep,
                        const Kiva::Foundation &foundation, const Kiva::BoundaryConditions &bcs,
                        std::vector<Eigen::Triplet<double>> &tripletList, Eigen::VectorXd &b,
                        std::vector<double> &a1, std::vector<double> &a2, std::vector<double> &a3,
                        std::vector<double> &b_) override;
    void calcCellADI(int dim, const Foundation &foundation, double timestep,
                     const BoundaryConditions &bcs,
                     double &Am, double &A, double &Ap, double &bVal) override;
  };

//
//  class ZeroThicknessCell : public Cell {
//
//  };
//
}

#endif
