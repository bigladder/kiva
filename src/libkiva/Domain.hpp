/* Copyright (c) 2012-2017 Big Ladder Software LLC. All rights reserved.
* See the LICENSE file for additional terms and conditions. */

#ifndef Domain_HPP
#define Domain_HPP

#include "Foundation.hpp"
#include "Mesher.hpp"
#include "Functions.hpp"

#include <fstream>
#include <memory>
#include <numeric>

namespace Kiva {

class Cell
{
public:

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
  // organizational properties
  enum CellType
  {
    EXTERIOR_AIR,  // 0
    INTERIOR_AIR,  // 1
    NORMAL,  // 2
    BOUNDARY,  // 3
    ZERO_THICKNESS  // 4
  };
  CellType cellType;

  Block* blockPtr;

  Surface* surfacePtr;
  Mesher *meshXptr, *meshYptr, *meshZptr;
  Cell* i_up_Ptr, *i_down_Ptr, *j_up_Ptr, *j_down_Ptr, *k_up_Ptr, *k_down_Ptr;

  double getKXP();
  double getKXM();
  double getKYP();
  double getKYM();
  double getKZP();
  double getKZM();
  int getNumZeroDims();
  void setZeroThicknessCellProperties(std::vector<Cell*> pointSet);
};

class Domain
{
public:

    // mesh
    Mesher meshX;
    Mesher meshY;
    Mesher meshZ;
    std::size_t nX, nY, nZ;
    std::size_t stepsize_i, stepsize_j, stepsize_k;

    std::vector<Cell> cell;
    std::vector< std::vector<std::size_t> > dest_index_vector;

public:

    Domain();
    Domain(Foundation &foundation);
    void setDomain(Foundation &foundation);
    double getDXP(std::size_t i);
    double getDXM(std::size_t i);
    double getDYP(std::size_t j);
    double getDYM(std::size_t j);
    double getDZP(std::size_t k);
    double getDZM(std::size_t k);
    void set2DZeroThicknessCellProperties(std::size_t index);
    void set3DZeroThicknessCellProperties(std::size_t index);
    void printCellTypes();
    std::tuple<std::size_t, std::size_t, std::size_t> get_coordinates(std::size_t index);
    std::tuple<std::size_t, std::size_t, std::size_t> get_step_size();
    std::vector<std::size_t> get_dest_index(std::size_t i, std::size_t j, std::size_t k);

};

}

#endif
