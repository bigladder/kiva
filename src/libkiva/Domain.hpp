/* Copyright (c) 2012-2017 Big Ladder Software LLC. All rights reserved.
* See the LICENSE file for additional terms and conditions. */

#ifndef Domain_HPP
#define Domain_HPP

#include "Foundation.hpp"
#include "Mesher.hpp"
#include "Functions.hpp"
#include "Cell.hpp"

#include <fstream>
#include <memory>
#include <numeric>

namespace Kiva {


class Domain
{
public:

    // mesh
    Mesher meshX;
    Mesher meshY;
    Mesher meshZ;
    std::size_t nX, nY, nZ;
    std::size_t stepsize_i, stepsize_j, stepsize_k;

    std::vector< std::shared_ptr<Cell> > cell;
    std::vector< std::vector<std::size_t> > dest_index_vector;

public:

    Domain();
    Domain(Foundation &foundation);
    void setDomain(Foundation &foundation);
    int getNumZeroDims(std::size_t i, std::size_t j, std::size_t k);
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
    void addCell(std::shared_ptr<Cell>);
};

}

#endif
