/* Copyright (c) 2012-2017 Big Ladder Software LLC. All rights reserved.
* See the LICENSE file for additional terms and conditions. */

#ifndef Domain_CPP
#define Domain_CPP

#include "Domain.hpp"

namespace Kiva {


Domain::Domain()
{

}

Domain::Domain(Foundation &foundation)
{
  setDomain(foundation);
}

// Having this separate from the constructor allows the correct resizing of
// multidimensional arrays on pre-existing initialized instances.
void Domain::setDomain(Foundation &foundation)
{

  {
    Mesher mX(foundation.xMeshData);
    meshX = mX;

    Mesher mY(foundation.yMeshData);
    meshY = mY;

    Mesher mZ(foundation.zMeshData);
    meshZ = mZ;
  }

  nX = meshX.centers.size();
  nY = meshY.centers.size();
  nZ = meshZ.centers.size();

  std::vector<double> dxp_vector, dxm_vector, dyp_vector, dym_vector, dzp_vector, dzm_vector;
  for (std::size_t i=0; i<nX; i++) {
      dxp_vector.emplace_back(getDXP(i));
      dxm_vector.emplace_back(getDXM(i));
  }
  for (std::size_t j=0; j<nY; j++) {
      dyp_vector.emplace_back(getDYP(j));
      dym_vector.emplace_back(getDYM(j));
  }
  for (std::size_t k=0; k<nZ; k++) {
      dzp_vector.emplace_back(getDZP(k));
      dzm_vector.emplace_back(getDZM(k));
  }

  stepsize[0] = 1;
  stepsize[1] = nX;
  stepsize[2] = nX*nY;
  cell.reserve(nX*nY*nZ);
  dest_index_vector.resize(3, std::vector<std::size_t>(nX*nY*nZ));
  std::vector<std::size_t> temp_di(3);
  std::size_t i, j, k;
  CellType cellType;

  for (std::size_t index = 0; index < nX*nY*nZ; index++)
  {
    std::tie(i, j, k) = get_coordinates(index);
    temp_di = get_dest_index(i, j, k);
    for (std::size_t d=0; d<3; d++) {
      dest_index_vector[d][index] = temp_di[d];
    }

    cellType = CellType::NORMAL;
    Surface *surfacePtr;

    for (auto &surface: foundation.surfaces)
    {
      if (pointOnPoly(Point(meshX.centers[i],meshY.centers[j]), surface.polygon))
      {
        if (isGreaterOrEqual(meshZ.centers[k], surface.zMin)
            &&  isLessOrEqual(meshZ.centers[k], surface.zMax))
        {
          cellType = CellType::BOUNDARY;
          surfacePtr = &surface;
          // Point/Line cells not on the boundary should be
          // zero-thickness cells
          int numZeroDims = getNumZeroDims(i, j, k);

          if (foundation.numberOfDimensions == 3)
          {
            if ((numZeroDims > 1) &&
                i != 0 && i != nX - 1 &&
                j != 0 && j != nY - 1 &&
                k != 0 && k != nZ - 1) {
              cellType = CellType::ZERO_THICKNESS;
            }
          }
          else if (foundation.numberOfDimensions == 2)
          {
            if ((numZeroDims > 1) &&
                i != 0 && i != nX - 1 &&
                k != 0 && k != nZ - 1) {
              cellType = CellType::ZERO_THICKNESS;
            }
          }
          else
          {
            if ((numZeroDims > 1) &&
                k != 0 && k != nZ - 1) {
              cellType = CellType::ZERO_THICKNESS;
            }
          }

          if (cellType == CellType::BOUNDARY) {
            surface.indices.push_back(index);
          }
        }
      }
    }

    // this cell creation needs to be separate from the previous for loop to prevent double-instantiation.
    if (cellType == CellType::ZERO_THICKNESS) {
      addCell(std::make_shared<ZeroThicknessCell>(index, cellType, i, j, k, stepsize, foundation, surfacePtr, nullptr,
                                          &meshX, &meshY, &meshZ));
    } else if (cellType == CellType::BOUNDARY) {
      addCell(std::make_shared<BoundaryCell>(index, cellType, i, j, k, stepsize, foundation, surfacePtr, nullptr,
                                          &meshX, &meshY, &meshZ));
    }

    if (cellType == CellType::NORMAL) {
      for (auto block: foundation.blocks) {
        if (boost::geometry::within(Point(meshX.centers[i], meshY.centers[j]), block.polygon) &&
            isGreaterThan(meshZ.centers[k], block.zMin) &&
            isLessThan(meshZ.centers[k], block.zMax)) {
          if (block.blockType == Block::INTERIOR_AIR) {
            cellType = CellType::INTERIOR_AIR;
            addCell(std::make_shared<InteriorAirCell>(index, cellType, i, j, k, stepsize, foundation, nullptr, &block,
                                                &meshX, &meshY, &meshZ));
          } else if (block.blockType == Block::EXTERIOR_AIR) {
            cellType = CellType::EXTERIOR_AIR;
            addCell(std::make_shared<ExteriorAirCell>(index, cellType, i, j, k, stepsize, foundation, nullptr, &block,
                                                           &meshX, &meshY, &meshZ));
          }
        }
      }
    }

    if (cellType == CellType::NORMAL) {
      // If not surface or block, find interior zero-width cells
      if (foundation.numberOfDimensions == 3) {
        if (isEqual(meshX.deltas[i], 0.0) ||
            isEqual(meshZ.deltas[k], 0.0) ||
            isEqual(meshY.deltas[j], 0.0)) {
          cellType = CellType::ZERO_THICKNESS;
        }
      } else if (foundation.numberOfDimensions == 2) {
        if (isEqual(meshX.deltas[i], 0.0) ||
            isEqual(meshZ.deltas[k], 0.0)) {
          cellType = CellType::ZERO_THICKNESS;
        }
      } else {
        if (isEqual(meshZ.deltas[k], 0.0)) {
          cellType = CellType::ZERO_THICKNESS;
        }
      }

      if (cellType == CellType::ZERO_THICKNESS) {
        addCell(std::make_shared<ZeroThicknessCell>(index, cellType, i, j, k, stepsize, foundation, nullptr, nullptr,
                                            &meshX, &meshY, &meshZ));
      } else {
        addCell(std::make_shared<Cell>(index, cellType, i, j, k, stepsize, foundation, nullptr, nullptr,
                                            &meshX, &meshY, &meshZ));
      }
    }
  }

  // Set effective properties of zero-thickness cells
  // based on other cells
  for (auto this_cell: cell)
  {
    std::size_t index = this_cell->index;
    std::tie(i, j, k) = get_coordinates(index);

    int numZeroDims = getNumZeroDims(i, j, k);

    if (numZeroDims > 0
        && this_cell->cellType != CellType::INTERIOR_AIR
        && this_cell->cellType != CellType::EXTERIOR_AIR)
    {
      if (foundation.numberOfDimensions == 3)
      {
        if (i != 0 && i != nX - 1 &&
          j != 0 && j != nY - 1 &&
          k != 0 && k != nZ - 1)
          set3DZeroThicknessCellProperties(index);
      }
      else if (foundation.numberOfDimensions == 2)
      {
        if (i != 0 && i != nX - 1 && k != 0 && k != nZ - 1)
          set2DZeroThicknessCellProperties(index);
      }
      else
      {
        if (k != 0 && k != nZ - 1)
        {
          if (isEqual(meshZ.deltas[k], 0.0))
          {
            std::vector< std::shared_ptr<Cell> > pointSet =
              {cell[index - stepsize[2]],
               cell[index + stepsize[2]]};

            this_cell->setZeroThicknessCellProperties(pointSet);
          }
        }
      }
    }
  }

  // Calculate matrix coefficients
  for (auto this_cell: cell)
  {
    this_cell->setNeighbors(cell, nX, nY, nZ);

    // PDE Coefficients
    this_cell->setDistances(dxp_vector[this_cell->i], dxm_vector[this_cell->i],
                           dyp_vector[this_cell->j], dym_vector[this_cell->j],
                           dzp_vector[this_cell->k], dzm_vector[this_cell->k]);
    this_cell->setConductivities();
    this_cell->setPDEcoefficients(foundation.numberOfDimensions,
                                  foundation.coordinateSystem == Foundation::CS_CYLINDRICAL);
  }

  for (auto &surface: foundation.surfaces)
  {
    surface.calcTilt();
    surface.area = 0;
    for (auto index: surface.indices)
    {
      surface.area += cell[index]->area;
    }
  }
}

double Domain::getDXP(std::size_t i)
{
  if (i == nX - 1)
  {
    // For boundary cells assume that the cell on the other side of the
    // boundary is the same as the previous cell
    return (meshX.deltas[i] + meshX.deltas[i - 1])/2.0;
  }
  else
  {
    return (meshX.deltas[i] + meshX.deltas[i + 1])/2.0;
  }
}

double Domain::getDXM(std::size_t i)
{
  if (i == 0)
  {
    // For boundary cells assume that the cell on the other side of the
    // boundary is the same as the previous cell
    return (meshX.deltas[i] + meshX.deltas[i + 1])/2.0;
  }
  else
  {
    return (meshX.deltas[i] + meshX.deltas[i - 1])/2.0;
  }
}

double Domain::getDYP(std::size_t j)
{
  if (j == nY - 1)
  {
    // For boundary cells assume that the cell on the other side of the
    // boundary is the same as the previous cell
    return (meshY.deltas[j] + meshY.deltas[j - 1])/2.0;
  }
  else
  {
    return (meshY.deltas[j] + meshY.deltas[j + 1])/2.0;
  }
}

double Domain::getDYM(std::size_t j)
{
  if (j == 0)
  {
    // For boundary cells assume that the cell on the other side of the
    // boundary is the same as the previous cell
    return (meshY.deltas[j] + meshY.deltas[j + 1])/2.0;
  }
  else
  {
    return (meshY.deltas[j] + meshY.deltas[j - 1])/2.0;
  }
}

double Domain::getDZP(std::size_t k)
{
  if (k == nZ - 1)
  {
    // For boundary cells assume that the cell on the other side of the
    // boundary is the same as the previous cell
    return (meshZ.deltas[k] + meshZ.deltas[k - 1])/2.0;
  }
  else
  {
    return (meshZ.deltas[k] + meshZ.deltas[k + 1])/2.0;
  }
}

double Domain::getDZM(std::size_t k)
{
  if (k == 0)
  {
    // For boundary cells assume that the cell on the other side of the
    // boundary is the same as the previous cell
    return (meshZ.deltas[k] + meshZ.deltas[k + 1])/2.0;
  }
  else
  {
    return (meshZ.deltas[k] + meshZ.deltas[k - 1])/2.0;
  }
}

void Domain::set2DZeroThicknessCellProperties(std::size_t index)
{
  if (isEqual(meshX.deltas[cell[index]->i], 0.0) &&
    isEqual(meshZ.deltas[cell[index]->k], 0.0))
  {
    std::vector< std::shared_ptr<Cell> > pointSet =
            {cell[index - stepsize[0] + stepsize[2]],
             cell[index + stepsize[0] + stepsize[2]],
             cell[index - stepsize[0] - stepsize[2]],
             cell[index + stepsize[0] - stepsize[2]]};
    cell[index]->setZeroThicknessCellProperties(pointSet);
  }
  else if (isEqual(meshX.deltas[cell[index]->i], 0.0))
  {
    std::vector< std::shared_ptr<Cell> > pointSet =
            {cell[index - stepsize[0]],
             cell[index + stepsize[0]]};
    cell[index]->setZeroThicknessCellProperties(pointSet);
  }
  else if (isEqual(meshZ.deltas[cell[index]->k], 0.0))
  {
    std::vector< std::shared_ptr<Cell> > pointSet =
            {cell[index - stepsize[2]],
             cell[index + stepsize[2]]};
    cell[index]->setZeroThicknessCellProperties(pointSet);
  }

}

void Domain::set3DZeroThicknessCellProperties(std::size_t index)
{
  if (isEqual(meshX.deltas[cell[index]->i], 0.0) &&
    isEqual(meshY.deltas[cell[index]->j], 0.0) &&
    isEqual(meshZ.deltas[cell[index]->k], 0.0))
  {
    // Use all 8 full volume cells
    std::vector< std::shared_ptr<Cell> > pointSet =
            {cell[index - stepsize[0] - stepsize[1] + stepsize[2]],
             cell[index + stepsize[0] - stepsize[1] + stepsize[2]],
             cell[index - stepsize[0] - stepsize[1] - stepsize[2]],
             cell[index + stepsize[0] - stepsize[1] - stepsize[2]],
             cell[index - stepsize[0] + stepsize[1] + stepsize[2]],
             cell[index + stepsize[0] + stepsize[1] + stepsize[2]],
             cell[index - stepsize[0] + stepsize[1] - stepsize[2]],
             cell[index + stepsize[0] + stepsize[1] - stepsize[2]]};

    cell[index]->setZeroThicknessCellProperties(pointSet);
  }
  else if (isEqual(meshX.deltas[cell[index]->i], 0.0) &&
    isEqual(meshY.deltas[cell[index]->j], 0.0))
  {
    std::vector< std::shared_ptr<Cell> > pointSet =
            {cell[index - stepsize[0] - stepsize[1]],
             cell[index + stepsize[0] - stepsize[1]],
             cell[index - stepsize[0] + stepsize[1]],
             cell[index + stepsize[0] + stepsize[1]]};

    cell[index]->setZeroThicknessCellProperties(pointSet);
  }
  else if (isEqual(meshX.deltas[cell[index]->i], 0.0) &&
    isEqual(meshZ.deltas[cell[index]->k], 0.0))
  {
    std::vector< std::shared_ptr<Cell> > pointSet =
            {cell[index - stepsize[0] + stepsize[2]],
             cell[index + stepsize[0] + stepsize[2]],
             cell[index - stepsize[0] - stepsize[2]],
             cell[index + stepsize[0] - stepsize[2]]};

    cell[index]->setZeroThicknessCellProperties(pointSet);
  }
  else if (isEqual(meshY.deltas[cell[index]->j], 0.0) &&
    isEqual(meshZ.deltas[cell[index]->k], 0.0))
  {
    std::vector< std::shared_ptr<Cell> > pointSet =
            {cell[index - stepsize[1] + stepsize[2]],
             cell[index + stepsize[1] + stepsize[2]],
             cell[index - stepsize[1] - stepsize[2]],
             cell[index + stepsize[1] - stepsize[2]],};

    cell[index]->setZeroThicknessCellProperties(pointSet);
  }
  else if (isEqual(meshX.deltas[cell[index]->i], 0.0))
  {
    std::vector< std::shared_ptr<Cell> > pointSet =
            {cell[index + stepsize[0]],
             cell[index - stepsize[0]]};

    cell[index]->setZeroThicknessCellProperties(pointSet);
  }
  else if (isEqual(meshY.deltas[cell[index]->j], 0.0))
  {
    std::vector< std::shared_ptr<Cell> > pointSet =
            {cell[index + stepsize[1]],
             cell[index - stepsize[1]]};

    cell[index]->setZeroThicknessCellProperties(pointSet);
  }
  else if (isEqual(meshZ.deltas[cell[index]->k], 0.0))
  {
    std::vector< std::shared_ptr<Cell> > pointSet =
            {cell[index + stepsize[2]],
             cell[index - stepsize[2]]};

    cell[index]->setZeroThicknessCellProperties(pointSet);
  }

}

int Domain::getNumZeroDims(std::size_t i, std::size_t j, std::size_t k)
{
  int numZeroDims = 0;
  if (isEqual(meshX.deltas[i], 0.0))
    numZeroDims += 1;
  if (isEqual(meshY.deltas[j], 0.0))
    numZeroDims += 1;
  if (isEqual(meshZ.deltas[k], 0.0))
    numZeroDims += 1;

  return numZeroDims;
}

void Domain::printCellTypes()
{
  // TODO: Make the ability to output a specific slice in i, j, or k
  std::ofstream output;
  output.open("Cells.csv");

  for (std::size_t i = 0; i < nX; i++)
  {

    output << ", " << i;

  }

  output << "\n";

  for (std::size_t k = nZ - 1; /* k >= 0 && */ k < nZ; k--)
  {

    output << k;

    for (std::size_t i = 0; i < nX; i++)
    {

      output << ", " << cell[i + (nY/2)*stepsize[1] + k*stepsize[2]]->cellType;

    }

    output << "\n";
  }
  output.close();

}

std::tuple<std::size_t, std::size_t, std::size_t> Domain::get_coordinates(std::size_t index){
    size_t i, j, k;
    i = index % nX;
    j = ((index - i) % nY) / nX;
    k = (index - i - nX*j) / (nX*nY);
    return std::make_tuple(i, j, k);
}

std::vector<std::size_t> Domain::get_dest_index(std::size_t i, std::size_t j, std::size_t k)
{
  std::vector<std::size_t> dest_index;
  dest_index.emplace_back(i + nX*j + nX*nY*k);
  dest_index.emplace_back(j + nY*i + nY*nX*k);
  dest_index.emplace_back(k + nZ*i + nZ*nX*j);
  return dest_index;
}

void Domain::addCell(std::shared_ptr<Cell> this_cell)
{
  cell.emplace_back(std::move(this_cell));
}

}

#endif
