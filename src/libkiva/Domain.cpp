/* Copyright (c) 2012-2017 Big Ladder Software LLC. All rights reserved.
* See the LICENSE file for additional terms and conditions. */

#ifndef Domain_CPP
#define Domain_CPP

#include "Domain.hpp"

namespace Kiva {

static const double PI = 4.0*atan(1.0);

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

  std::tie(stepsize_i, stepsize_j, stepsize_k) = get_step_size();

  cell.resize(nX*nY*nZ);
  dest_index_vector.resize(nX*nY*nZ, std::vector<std::size_t>(3));

  for (std::size_t index = 0; index < nX*nY*nZ; index++)
  {
        std::size_t i, j, k;
        std::tie(i, j, k) = get_coordinates(index);

        Cell* this_cell = &cell[index];
        this_cell->index = index;
        dest_index_vector[index] = get_dest_index(i, j, k);
        // Set Cell Properties
        this_cell->density = foundation.soil.density;
        this_cell->specificHeat = foundation.soil.specificHeat;
        this_cell->conductivity = foundation.soil.conductivity;
        this_cell->heatGain = 0.0;
        this_cell->i = i;
        this_cell->j = j;
        this_cell->k = k;

        this_cell->i_up_Ptr = &cell[index+stepsize_i];
        this_cell->i_down_Ptr = &cell[index-stepsize_i];
        this_cell->j_up_Ptr = &cell[index+stepsize_j];
        this_cell->j_down_Ptr = &cell[index-stepsize_j];
        this_cell->k_up_Ptr = &cell[index+stepsize_k];
        this_cell->k_down_Ptr = &cell[index-stepsize_k];
        this_cell->meshXptr = &meshX;
        this_cell->meshYptr = &meshY;
        this_cell->meshZptr = &meshZ;

        // Default to normal cells
        this_cell->cellType = Cell::NORMAL;

        // Next set interior zero-width cells
        if (foundation.numberOfDimensions == 3)
        {
          if (isEqual(meshX.deltas[i], 0.0) ||
            isEqual(meshZ.deltas[k], 0.0) ||
            isEqual(meshY.deltas[j], 0.0))
          {
            this_cell->cellType = Cell::ZERO_THICKNESS;
          }
        }
        else if (foundation.numberOfDimensions == 2)
        {
          if (isEqual(meshX.deltas[i], 0.0) ||
            isEqual(meshZ.deltas[k], 0.0))
          {
            this_cell->cellType = Cell::ZERO_THICKNESS;
          }
        }
        else
        {
          if (isEqual(meshZ.deltas[k], 0.0))
          {
            this_cell->cellType = Cell::ZERO_THICKNESS;
          }
        }

        for (std::size_t b = 0; b < foundation.blocks.size(); b++)
        {
          if (boost::geometry::within(Point(meshX.centers[i],meshY.centers[j]), foundation.blocks[b].polygon) &&
            isGreaterThan(meshZ.centers[k], foundation.blocks[b].zMin) &&
            isLessThan(meshZ.centers[k], foundation.blocks[b].zMax))
          {
            this_cell->density = foundation.blocks[b].material.density;
            this_cell->specificHeat = foundation.blocks[b].material.specificHeat;
            this_cell->conductivity = foundation.blocks[b].material.conductivity;

            this_cell->blockPtr = &foundation.blocks[b];

            if (foundation.blocks[b].blockType == Block::INTERIOR_AIR)
            {
              this_cell->cellType = Cell::INTERIOR_AIR;
            }
            else if (foundation.blocks[b].blockType == Block::EXTERIOR_AIR)
            {
              this_cell->cellType = Cell::EXTERIOR_AIR;
            } else {
              this_cell->cellType = Cell::NORMAL;
            }
          }
        }

        for (std::size_t s = 0; s < foundation.surfaces.size(); s++)
        {
          if (pointOnPoly(Point(meshX.centers[i],meshY.centers[j]), foundation.surfaces[s].polygon))
          {
            if (isGreaterOrEqual(meshZ.centers[k], foundation.surfaces[s].zMin)
            &&  isLessOrEqual(meshZ.centers[k], foundation.surfaces[s].zMax))
            {
              this_cell->cellType = Cell::BOUNDARY;

              this_cell->surfacePtr = &foundation.surfaces[s];

              // Point/Line cells not on the boundary should be
              // zero-thickness cells
              int numZeroDims = this_cell->getNumZeroDims();

              if (foundation.numberOfDimensions == 3)
              {
                if ((numZeroDims > 1) &&
                  i != 0 && i != nX - 1 &&
                  j != 0 && j != nY - 1 &&
                  k != 0 && k != nZ - 1)
                  this_cell->cellType = Cell::ZERO_THICKNESS;
              }
              else if (foundation.numberOfDimensions == 2)
              {
                if ((numZeroDims > 1) &&
                  i != 0 && i != nX - 1 &&
                  k != 0 && k != nZ - 1)
                  this_cell->cellType = Cell::ZERO_THICKNESS;
              }
              else
              {
                if ((numZeroDims > 1) &&
                  k != 0 && k != nZ - 1)
                  this_cell->cellType = Cell::ZERO_THICKNESS;
              }

              if (this_cell->cellType == Cell::BOUNDARY)
              {
                foundation.surfaces[s].indices.push_back(index);
              }

            }
          }
        }

        // Set cell volume
        this_cell->volume = meshX.deltas[i]*meshY.deltas[j]*meshZ.deltas[k];

        // for boundary cells, set cell area
        if (this_cell->cellType == Cell::BOUNDARY)
        {
          if (foundation.numberOfDimensions == 2 &&
              foundation.coordinateSystem == Foundation::CS_CYLINDRICAL)
          {
            if (this_cell->surfacePtr->orientation == Surface::X_POS ||
              this_cell->surfacePtr->orientation == Surface::X_NEG)
            {
              this_cell->area = 2.0*PI*meshX.centers[i]*meshZ.deltas[k];
            }
            else // if (surface.orientation == Surface::Z_POS ||
               // surface.orientation == Surface::Z_NEG)
            {
              this_cell->area = PI*(meshX.dividers[i+1]*meshX.dividers[i+1] -
            		  meshX.dividers[i]*meshX.dividers[i] );
            }
          }
          else if (foundation.numberOfDimensions == 2 &&
                   foundation.coordinateSystem == Foundation::CS_CARTESIAN)
          {
            if (this_cell->surfacePtr->orientation == Surface::X_POS ||
              this_cell->surfacePtr->orientation == Surface::X_NEG)
            {
              this_cell->area = 2.0*meshZ.deltas[k]*foundation.linearAreaMultiplier;
            }
            else // if (surface.orientation == Surface::Z_POS ||
               // surface.orientation == Surface::Z_NEG)
            {
              this_cell->area = 2.0*meshX.deltas[i]*foundation.linearAreaMultiplier;
            }
          }
          else if (foundation.numberOfDimensions == 3)
          {
            if (this_cell->surfacePtr->orientation == Surface::X_POS ||
              this_cell->surfacePtr->orientation == Surface::X_NEG)
            {
              this_cell->area = meshY.deltas[j]*meshZ.deltas[k];
            }
            else if (this_cell->surfacePtr->orientation == Surface::Y_POS ||
                 this_cell->surfacePtr->orientation == Surface::Y_NEG)
            {
              this_cell->area = meshX.deltas[i]*meshZ.deltas[k];
            }
            else // if (surface.orientation == Surface::Z_POS ||
               // surface.orientation == Surface::Z_NEG)
            {
              this_cell->area = meshX.deltas[i]*meshY.deltas[j];
            }

            if (foundation.useSymmetry)
            {
              if (foundation.isXSymm)
                this_cell->area = 2*this_cell->area;

              if (foundation.isYSymm)
                this_cell->area = 2*this_cell->area;
            }
          }
          else
          {
            this_cell->area = 1.0;
          }
        }
  }

  // Set effective properties of zero-thickness cells
  // based on other cells
  for (std::size_t index=0; index<nX*nY*nZ; index++)
  {
        Cell* this_cell = &cell[index];
        std::size_t i = this_cell->i;
        std::size_t j = this_cell->j;
        std::size_t k = this_cell->k;

        int numZeroDims = this_cell->getNumZeroDims();

        if (numZeroDims > 0
            && this_cell->cellType != Cell::INTERIOR_AIR
            && this_cell->cellType != Cell::EXTERIOR_AIR)
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
                std::vector<Cell*> pointSet =
                  {&cell[index - stepsize_k],
                   &cell[index + stepsize_k]};

                this_cell->setZeroThicknessCellProperties(pointSet);
              }
            }
          }
        }
  }

  // Calculate matrix coefficients
  for (std::size_t index = 0; index < nX*nY*nZ; index++)
  {
        Cell* this_cell = &cell[index];
        // PDE Coefficients

        this_cell->dxp = dxp_vector[this_cell->i];
        this_cell->dxm = dxm_vector[this_cell->i];
        this_cell->dyp = dyp_vector[this_cell->j];
        this_cell->dym = dym_vector[this_cell->j];
        this_cell->dzp = dzp_vector[this_cell->k];
        this_cell->dzm = dzm_vector[this_cell->k];
        this_cell->kxp = this_cell->getKXP();
        this_cell->kxm = this_cell->getKXM();
        this_cell->kyp = this_cell->getKYP();
        this_cell->kym = this_cell->getKYM();
        this_cell->kzp = this_cell->getKZP();
        this_cell->kzm = this_cell->getKZM();

        if (foundation.numberOfDimensions > 1) {
          // Radial X terms
          if (foundation.coordinateSystem == Foundation::CS_CYLINDRICAL)
          {
            this_cell->cxp_c = (this_cell->dxm*this_cell->kxp)/
                ((this_cell->dxm + this_cell->dxp)*this_cell->dxp);
            this_cell->cxm_c = (this_cell->dxp*this_cell->kxm)/
                ((this_cell->dxm + this_cell->dxp)*this_cell->dxm);
          }
          else
          {
            this_cell->cxp_c = 0.0;
            this_cell->cxm_c = 0.0;
          }

          // Cartesian X terms
          this_cell->cxp = (2*this_cell->kxp)/
              ((this_cell->dxm + this_cell->dxp)*this_cell->dxp);
          this_cell->cxm = -1*(2*this_cell->kxm)/
              ((this_cell->dxm + this_cell->dxp)*this_cell->dxm);
        }

        // Cartesian Z terms
        this_cell->czp = (2*this_cell->kzp)/
            ((this_cell->dzm + this_cell->dzp)*this_cell->dzp);
        this_cell->czm = -1*(2*this_cell->kzm)/
            ((this_cell->dzm + this_cell->dzp)*this_cell->dzm);

        // Cartesian Y terms
        if (foundation.numberOfDimensions == 3)
        {
          this_cell->cyp = (2*this_cell->kyp)/
              ((this_cell->dym + this_cell->dyp)*this_cell->dyp);
          this_cell->cym = -1*(2*this_cell->kym)/
              ((this_cell->dym + this_cell->dyp)*this_cell->dym);
        }
        else
        {
          this_cell->cyp = 0.0;
          this_cell->cym = 0.0;
        }
  }

  for (std::size_t s = 0; s < foundation.surfaces.size(); s++)
  {
    foundation.surfaces[s].area = 0;
    for (auto index: foundation.surfaces[s].indices)
    {
      foundation.surfaces[s].area += cell[index].area;
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

double Cell::getKXP()
{
  if (i == meshXptr->centers.size() - 1)
  {
    // For boundary cells assume that the cell on the other side of the
    // boundary is the same as the current cell
    return conductivity;
  }
  else
  {
    return 1/(meshXptr->deltas[i]/(2*dxp*conductivity) +
        meshXptr->deltas[i + 1]/(2*dxp*i_up_Ptr->conductivity));
  }
}

double Cell::getKXM()
{
  if (i == 0)
  {
    // For boundary cells assume that the cell on the other side of the
    // boundary is the same as the current cell
    return conductivity;
  }
  else
  {
    return 1/(meshXptr->deltas[i]/(2*dxm*conductivity) +
        meshXptr->deltas[i - 1]/(2*dxm*i_down_Ptr->conductivity));
  }
}

double Cell::getKYP()
{
  if (j == meshYptr->centers.size() - 1)
  {
    // For boundary cells assume that the cell on the other side of the
    // boundary is the same as the current cell
    return conductivity;
  }
  else
  {
    return 1/(meshYptr->deltas[j]/(2*dyp*conductivity) +
        meshYptr->deltas[j + 1]/(2*dyp*j_up_Ptr->conductivity));
  }
}

double Cell::getKYM()
{
  if (j == 0)
  {
    // For boundary cells assume that the cell on the other side of the
    // boundary is the same as the current cell
    return conductivity;
  }
  else
  {
    return 1/(meshYptr->deltas[j]/(2*dym*conductivity) +
        meshYptr->deltas[j - 1]/(2*dym*j_down_Ptr->conductivity));
  }
}

double Cell::getKZP()
{
  if (k == meshZptr->centers.size() - 1)
  {
    // For boundary cells assume that the cell on the other side of the
    // boundary is the same as the current cell
    return conductivity;
  }
  else
  {
    return 1/(meshZptr->deltas[k]/(2*dzp*conductivity) +
        meshZptr->deltas[k + 1]/(2*dzp*k_up_Ptr->conductivity));
  }
}

double Cell::getKZM()
{
  if (k == 0)
  {
    // For boundary cells assume that the cell on the other side of the
    // boundary is the same as the current cell
    return conductivity;
  }
  else
  {
    return 1/(meshZptr->deltas[k]/(2*dzm*conductivity) +
        meshZptr->deltas[k - 1]/(2*dzm*k_down_Ptr->conductivity));
  }
}

void Domain::set2DZeroThicknessCellProperties(std::size_t index)
{
  if (isEqual(meshX.deltas[cell[index].i], 0.0) &&
    isEqual(meshZ.deltas[cell[index].k], 0.0))
  {
    std::vector<Cell*> pointSet =
            {&cell[index - stepsize_i + stepsize_k],
             &cell[index + stepsize_i + stepsize_k],
             &cell[index - stepsize_i - stepsize_k],
             &cell[index + stepsize_i - stepsize_k]};
    cell[index].setZeroThicknessCellProperties(pointSet);
  }
  else if (isEqual(meshX.deltas[cell[index].i], 0.0))
  {
    std::vector<Cell*> pointSet =
            {&cell[index - stepsize_i],
             &cell[index + stepsize_i]};
    cell[index].setZeroThicknessCellProperties(pointSet);
  }
  else if (isEqual(meshZ.deltas[cell[index].k], 0.0))
  {
    std::vector<Cell*> pointSet =
            {&cell[index - stepsize_k],
             &cell[index + stepsize_k]};
    cell[index].setZeroThicknessCellProperties(pointSet);
  }

}

void Domain::set3DZeroThicknessCellProperties(std::size_t index)
{
  if (isEqual(meshX.deltas[cell[index].i], 0.0) &&
    isEqual(meshY.deltas[cell[index].j], 0.0) &&
    isEqual(meshZ.deltas[cell[index].k], 0.0))
  {
    // Use all 8 full volume cells
    std::vector<Cell*> pointSet =
            {&cell[index - stepsize_i - stepsize_j + stepsize_k],
             &cell[index + stepsize_i - stepsize_j + stepsize_k],
             &cell[index - stepsize_i - stepsize_j - stepsize_k],
             &cell[index + stepsize_i - stepsize_j - stepsize_k],
             &cell[index - stepsize_i + stepsize_j + stepsize_k],
             &cell[index + stepsize_i + stepsize_j + stepsize_k],
             &cell[index - stepsize_i + stepsize_j - stepsize_k],
             &cell[index + stepsize_i + stepsize_j - stepsize_k]};

    cell[index].setZeroThicknessCellProperties(pointSet);
  }
  else if (isEqual(meshX.deltas[cell[index].i], 0.0) &&
    isEqual(meshY.deltas[cell[index].j], 0.0))
  {
    std::vector<Cell*> pointSet =
            {&cell[index - stepsize_i - stepsize_j],
             &cell[index + stepsize_i - stepsize_j],
             &cell[index - stepsize_i + stepsize_j],
             &cell[index + stepsize_i + stepsize_j]};

    cell[index].setZeroThicknessCellProperties(pointSet);
  }
  else if (isEqual(meshX.deltas[cell[index].i], 0.0) &&
    isEqual(meshZ.deltas[cell[index].k], 0.0))
  {
    std::vector<Cell*> pointSet =
            {&cell[index - stepsize_i + stepsize_k],
             &cell[index + stepsize_i + stepsize_k],
             &cell[index - stepsize_i - stepsize_k],
             &cell[index + stepsize_i - stepsize_k]};

    cell[index].setZeroThicknessCellProperties(pointSet);
  }
  else if (isEqual(meshY.deltas[cell[index].j], 0.0) &&
    isEqual(meshZ.deltas[cell[index].k], 0.0))
  {
    std::vector<Cell*> pointSet =
            {&cell[index - stepsize_j + stepsize_k],
             &cell[index + stepsize_j + stepsize_k],
             &cell[index - stepsize_j - stepsize_k],
             &cell[index + stepsize_j - stepsize_k],};

    cell[index].setZeroThicknessCellProperties(pointSet);
  }
  else if (isEqual(meshX.deltas[cell[index].i], 0.0))
  {
    std::vector<Cell*> pointSet =
            {&cell[index + stepsize_i],
             &cell[index - stepsize_i]};

    cell[index].setZeroThicknessCellProperties(pointSet);
  }
  else if (isEqual(meshY.deltas[cell[index].j], 0.0))
  {
    std::vector<Cell*> pointSet =
            {&cell[index + stepsize_j],
             &cell[index - stepsize_j]};

    cell[index].setZeroThicknessCellProperties(pointSet);
  }
  else if (isEqual(meshZ.deltas[cell[index].k], 0.0))
  {
    std::vector<Cell*> pointSet =
            {&cell[index + stepsize_k],
             &cell[index - stepsize_k]};

    cell[index].setZeroThicknessCellProperties(pointSet);
  }

}

void Cell::setZeroThicknessCellProperties(std::vector<Cell*> pointSet)
{
  std::vector<double> volumes;
  std::vector<double> densities;
  std::vector<double> specificHeats;
  std::vector<double> conductivities;

  std::vector<double> masses;
  std::vector<double> capacities;
  std::vector<double> weightedConductivity;

  for (auto p_cell : pointSet)
  {
    // Do not add air cell properties into the weighted average
    if (p_cell->cellType != Cell::INTERIOR_AIR &&
      p_cell->cellType != Cell::EXTERIOR_AIR)
    {
    double vol = p_cell->volume;
    double rho = p_cell->density;
    double cp = p_cell->specificHeat;
    double kth = p_cell->conductivity;

    volumes.push_back(vol);
    masses.push_back(vol*rho);
    capacities.push_back(vol*rho*cp);
    weightedConductivity.push_back(vol*kth);
    }
  }

  // if the neighboring cells are all air cells set properties to air properties
  if (volumes.size() == 0)
  {
    volumes.push_back(1.0);
    masses.push_back(1.275);
    capacities.push_back(1.275*1007);
    weightedConductivity.push_back(0.02587);

  }

  double totalVolume = std::accumulate(volumes.begin(), volumes.end(), 0.0);

  density = std::accumulate(masses.begin(), masses.end(), 0.0) /
      totalVolume;

  specificHeat = std::accumulate(capacities.begin(), capacities.end(), 0.0) /
      (totalVolume*density);

  conductivity = std::accumulate(weightedConductivity.begin(), weightedConductivity.end(), 0.0) /
      totalVolume;
}

int Cell::getNumZeroDims()
{
  int numZeroDims = 0;
  if (isEqual(meshXptr->deltas[i], 0.0))
    numZeroDims += 1;
  if (isEqual(meshYptr->deltas[j], 0.0))
    numZeroDims += 1;
  if (isEqual(meshZptr->deltas[k], 0.0))
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

      output << ", " << cell[i + (nY/2)*stepsize_j + k*stepsize_k].cellType;

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

std::tuple<std::size_t, std::size_t, std::size_t> Domain::get_step_size()
{
  size_t i_step, j_step, k_step;
  i_step = 1;
  j_step = nX;
  k_step = nX*nY;
  return std::make_tuple(i_step, j_step, k_step);
}

std::vector<std::size_t> Domain::get_dest_index(std::size_t i, std::size_t j, std::size_t k)
{
  std::vector<std::size_t> dest_index;
  dest_index.emplace_back(i + nX*j + nX*nY*k);
  dest_index.emplace_back(j + nY*i + nY*nX*k);
  dest_index.emplace_back(k + nZ*i + nZ*nX*j);
  return dest_index;
}

}

#endif
