/* Copyright (c) 2012-2017 Big Ladder Software LLC. All rights reserved.
* See the LICENSE file for additional terms and conditions. */

#ifndef Cell_CPP
#define Cell_CPP

#include "Cell.hpp"

namespace Kiva {

static const double PI = 4.0*atan(1.0);


Cell::Cell(const std::size_t &index, const CellType cellType,
           const std::size_t &i, const std::size_t &j, const std::size_t &k,
           const Foundation &foundation, Surface *surfacePtr, Block *blockPtr,
           Mesher *meshXptr, Mesher *meshYptr, Mesher *meshZptr):
        index(index),
        i(i), j(j), k(k),
        cellType(cellType),
        blockPtr(blockPtr),
        surfacePtr(surfacePtr),
        meshXptr(meshXptr),
        meshYptr(meshYptr),
        meshZptr(meshZptr)
{
  if (blockPtr) {
    density = blockPtr->material.density;
    specificHeat = blockPtr->material.specificHeat;
    conductivity = blockPtr->material.conductivity;
  } else {
    density = foundation.soil.density;
    specificHeat = foundation.soil.specificHeat;
    conductivity = foundation.soil.conductivity;
  }
  heatGain = 0.0;
  volume = meshXptr->deltas[i]*meshYptr->deltas[j]*meshZptr->deltas[k];

  // for boundary cells, set cell area
  if (cellType == CellType::BOUNDARY)
  {
    if (foundation.numberOfDimensions == 2 &&
        foundation.coordinateSystem == Foundation::CS_CYLINDRICAL)
    {
      if (surfacePtr->orientation == Surface::X_POS ||
          surfacePtr->orientation == Surface::X_NEG)
      {
        area = 2.0*PI*meshXptr->centers[i]*meshZptr->deltas[k];
      }
      else // if (surface.orientation == Surface::Z_POS ||
        // surface.orientation == Surface::Z_NEG)
      {
        area = PI*(meshXptr->dividers[i+1]*meshXptr->dividers[i+1] -
                   meshXptr->dividers[i]*meshXptr->dividers[i] );
      }
    }
    else if (foundation.numberOfDimensions == 2 &&
             foundation.coordinateSystem == Foundation::CS_CARTESIAN)
    {
      if (surfacePtr->orientation == Surface::X_POS ||
          surfacePtr->orientation == Surface::X_NEG)
      {
        area = 2.0*meshZptr->deltas[k]*foundation.linearAreaMultiplier;
      }
      else // if (surface.orientation == Surface::Z_POS ||
        // surface.orientation == Surface::Z_NEG)
      {
        area = 2.0*meshXptr->deltas[i]*foundation.linearAreaMultiplier;
      }
    }
    else if (foundation.numberOfDimensions == 3)
    {
      if (surfacePtr->orientation == Surface::X_POS ||
          surfacePtr->orientation == Surface::X_NEG)
      {
        area = meshYptr->deltas[j]*meshZptr->deltas[k];
      }
      else if (surfacePtr->orientation == Surface::Y_POS ||
               surfacePtr->orientation == Surface::Y_NEG)
      {
        area = meshXptr->deltas[i]*meshZptr->deltas[k];
      }
      else // if (surface.orientation == Surface::Z_POS ||
        // surface.orientation == Surface::Z_NEG)
      {
        area = meshXptr->deltas[i]*meshYptr->deltas[j];
      }

      if (foundation.useSymmetry)
      {
        if (foundation.isXSymm)
          area = 2*area;

        if (foundation.isYSymm)
          area = 2*area;
      }
    }
    else
    {
      area = 1.0;
    }
  }

  if (cellType == CellType::ZERO_THICKNESS & foundation.numberOfDimensions == 2) {
    if (i != 0)
      r = meshXptr->centers[i];
  }
}



void Cell::setDistances(double &dxp_in, double &dxm_in, double &dyp_in, double &dym_in,
                        double &dzp_in, double &dzm_in) {
  dxp = dxp_in;
  dxm = dxm_in;
  dyp = dyp_in;
  dym = dym_in;
  dzp = dzp_in;
  dzm = dzm_in;
}

void Cell::setConductivities() {
  getKXP();
  getKXM();
  getKYP();
  getKYM();
  getKZP();
  getKZM();
}

void Cell::getKXP()
{
  if (i == meshXptr->centers.size() - 1)
  {
    // For boundary cells assume that the cell on the other side of the
    // boundary is the same as the current cell
    kxp = conductivity;
  }
  else
  {
    kxp = 1/(meshXptr->deltas[i]/(2*dxp*conductivity) +
             meshXptr->deltas[i + 1]/(2*dxp*i_up_Ptr->conductivity));
  }
}

void Cell::getKXM()
{
  if (i == 0)
  {
    // For boundary cells assume that the cell on the other side of the
    // boundary is the same as the current cell
    kxm = conductivity;
  }
  else
  {
    kxm = 1/(meshXptr->deltas[i]/(2*dxm*conductivity) +
             meshXptr->deltas[i - 1]/(2*dxm*i_down_Ptr->conductivity));
  }
}

void Cell::getKYP()
{
  if (j == meshYptr->centers.size() - 1)
  {
    // For boundary cells assume that the cell on the other side of the
    // boundary is the same as the current cell
    kyp = conductivity;
  }
  else
  {
    kyp = 1/(meshYptr->deltas[j]/(2*dyp*conductivity) +
             meshYptr->deltas[j + 1]/(2*dyp*j_up_Ptr->conductivity));
  }
}

void Cell::getKYM()
{
  if (j == 0)
  {
    // For boundary cells assume that the cell on the other side of the
    // boundary is the same as the current cell
    kym = conductivity;
  }
  else
  {
    kym = 1/(meshYptr->deltas[j]/(2*dym*conductivity) +
             meshYptr->deltas[j - 1]/(2*dym*j_down_Ptr->conductivity));
  }
}

void Cell::getKZP()
{
  if (k == meshZptr->centers.size() - 1)
  {
    // For boundary cells assume that the cell on the other side of the
    // boundary is the same as the current cell
    kzp = conductivity;
  }
  else
  {
    kzp = 1/(meshZptr->deltas[k]/(2*dzp*conductivity) +
             meshZptr->deltas[k + 1]/(2*dzp*k_up_Ptr->conductivity));
  }
}

void Cell::getKZM()
{
  if (k == 0)
  {
    // For boundary cells assume that the cell on the other side of the
    // boundary is the same as the current cell
    kzm = conductivity;
  }
  else
  {
    kzm = 1/(meshZptr->deltas[k]/(2*dzm*conductivity) +
             meshZptr->deltas[k - 1]/(2*dzm*k_down_Ptr->conductivity));
  }
}

void Cell::setWhatever(int ndims, bool cylindrical) {
  if (ndims > 1) {
    // Radial X terms
    if (cylindrical)
    {
      cxp_c = (dxm*kxp) / ((dxm + dxp)*dxp);
      cxm_c = (dxp*kxm) / ((dxm + dxp)*dxm);
    }
    else
    {
      cxp_c = 0.0;
      cxm_c = 0.0;
    }

    // Cartesian X terms
    cxp = (2*kxp) /  ((dxm + dxp)*dxp);
    cxm = -1*(2*kxm) / ((dxm + dxp)*dxm);
  }

  // Cartesian Z terms
  czp = (2*kzp) / ((dzm + dzp)*dzp);
  czm = -1*(2*kzm) / ((dzm + dzp)*dzm);

  // Cartesian Y terms
  if (ndims == 3)
  {
    cyp = (2*kyp) / ((dym + dyp)*dyp);
    cym = -1*(2*kym) / ((dym + dyp)*dym);
  }
  else
  {
    cyp = 0.0;
    cym = 0.0;
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
    if (p_cell->cellType != CellType::INTERIOR_AIR &&
        p_cell->cellType != CellType::EXTERIOR_AIR)
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


}

#endif
