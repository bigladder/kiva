/* Copyright (c) 2012-2017 Big Ladder Software LLC. All rights reserved.
* See the LICENSE file for additional terms and conditions. */

#ifndef Cell_CPP
#define Cell_CPP

#include "Cell.hpp"

namespace Kiva {

static const double PI = 4.0*atan(1.0);

Cell::Cell(const std::size_t &index, const CellType cellType,
           const std::size_t &i, const std::size_t &j, const std::size_t &k,
           std::size_t *stepsize,
           const Foundation &foundation, Surface *surfacePtr, Block *blockPtr,
           Mesher *meshPtr):
        coords{i, j, k},
        index(index),
        stepsize(stepsize),
        cellType(cellType),
        blockPtr(blockPtr),
        surfacePtr(surfacePtr),
        meshPtr(meshPtr)
{
  Assemble(foundation);
}

void Cell::Assemble(const Foundation &foundation) {
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
  volume = meshPtr[0].deltas[coords[0]] * meshPtr[1].deltas[coords[1]] * meshPtr[2].deltas[coords[2]];

  if (foundation.numberOfDimensions == 2) {
      r = meshPtr[0].centers[coords[0]];
  }
}

void Cell::setDistances(const double &dxp_in, const double &dxm_in, const double &dyp_in, const double &dym_in,
                        const double &dzp_in, const double &dzm_in) {
  dist[0][1] = dxp_in;
  dist[0][0] = dxm_in;
  dist[1][1] = dyp_in;
  dist[1][0] = dym_in;
  dist[2][1] = dzp_in;
  dist[2][0] = dzm_in;
}

  void Cell::setConductivities(const std::vector< std::shared_ptr<Cell> > &cell_v)
  {
    for (std::size_t dim=0; dim<3; ++dim) {
      //  dir == 0
      if (coords[dim] == 0) {
        // For boundary cells assume that the cell on the other side of the
        // boundary is the same as the current cell
        kcoeff[dim][0] = conductivity;
      } else {
        kcoeff[dim][0] = 1 / (meshPtr[dim].deltas[coords[dim]] / (2 * dist[dim][0] * conductivity) +
                              meshPtr[dim].deltas[coords[dim] - 1] /
                              (2 * dist[dim][0] * cell_v[index - stepsize[dim]]->conductivity));
      }

      //  dir == 1
      if (coords[dim] == meshPtr[dim].centers.size() - 1) {
        kcoeff[dim][1] = conductivity;
      } else {
        kcoeff[dim][1] = 1 / (meshPtr[dim].deltas[coords[dim]] / (2 * dist[dim][1] * conductivity) +
                              meshPtr[dim].deltas[coords[dim] + 1] /
                              (2 * dist[dim][1] * cell_v[index + stepsize[dim]]->conductivity));
      }
    }
  }

void Cell::setPDEcoefficients(int ndims, bool cylindrical) {
  if (ndims > 1) {
    // Radial X terms
    if (cylindrical)
    {
      pde_c[1] = (dist[0][0]*kcoeff[0][1]) / ((dist[0][0] + dist[0][1])*dist[0][1]);
      pde_c[0] = (dist[0][1]*kcoeff[0][0]) / ((dist[0][0] + dist[0][1])*dist[0][0]);
    }
    else
    {
      pde_c[1] = 0.0;
      pde_c[0] = 0.0;
    }

    // Cartesian X terms
    pde[0][1] = onePDEcoefficient(0, 1);
    pde[0][0] = onePDEcoefficient(0, 0);
  }

  // Cartesian Z terms
  pde[2][1] = onePDEcoefficient(2, 1);
  pde[2][0] = onePDEcoefficient(2, 0);

  // Cartesian Y terms
  if (ndims == 3)
  {
    pde[1][1] = onePDEcoefficient(1, 1);
    pde[1][0] = onePDEcoefficient(1, 0);
  }
  else
  {
    pde[1][1] = 0.0;
    pde[1][0] = 0.0;
  }
}

double Cell::onePDEcoefficient(std::size_t dim, std::size_t dir) {
  int sign = dir == 0? -1 : 1;
  double c = sign * (2*kcoeff[dim][dir]) /  ((dist[dim][0] + dist[dim][1])*dist[dim][dir]);
  return c;
}

void Cell::setZeroThicknessCellProperties(std::vector< std::shared_ptr<Cell> > pointSet)
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



void Cell::calcCellADEUp(double timestep, const Foundation &foundation, const BoundaryConditions &/*bcs*/,
                         double &U)
{
  double theta = timestep/
                 (density*specificHeat);

  double CXP = pde[0][1]*theta;
  double CXM = pde[0][0]*theta;
  double CZP = pde[2][1]*theta;
  double CZM = pde[2][0]*theta;
  double CYP = pde[1][1]*theta;
  double CYM = pde[1][0]*theta;
  double Q = heatGain*theta;

  if (foundation.numberOfDimensions == 3)
    U = (*told_ptr*(1.0 - CXP - CZP - CYP)
            - *(&U - stepsize[0])*CXM
            + *(told_ptr + stepsize[0])*CXP
            - *(&U - stepsize[2])*CZM
            + *(told_ptr + stepsize[2])*CZP
            - *(&U - stepsize[1])*CYM
            + *(told_ptr + stepsize[1])*CYP
            + Q) /
        (1.0 - CXM - CZM - CYM);
  else if (foundation.numberOfDimensions == 2)
  {
    double CXPC = 0;
    double CXMC = 0;

    if (coords[0] != 0)
    {
      CXPC = pde_c[1]*theta/r;
      CXMC = pde_c[0]*theta/r;
    }
    U = (*told_ptr*(1.0 - CXPC - CXP - CZP)
            - *(&U - stepsize[0])*(CXMC + CXM)
            + *(told_ptr + stepsize[0])*(CXPC + CXP)
            - *(&U - stepsize[2])*CZM
            + *(told_ptr + stepsize[2])*CZP
            + Q) /
         (1.0 - CXMC - CXM - CZM);
  }
  else
  {
    U = (*told_ptr*(1.0 - CZP)
            - *(&U - stepsize[2])*CZM
            + *(told_ptr + stepsize[2])*CZP
            + Q) /
        (1.0 - CZM);
  }
}

void Cell::calcCellADEDown(double timestep, const Foundation &foundation, const BoundaryConditions &/*bcs*/,
                           double &V)
{
  double theta = timestep/
                 (density*specificHeat);

  double CXP = pde[0][1]*theta;
  double CXM = pde[0][0]*theta;
  double CZP = pde[2][1]*theta;
  double CZM = pde[2][0]*theta;
  double Q = heatGain*theta;

  if (foundation.numberOfDimensions == 3) {
    double CYP = pde[1][1] * theta;
    double CYM = pde[1][0] * theta;
    V = (*told_ptr * (1.0 + CXM + CZM + CYM)
                - *(told_ptr - stepsize[0]) * CXM
                + *(&V + stepsize[0]) * CXP
                - *(told_ptr - stepsize[2]) * CZM
                + *(&V + stepsize[2]) * CZP
                - *(told_ptr - stepsize[1]) * CYM
                + *(&V + stepsize[1]) * CYP
                + Q) /
               (1.0 + CXP + CZP + CYP);
  } else if (foundation.numberOfDimensions == 2)
  {
    double CXPC = 0;
    double CXMC = 0;

    if (coords[0] != 0)
    {
      CXPC = pde_c[1]*theta/r;
      CXMC = pde_c[0]*theta/r;
    }
    V = (*told_ptr*(1.0 + CXMC + CXM + CZM)
                - *(told_ptr - stepsize[0])*(CXMC + CXM)
                + *(&V + stepsize[0])*(CXPC + CXP)
                - *(told_ptr - stepsize[2])*CZM
                + *(&V + stepsize[2])*CZP
                + Q) /
               (1.0 + CXPC + CXP + CZP);
  }
  else
  {
    V = (*told_ptr*(1.0 + CZM)
                - *(told_ptr - stepsize[2])*CZM
                + *(&V + stepsize[2])*CZP
                + Q) /
               (1.0 + CZP);
  }
}

double Cell::calcCellExplicit(double timestep, const Foundation &foundation,
                              const BoundaryConditions &/*bcs*/)
{
  double theta = timestep/
                 (density*specificHeat);

  double CXP = pde[0][1]*theta;
  double CXM = pde[0][0]*theta;
  double CZP = pde[2][1]*theta;
  double CZM = pde[2][0]*theta;
  double Q = heatGain*theta;

  if (foundation.numberOfDimensions == 3) {
    double CYP = pde[1][1]*theta;
    double CYM = pde[1][0]*theta;
    double TNew = *told_ptr * (1.0 + CXM + CZM + CYM - CXP - CZP - CYP)
                  - *(told_ptr - stepsize[0]) * CXM
                  + *(told_ptr + stepsize[0]) * CXP
                  - *(told_ptr - stepsize[2]) * CZM
                  + *(told_ptr + stepsize[2]) * CZP
                  - *(told_ptr - stepsize[1]) * CYM
                  + *(told_ptr + stepsize[1]) * CYP
                  + Q;
    return TNew;
  }
  if (foundation.numberOfDimensions == 2)
  {
    double CXPC = 0;
    double CXMC = 0;

    if (coords[0] != 0)
    {
      CXPC = pde_c[1]*theta/r;
      CXMC = pde_c[0]*theta/r;
    }

    double TNew = *told_ptr*(1.0 + CXMC + CXM + CZM - CXPC - CXP - CZP)
                  - *(told_ptr - stepsize[0])*(CXMC + CXM)
                  + *(told_ptr + stepsize[0])*(CXPC + CXP)
                  - *(told_ptr - stepsize[2])*CZM
                  + *(told_ptr + stepsize[2])*CZP
                  + Q;
    return TNew;
  }

  double TNew = *told_ptr*(1.0 + CZM - CZP)
                - *(told_ptr - stepsize[2])*CZM
                + *(told_ptr + stepsize[2])*CZP
                + Q;
  return TNew;
}

void Cell::calcCellMatrix(Foundation::NumericalScheme scheme, const double &timestep,
                          const Foundation &foundation, const BoundaryConditions &/*bcs*/,
                          double &A, double (&Alt)[3][2], double &bVal) {
  if (scheme == Foundation::NS_STEADY_STATE) {
    calcCellSteadyState(foundation, A, Alt, bVal);
  } else {
    double theta = timestep/(density*specificHeat);

    double f = scheme == Foundation::NS_IMPLICIT? 1.0 : 0.5;

    std::size_t dims[3]{0, 1, 2};
    if (foundation.numberOfDimensions == 2) {
      dims[1] = 5;
    } else if (foundation.numberOfDimensions == 1) {
      dims[0] = 5;
      dims[1] = 5;
    }

    bool cylindrical = (foundation.coordinateSystem == Foundation::CS_CYLINDRICAL);
    double C[3][2]{{0}};
    gatherCCoeffs(dims, theta, cylindrical, C);

    double bit{0};
    bVal = heatGain * theta;
    for (auto dim : dims) {
      if (dim < 5) {
        bit += C[dim][1] - C[dim][0];
        Alt[dim][1] = -f * C[dim][1];
        Alt[dim][0] = f * C[dim][0];
        bVal += *(told_ptr + stepsize[dim]) * (1 - f) * C[dim][1] -
                *(told_ptr - stepsize[dim]) * (1 - f) * C[dim][0];
      }
    }
    A = (1.0 + f*bit);
    bVal += *told_ptr*(1.0 - (1-f)*bit);
  }
}

void Cell::calcCellSteadyState(const Foundation &foundation,
                               double &A, double (&Alt)[3][2], double &bVal)
{
  std::size_t dims[3]{0, 1, 2};
  if (foundation.numberOfDimensions == 2) {
    dims[1] = 5;
  } else if (foundation.numberOfDimensions == 1) {
    dims[0] = 5;
    dims[1] = 5;
  }

  A = 0;
  for (auto dim : dims) {
    if (dim < 5) {
      Alt[dim][1] = pde[dim][1];
      Alt[dim][0] = -pde[dim][0];
      A += Alt[dim][1] + Alt[dim][0];
    }
  }
  if (foundation.coordinateSystem == Foundation::CS_CYLINDRICAL &&
      coords[0] != 0) {
    Alt[0][1] += pde_c[1]/r;
    Alt[0][0] += -pde_c[0]/r;
    A += pde_c[1]/r - pde_c[0]/r;
  }
  A *= -1;
  bVal = -heatGain;
}

void Cell::calcCellADI(int dim, const double &timestep,
                       const Foundation &foundation, const BoundaryConditions &/*bcs*/,
                       double &A, double (&Alt)[2], double &bVal) {
  double theta = timestep / (foundation.numberOfDimensions
                             *density*specificHeat);
  double Q = heatGain*theta;

  if (foundation.numberOfDimensions == 1) {
    A = 1.0 + (pde[2][1] - pde[2][0])*theta;
    Alt[0] = pde[2][0]*theta;
    Alt[1] = -pde[2][1]*theta;

    bVal = *told_ptr + Q;
    return;
  }

  size_t dims[3]{0, 1, 2};
  if (foundation.numberOfDimensions == 2) {
    dims[1] = 5;
  }

  double f = foundation.fADI;
  double multiplier = foundation.numberOfDimensions == 2 ? (2.0 - f) : (3.0-2.0*f);
  double C[3][2]{{0}};
  gatherCCoeffs(dims, theta, foundation.coordinateSystem == Foundation::CS_CYLINDRICAL, C);

  ADImath(dim, dims, Q, f, multiplier, C, A, Alt, bVal);
}

void Cell::ADImath(int dim, const std::size_t (&dims)[3],
                   const double Q, const double f, const double multiplier, const double (&C)[3][2],
                   double &A, double (&Alt)[2], double &bVal) {
  bVal = Q;
  double bit{0};
  for (auto sdim : dims) {
    if (sdim == dim) {
      Alt[0] = multiplier*C[sdim][0];
      Alt[1] = -multiplier*C[sdim][1];
      A = 1.0 - (Alt[0] + Alt[1]);
    } else if (sdim < 5) {
      bit += C[sdim][0] - C[sdim][1];
      bVal += *(told_ptr + stepsize[sdim]) * f * C[sdim][1] -
              *(told_ptr - stepsize[sdim]) * f * C[sdim][0];
    }
  }
  bVal += *told_ptr * (1.0 + f * bit);
}


void Cell::gatherCCoeffs(const std::size_t (&dims)[3], const double &theta,
                         bool cylindrical, double (&C)[3][2])
{
  for (auto dim : dims) {
    if (dim < 5) {
      C[dim][0] = pde[dim][0] * theta;
      C[dim][1] = pde[dim][1] * theta;
    }
  }
  if (cylindrical && coords[0] != 0) {
    C[0][0] += pde_c[0]*theta/r;
    C[0][1] += pde_c[1]*theta/r;
  }
}

std::vector<double> Cell::calculateHeatFlux(int ndims, double &TNew,
                                            std::size_t nX, std::size_t nY, std::size_t nZ,
                                            const std::vector< std::shared_ptr<Cell> > &/*cell_v*/)
{
  std::vector<double> Qflux;
  double Qx = 0;
  double Qy = 0;
  double Qz = 0;

  double CXP = 0;
  double CXM = 0;
  double CYP = 0;
  double CYM = 0;
  double CZP = -kcoeff[2][1]*dist[2][0]/(dist[2][1]+dist[2][0])/dist[2][1];
  double CZM = -kcoeff[2][0]*dist[2][1]/(dist[2][1]+dist[2][0])/dist[2][0];

  if (ndims > 1)
  {
    CXP = -kcoeff[0][1]*dist[0][0]/(dist[0][1]+dist[0][0])/dist[0][1];
    CXM = -kcoeff[0][0]*dist[0][1]/(dist[0][1]+dist[0][0])/dist[0][0];
  }


  if (ndims == 3)
  {
    CYP = -kcoeff[1][1]*dist[1][0]/(dist[1][1]+dist[1][0])/dist[1][1];
    CYM = -kcoeff[1][0]*dist[1][1]/(dist[1][1]+dist[1][0])/dist[1][0];
  }

  double DTXP = 0;
  double DTXM = 0;
  double DTYP = 0;
  double DTYM = 0;
  double DTZP = 0;
  double DTZM = 0;

  if (coords[0] != nX - 1)
    DTXP = *(&TNew + stepsize[0])-TNew;

  if (coords[0] != 0)
    DTXM = TNew-*(&TNew - stepsize[0]);

  if (coords[1] != nY - 1)
    DTYP = *(&TNew + stepsize[1])-TNew;

  if (coords[1] != 0)
    DTYM = TNew-*(&TNew - stepsize[1]);

  if (coords[2] != nZ - 1)
    DTZP = *(&TNew + stepsize[2])-TNew;

  if (coords[2] != 0)
    DTZM = TNew-*(&TNew - stepsize[2]);

  Qx = CXP*DTXP + CXM*DTXM;
  Qy = CYP*DTYP + CYM*DTYM;
  Qz = CZP*DTZP + CZM*DTZM;

  Qflux.push_back(Qx);
  Qflux.push_back(Qy);
  Qflux.push_back(Qz);

  return Qflux;
}


void Cell::doOutdoorTemp(const BoundaryConditions &bcs, double &A, double &bVal)
{
  A = 1.0;
  bVal = bcs.outdoorTemp;
}

void Cell::doIndoorTemp(const BoundaryConditions &bcs, double &A, double &bVal)
{
  A = 1.0;
  bVal = bcs.indoorTemp;
}


ExteriorAirCell::ExteriorAirCell(const std::size_t &index, const CellType cellType,
                                 const std::size_t &i, const std::size_t &j, const std::size_t &k,
                                 std::size_t *stepsize,
                                 const Foundation &foundation, Surface *surfacePtr, Block *blockPtr,
                                 Mesher *meshPtr):
        Cell(index, cellType, i, j, k, stepsize, foundation, surfacePtr, blockPtr, meshPtr)
{}

void ExteriorAirCell::calcCellADEUp(double /*timestep*/, const Foundation &/*foundation*/, const BoundaryConditions &bcs,
                                    double &U)
{
  U = bcs.outdoorTemp;
}

void ExteriorAirCell::calcCellADEDown(double /*timestep*/, const Foundation &/*foundation*/, const BoundaryConditions &bcs,
                                     double &V)
{
  V = bcs.outdoorTemp;
}

double ExteriorAirCell::calcCellExplicit(double /*timestep*/, const Foundation &/*foundation*/,
                                         const BoundaryConditions &bcs)
{
  return bcs.outdoorTemp;
}

void ExteriorAirCell::calcCellADI(int /*dim*/, const double &/*timestep*/,
                                  const Foundation &, const BoundaryConditions &bcs,
                                  double &A, double (&)[2], double &bVal)
{
  doOutdoorTemp(bcs, A, bVal);
};


void ExteriorAirCell::calcCellMatrix(Foundation::NumericalScheme, const double &/*timestep*/,
                                     const Foundation &, const BoundaryConditions &bcs,
                                     double &A, double (&)[3][2], double &bVal)
{
  doOutdoorTemp(bcs, A, bVal);
}


std::vector<double> ExteriorAirCell::calculateHeatFlux(int /*ndims*/, double &/*TNew*/,
                                            std::size_t /*nX*/, std::size_t /*nY*/, std::size_t /*nZ*/,
                                                       const std::vector< std::shared_ptr<Cell> > &/*cell_v*/)
{
  std::vector<double> Qflux;
  Qflux.push_back(0);
  Qflux.push_back(0);
  Qflux.push_back(0);
  return Qflux;
}



InteriorAirCell::InteriorAirCell(const std::size_t &index, const CellType cellType,
                                 const std::size_t &i, const std::size_t &j, const std::size_t &k,
                                 std::size_t *stepsize,
                                 const Foundation &foundation, Surface *surfacePtr, Block *blockPtr,
                                 Mesher *meshPtr):
        Cell(index, cellType, i, j, k, stepsize, foundation, surfacePtr, blockPtr, meshPtr)
{}

void InteriorAirCell::calcCellADEUp(double /*timestep*/, const Foundation &/*foundation*/, const BoundaryConditions &bcs,
                                    double &U)
{
  U = bcs.indoorTemp;
}

void InteriorAirCell::calcCellADEDown(double /*timestep*/, const Foundation &/*foundation*/, const BoundaryConditions &bcs,
                                      double &V)
{
  V = bcs.indoorTemp;
}

double InteriorAirCell::calcCellExplicit(double /*timestep*/, const Foundation &/*foundation*/,
                                         const BoundaryConditions &bcs)
{
  return bcs.indoorTemp;
}

void InteriorAirCell::calcCellMatrix(Foundation::NumericalScheme, const double &/*timestep*/,
                                     const Foundation &, const BoundaryConditions &bcs,
                                     double &A, double (&)[3][2], double &bVal)
{
  doIndoorTemp(bcs, A, bVal);
}

void InteriorAirCell::calcCellADI(int /*dim*/, const double &/*timestep*/,
                                  const Foundation &, const BoundaryConditions &bcs,
                                  double &A, double (&)[2], double &bVal)
{
  doIndoorTemp(bcs, A, bVal);
};

std::vector<double> InteriorAirCell::calculateHeatFlux(int /*ndims*/, double &/*TNew*/,
                                                       std::size_t /*nX*/, std::size_t /*nY*/, std::size_t /*nZ*/,
                                                       const std::vector< std::shared_ptr<Cell> > &/*cell_v*/)
{
  std::vector<double> Qflux;
  Qflux.push_back(0);
  Qflux.push_back(0);
  Qflux.push_back(0);
  return Qflux;
}


BoundaryCell::BoundaryCell(const std::size_t &index, const CellType cellType,
                           const std::size_t &i, const std::size_t &j, const std::size_t &k,
                           std::size_t *stepsize,
                           const Foundation &foundation, Surface *surfacePtr, Block *blockPtr,
                           Mesher *meshPtr):
        Cell(index, cellType, i, j, k, stepsize, foundation, surfacePtr, blockPtr, meshPtr)
{
    if (foundation.numberOfDimensions == 2 &&
        foundation.coordinateSystem == Foundation::CS_CYLINDRICAL)
    {
      if (surfacePtr->orientation == Surface::X_POS ||
          surfacePtr->orientation == Surface::X_NEG)
      {
        area = 2.0 * PI * meshPtr[0].centers[coords[0]] * meshPtr[2].deltas[coords[2]];
      }
      else // if (surface.orientation == Surface::Z_POS ||
        // surface.orientation == Surface::Z_NEG)
      {
        area = PI * (meshPtr[0].dividers[coords[0] + 1] * meshPtr[0].dividers[coords[0] + 1] -
                     meshPtr[0].dividers[coords[0]]*meshPtr[0].dividers[coords[0]] );
      }
    }
    else if (foundation.numberOfDimensions == 2 &&
             foundation.coordinateSystem == Foundation::CS_CARTESIAN)
    {
      if (surfacePtr->orientation == Surface::X_POS ||
          surfacePtr->orientation == Surface::X_NEG)
      {
        area = 2.0 * meshPtr[2].deltas[coords[2]] * foundation.linearAreaMultiplier;
      }
      else // if (surface.orientation == Surface::Z_POS ||
        // surface.orientation == Surface::Z_NEG)
      {
        area = 2.0 * meshPtr[0].deltas[coords[0]] * foundation.linearAreaMultiplier;
      }
    }
    else if (foundation.numberOfDimensions == 3)
    {
      if (surfacePtr->orientation == Surface::X_POS ||
          surfacePtr->orientation == Surface::X_NEG)
      {
        area = meshPtr[1].deltas[coords[1]] * meshPtr[2].deltas[coords[2]];
      }
      else if (surfacePtr->orientation == Surface::Y_POS ||
               surfacePtr->orientation == Surface::Y_NEG)
      {
        area = meshPtr[0].deltas[coords[0]] * meshPtr[2].deltas[coords[2]];
      }
      else // if (surface.orientation == Surface::Z_POS ||
        // surface.orientation == Surface::Z_NEG)
      {
        area = meshPtr[0].deltas[coords[0]] * meshPtr[1].deltas[coords[1]];
      }

      if (foundation.useSymmetry)
      {
        if (foundation.isXSymm)
          area = 2 * area;

        if (foundation.isYSymm)
          area = 2 * area;
      }
    }
    else  /* if (foundation.numberOfDimensions == 1) */
    {
      area = 1.0;
    }
}

void BoundaryCell::calcCellADEUp(double /*timestep*/, const Foundation &foundation, const BoundaryConditions &bcs,
                                    double &U)
{
  switch (surfacePtr->boundaryConditionType)
  {
    case Surface::ZERO_FLUX:
    {
      switch (surfacePtr->orientation)
      {
        case Surface::X_NEG:
          U = *(told_ptr + stepsize[0]);
          break;
        case Surface::X_POS:
          U = *(&U - stepsize[0]);
          break;
        case Surface::Y_NEG:
          U = *(told_ptr + stepsize[1]);
          break;
        case Surface::Y_POS:
          U = *(&U - stepsize[1]);
          break;
        case Surface::Z_NEG:
          U = *(told_ptr + stepsize[2]);
          break;
        case Surface::Z_POS:
          U = *(&U - stepsize[2]);
          break;
      }
    }
      break;

    case Surface::CONSTANT_TEMPERATURE:

      U = surfacePtr->temperature;
      break;

    case Surface::INTERIOR_TEMPERATURE:

      U = bcs.indoorTemp;
      break;

    case Surface::EXTERIOR_TEMPERATURE:

      U = bcs.outdoorTemp;
      break;

    case Surface::INTERIOR_FLUX:
    {
      double Tair = bcs.indoorTemp;
      double q = heatGain;

      double hc = foundation.getConvectionCoeff(*told_ptr,
                                     Tair,0.0,0.00208,false,surfacePtr->tilt);  // TODO Make roughness a property of the interior surfaces
      double hr = getSimpleInteriorIRCoeff(surfacePtr->emissivity,
                                           *told_ptr,Tair);

      switch (surfacePtr->orientation)
      {
        case Surface::X_NEG:
          U = (kcoeff[0][1] * *(told_ptr + stepsize[0]) / dist[0][1] +
                      (hc + hr)*Tair + q)/(kcoeff[0][1]/dist[0][1] + (hc + hr));
          break;
        case Surface::X_POS:
          U = (kcoeff[0][0] * *(&U - stepsize[0])/dist[0][0] +
                      (hc + hr)*Tair + q)/(kcoeff[0][0]/dist[0][0] + (hc + hr));
          break;
        case Surface::Y_NEG:
          U = (kcoeff[1][1] * *(told_ptr + stepsize[1]) / dist[1][1] +
                      (hc + hr)*Tair + q)/(kcoeff[1][1]/dist[1][1] + (hc + hr));
          break;
        case Surface::Y_POS:
          U = (kcoeff[1][0] * *(&U - stepsize[1])/dist[1][0] +
                      (hc + hr)*Tair + q)/(kcoeff[1][0]/dist[1][0] + (hc + hr));
          break;
        case Surface::Z_NEG:
          U = (kcoeff[2][1] * *(told_ptr + stepsize[2]) / dist[2][1] +
                      (hc + hr)*Tair + q)/(kcoeff[2][1]/dist[2][1] + (hc + hr));
          break;
        case Surface::Z_POS:
          U = (kcoeff[2][0] * *(&U - stepsize[2])/dist[2][0] +
                      (hc + hr)*Tair + q)/(kcoeff[2][0]/dist[2][0] + (hc + hr));
          break;
      }
    }
      break;

    case Surface::EXTERIOR_FLUX:
    {
      double Tair = bcs.outdoorTemp;
      double v = bcs.localWindSpeed;
      double eSky = bcs.skyEmissivity;
      double tilt = surfacePtr->tilt;
      double F = getEffectiveExteriorViewFactor(eSky,tilt);
      double hc = foundation.getConvectionCoeff(*told_ptr,Tair,v,foundation.surfaceRoughness,true,tilt);
      double hr = getExteriorIRCoeff(surfacePtr->emissivity,*told_ptr,Tair,eSky,tilt);
      double q = heatGain;

      switch (surfacePtr->orientation)
      {
        case Surface::X_NEG:
          U = (kcoeff[0][1] * *(told_ptr + stepsize[0]) / dist[0][1] +
                      (hc + hr*pow(F,0.25))*Tair + q)/(kcoeff[0][1]/dist[0][1] + (hc + hr));
          break;
        case Surface::X_POS:
          U = (kcoeff[0][0] * *(&U - stepsize[0])/dist[0][0] +
                      (hc + hr*pow(F,0.25))*Tair + q)/(kcoeff[0][0]/dist[0][0] + (hc + hr));
          break;
        case Surface::Y_NEG:
          U = (kcoeff[1][1] * *(told_ptr + stepsize[1]) / dist[1][1] +
                      (hc + hr*pow(F,0.25))*Tair + q)/(kcoeff[1][1]/dist[1][1] + (hc + hr));
          break;
        case Surface::Y_POS:
          U = (kcoeff[1][0] * *(&U - stepsize[1])/dist[1][0] +
                      (hc + hr*pow(F,0.25))*Tair + q)/(kcoeff[1][0]/dist[1][0] + (hc + hr));
          break;
        case Surface::Z_NEG:
          U = (kcoeff[2][1] * *(told_ptr + stepsize[2]) / dist[2][1] +
                      (hc + hr*pow(F,0.25))*Tair + q)/(kcoeff[2][1]/dist[2][1] + (hc + hr));
          break;
        case Surface::Z_POS:
          U = (kcoeff[2][0] * *(&U - stepsize[2])/dist[2][0] +
                      (hc + hr*pow(F,0.25))*Tair + q)/(kcoeff[2][0]/dist[2][0] + (hc + hr));
          break;
      }
    }
      break;
  }

}

void BoundaryCell::calcCellADEDown(double /*timestep*/, const Foundation &foundation, const BoundaryConditions &bcs,
                                 double &V)
{
  switch (surfacePtr->boundaryConditionType)
  {
    case Surface::ZERO_FLUX:
    {
      switch (surfacePtr->orientation)
      {
        case Surface::X_NEG:
          V = *(&V + stepsize[0]);
          break;
        case Surface::X_POS:
          V = *(told_ptr - stepsize[0]);
          break;
        case Surface::Y_NEG:
          V = *(&V + stepsize[1]);
          break;
        case Surface::Y_POS:
          V = *(told_ptr - stepsize[1]);
          break;
        case Surface::Z_NEG:
          V = *(&V + stepsize[2]);
          break;
        case Surface::Z_POS:
          V = *(told_ptr - stepsize[2]);
          break;
      }
    }
      break;

    case Surface::CONSTANT_TEMPERATURE:

      V = surfacePtr->temperature;
      break;

    case Surface::INTERIOR_TEMPERATURE:

      V = bcs.indoorTemp;
      break;

    case Surface::EXTERIOR_TEMPERATURE:

      V = bcs.outdoorTemp;
      break;

    case Surface::INTERIOR_FLUX:
    {
      double Tair = bcs.indoorTemp;
      double q = heatGain;

      double hc = foundation.getConvectionCoeff(*told_ptr,
                                     Tair,0.0,0.00208,false,surfacePtr->tilt);
      double hr = getSimpleInteriorIRCoeff(surfacePtr->emissivity,
                                           *told_ptr,Tair);

      switch (surfacePtr->orientation)
      {
        case Surface::X_NEG:
          V = (kcoeff[0][1] * *(&V + stepsize[0])/dist[0][1] +
                      (hc + hr)*Tair + q)/(kcoeff[0][1]/dist[0][1] + (hc + hr));
          break;
        case Surface::X_POS:
          V = (kcoeff[0][0] * *(told_ptr - stepsize[0]) / dist[0][0] +
                      (hc + hr)*Tair + q)/(kcoeff[0][0]/dist[0][0] + (hc + hr));
          break;
        case Surface::Y_NEG:
          V = (kcoeff[1][1] * *(&V + stepsize[1])/dist[1][1] +
                      (hc + hr)*Tair + q)/(kcoeff[1][1]/dist[1][1] + (hc + hr));
          break;
        case Surface::Y_POS:
          V = (kcoeff[1][0] * *(told_ptr - stepsize[1]) / dist[1][0] +
                      (hc + hr)*Tair + q)/(kcoeff[1][0]/dist[1][0] + (hc + hr));
          break;
        case Surface::Z_NEG:
          V = (kcoeff[2][1] * *(&V + stepsize[2])/dist[2][1] +
                      (hc + hr)*Tair + q)/(kcoeff[2][1]/dist[2][1] + (hc + hr));
          break;
        case Surface::Z_POS:
          V = (kcoeff[2][0] * *(told_ptr - stepsize[2]) / dist[2][0] +
                      (hc + hr)*Tair + q)/(kcoeff[2][0]/dist[2][0] + (hc + hr));
          break;
      }
    }
      break;

    case Surface::EXTERIOR_FLUX:
    {
      double Tair = bcs.outdoorTemp;
      double v = bcs.localWindSpeed;
      double eSky = bcs.skyEmissivity;
      double tilt = surfacePtr->tilt;
      double F = getEffectiveExteriorViewFactor(eSky,tilt);
      double hc = foundation.getConvectionCoeff(*told_ptr,Tair,v,foundation.surfaceRoughness,true,tilt);
      double hr = getExteriorIRCoeff(surfacePtr->emissivity,*told_ptr,Tair,eSky,tilt);
      double q = heatGain;

      switch (surfacePtr->orientation)
      {
        case Surface::X_NEG:
          V = (kcoeff[0][1] * *(&V + stepsize[0])/dist[0][1] +
                      (hc + hr*pow(F,0.25))*Tair + q)/(kcoeff[0][1]/dist[0][1] + (hc + hr));
          break;
        case Surface::X_POS:
          V = (kcoeff[0][0] * *(told_ptr - stepsize[0]) / dist[0][0] +
                      (hc + hr*pow(F,0.25))*Tair + q)/(kcoeff[0][0]/dist[0][0] + (hc + hr));
          break;
        case Surface::Y_NEG:
          V = (kcoeff[1][1] * *(&V + stepsize[1])/dist[1][1] +
                      (hc + hr*pow(F,0.25))*Tair + q)/(kcoeff[1][1]/dist[1][1] + (hc + hr));
          break;
        case Surface::Y_POS:
          V = (kcoeff[1][0] * *(told_ptr - stepsize[1]) / dist[1][0] +
                      (hc + hr*pow(F,0.25))*Tair + q)/(kcoeff[1][0]/dist[1][0] + (hc + hr));
          break;
        case Surface::Z_NEG:
          V = (kcoeff[2][1] * *(&V + stepsize[2])/dist[2][1] +
                      (hc + hr*pow(F,0.25))*Tair + q)/(kcoeff[2][1]/dist[2][1] + (hc + hr));
          break;
        case Surface::Z_POS:
          V = (kcoeff[2][0] * *(told_ptr - stepsize[2]) / dist[2][0] +
                      (hc + hr*pow(F,0.25))*Tair + q)/(kcoeff[2][0]/dist[2][0] + (hc + hr));
          break;
      }
    }
      break;
  }
}

double BoundaryCell::calcCellExplicit(double /*timestep*/, const Foundation &foundation,
                                      const BoundaryConditions &bcs)
{
  switch (surfacePtr->boundaryConditionType)
  {
    case Surface::ZERO_FLUX:
    {
      switch (surfacePtr->orientation)
      {
        case Surface::X_NEG:
          return *(told_ptr + stepsize[0]);
        case Surface::X_POS:
          return *(told_ptr - stepsize[0]);
        case Surface::Y_NEG:
          return *(told_ptr + stepsize[1]);
        case Surface::Y_POS:
          return *(told_ptr - stepsize[1]);
        case Surface::Z_NEG:
          return *(told_ptr + stepsize[2]);
        case Surface::Z_POS:
          return *(told_ptr - stepsize[2]);
      }
    }
    break;

    case Surface::CONSTANT_TEMPERATURE:
      return surfacePtr->temperature;

    case Surface::INTERIOR_TEMPERATURE:
      return bcs.indoorTemp;

    case Surface::EXTERIOR_TEMPERATURE:
      return bcs.outdoorTemp;

    case Surface::INTERIOR_FLUX:
    {
      double Tair = bcs.indoorTemp;
      double& q = heatGain;

      double hc = foundation.getConvectionCoeff(*told_ptr,
                                     Tair,0.0,0.00208,false,surfacePtr->tilt);
      double hr = getSimpleInteriorIRCoeff(surfacePtr->emissivity,
                                           *told_ptr,Tair);

      switch (surfacePtr->orientation)
      {
        case Surface::X_NEG:
          return (kcoeff[0][1] * *(told_ptr + stepsize[0])/dist[0][1] +
                         (hc + hr)*Tair + q)/(kcoeff[0][1]/dist[0][1] + (hc + hr));
        case Surface::X_POS:
          return (kcoeff[0][0] * *(told_ptr - stepsize[0])/dist[0][0] +
                         (hc + hr)*Tair + q)/(kcoeff[0][0]/dist[0][0] + (hc + hr));
        case Surface::Y_NEG:
          return (kcoeff[1][1] * *(told_ptr + stepsize[1])/dist[1][1] +
                         (hc + hr)*Tair + q)/(kcoeff[1][1]/dist[1][1] + (hc + hr));
        case Surface::Y_POS:
          return (kcoeff[1][0] * *(told_ptr - stepsize[1])/dist[1][0] +
                         (hc + hr)*Tair + q)/(kcoeff[1][0]/dist[1][0] + (hc + hr));
        case Surface::Z_NEG:
          return (kcoeff[2][1] * *(told_ptr + stepsize[2])/dist[2][1] +
                         (hc + hr)*Tair + q)/(kcoeff[2][1]/dist[2][1] + (hc + hr));
        case Surface::Z_POS:
          return (kcoeff[2][0] * *(told_ptr - stepsize[2])/dist[2][0] +
                         (hc + hr)*Tair + q)/(kcoeff[2][0]/dist[2][0] + (hc + hr));
      }
    }
    break;

    case Surface::EXTERIOR_FLUX:
    {
      double Tair = bcs.outdoorTemp;
      double v = bcs.localWindSpeed;
      double eSky = bcs.skyEmissivity;
      double tilt = surfacePtr->tilt;
      double F = getEffectiveExteriorViewFactor(eSky,tilt);
      double hc = foundation.getConvectionCoeff(*told_ptr,Tair,v,foundation.surfaceRoughness,true,tilt);
      double hr = getExteriorIRCoeff(surfacePtr->emissivity,*told_ptr,Tair,eSky,tilt);
      double& q = heatGain;

      switch (surfacePtr->orientation)
      {
        case Surface::X_NEG:
          return (kcoeff[0][1] * *(told_ptr + stepsize[0])/dist[0][1] +
                         (hc + hr*pow(F,0.25))*Tair + q)/(kcoeff[0][1]/dist[0][1] + (hc + hr));
        case Surface::X_POS:
          return (kcoeff[0][0] * *(told_ptr - stepsize[0])/dist[0][0] +
                         (hc + hr*pow(F,0.25))*Tair + q)/(kcoeff[0][0]/dist[0][0] + (hc + hr));
        case Surface::Y_NEG:
          return (kcoeff[1][1] * *(told_ptr + stepsize[1])/dist[1][1] +
                         (hc + hr*pow(F,0.25))*Tair + q)/(kcoeff[1][1]/dist[1][1] + (hc + hr));
        case Surface::Y_POS:
          return (kcoeff[1][0] * *(told_ptr - stepsize[1])/dist[1][0] +
                         (hc + hr*pow(F,0.25))*Tair + q)/(kcoeff[1][0]/dist[1][0] + (hc + hr));
        case Surface::Z_NEG:
          return (kcoeff[2][1] * *(told_ptr + stepsize[2])/dist[2][1] +
                         (hc + hr*pow(F,0.25))*Tair + q)/(kcoeff[2][1]/dist[2][1] + (hc + hr));
        case Surface::Z_POS:
          return (kcoeff[2][0] * *(told_ptr - stepsize[2])/dist[2][0] +
                         (hc + hr*pow(F,0.25))*Tair + q)/(kcoeff[2][0]/dist[2][0] + (hc + hr));
      }
    }
    break;
  }
}

void BoundaryCell::calcCellADI(int dim, const double &/*timestep*/,
                               const Foundation &foundation, const BoundaryConditions &bcs,
                               double &A, double (&Alt)[2], double &bVal)
{
  std::size_t sdim = surfacePtr->orientation_dim;
  std::size_t dir = surfacePtr->orientation_dir;

  switch (surfacePtr->boundaryConditionType) {
    case Surface::ZERO_FLUX:
      zfCellADI(dim, sdim, dir, A, Alt[dir], bVal);
      break;
    case Surface::CONSTANT_TEMPERATURE:
      A = 1.0;
      bVal = surfacePtr->temperature;
      break;
    case Surface::INTERIOR_TEMPERATURE:
      doIndoorTemp(bcs, A, bVal);
      break;
    case Surface::EXTERIOR_TEMPERATURE:
      doOutdoorTemp(bcs, A, bVal);
      break;
    case Surface::INTERIOR_FLUX:
      ifCellADI(dim, sdim, dir, foundation, bcs, A, Alt[dir], bVal);
      break;
    case Surface::EXTERIOR_FLUX:
      efCellADI(dim, sdim, dir, foundation, bcs, A, Alt[dir], bVal);
      break;
  }
}

void BoundaryCell::calcCellMatrix(Foundation::NumericalScheme, const double &/*timestep*/,
                                  const Foundation &foundation, const BoundaryConditions &bcs,
                                  double &A, double (&Alt)[3][2], double &bVal)
{
  std::size_t dim = surfacePtr->orientation_dim;
  std::size_t dir = surfacePtr->orientation_dir;

  switch (surfacePtr->boundaryConditionType)
  {
    case Surface::ZERO_FLUX: {
      zfCellMatrix(A, Alt[dim][dir], bVal);
      break;
    }
    case Surface::CONSTANT_TEMPERATURE: {
      A = 1.0;
      bVal = surfacePtr->temperature;
      break;
    }
    case Surface::INTERIOR_TEMPERATURE: {
      doIndoorTemp(bcs, A, bVal);
      break;
    }
    case Surface::EXTERIOR_TEMPERATURE: {
      doOutdoorTemp(bcs, A, bVal);
      break;
    }
    case Surface::INTERIOR_FLUX: {
      ifCellMatrix(dim, dir, foundation, bcs, A, Alt[dim][dir], bVal);
      break;
    }
    case Surface::EXTERIOR_FLUX: {
      efCellMatrix(dim, dir, foundation, bcs, A, Alt[dim][dir], bVal);
      break;
    }
  }
}

std::vector<double> BoundaryCell::calculateHeatFlux(int ndims, double &TNew,
                                                    std::size_t nX, std::size_t nY, std::size_t nZ,
                                                    const std::vector< std::shared_ptr<Cell> > &/*cell_v*/)
{
  std::vector<double> Qflux;
  double Qx = 0;
  double Qy = 0;
  double Qz = 0;

  double CXP = 0;
  double CXM = 0;
  double CYP = 0;
  double CYM = 0;
  double CZP = -kcoeff[2][1]*dist[2][0]/(dist[2][1]+dist[2][0])/dist[2][1];
  double CZM = -kcoeff[2][0]*dist[2][1]/(dist[2][1]+dist[2][0])/dist[2][0];

  if (ndims > 1)
  {
    CXP = -kcoeff[0][1]*dist[0][0]/(dist[0][1]+dist[0][0])/dist[0][1];
    CXM = -kcoeff[0][0]*dist[0][1]/(dist[0][1]+dist[0][0])/dist[0][0];
  }


  if (ndims == 3)
  {
    CYP = -kcoeff[1][1]*dist[1][0]/(dist[1][1]+dist[1][0])/dist[1][1];
    CYM = -kcoeff[1][0]*dist[1][1]/(dist[1][1]+dist[1][0])/dist[1][0];
  }

  double DTXP = 0;
  double DTXM = 0;
  double DTYP = 0;
  double DTYM = 0;
  double DTZP = 0;
  double DTZM = 0;

  if (coords[0] != nX - 1)
    DTXP = *(&TNew + stepsize[0])-TNew;

  if (coords[0] != 0)
    DTXM = TNew-*(&TNew - stepsize[0]);

  if (coords[1] != nY - 1)
    DTYP = *(&TNew + stepsize[1])-TNew;

  if (coords[1] != 0)
    DTYM = TNew-*(&TNew - stepsize[1]);

  if (coords[2] != nZ - 1)
    DTZP = *(&TNew + stepsize[2])-TNew;

  if (coords[2] != 0)
    DTZM = TNew-*(&TNew - stepsize[2]);

  switch (surfacePtr->orientation)
  {
    case Surface::X_NEG:
    {
      CXP = -kcoeff[0][1]/dist[0][1];
      CXM = 0;
    }
      break;
    case Surface::X_POS:
    {
      CXP = 0;
      CXM = -kcoeff[0][0]/dist[0][0];
    }
      break;
    case Surface::Y_NEG:
    {
      CYP = -kcoeff[1][1]/dist[1][1];
      CYM = 0;
    }
      break;
    case Surface::Y_POS:
    {
      CYP = 0;
      CYM = -kcoeff[1][0]/dist[1][0];
    }
      break;
    case Surface::Z_NEG:
    {
      CZP = -kcoeff[2][1]/dist[2][1];
      CZM = 0;
    }
      break;
    case Surface::Z_POS:
    {
      CZP = 0;
      CZM = -kcoeff[2][0]/dist[2][0];
    }
      break;
  }
  Qx = CXP*DTXP + CXM*DTXM;
  Qy = CYP*DTYP + CYM*DTYM;
  Qz = CZP*DTZP + CZM*DTZM;

  Qflux.push_back(Qx);
  Qflux.push_back(Qy);
  Qflux.push_back(Qz);

  return Qflux;
}

void BoundaryCell::zfCellADI(const int &dim, const int &sdim, const int &sign,
                             double &A, double &Alt, double &bVal)
{
  A = 1.0;
  if (dim == sdim) {
    Alt = -1.0;
    bVal = 0;
  } else {
    Alt = 0.0;
    bVal = *(told_ptr + sign*stepsize[sdim]);
  }
}

void BoundaryCell::ifCellADI(const int &dim, const int &sdim, const int &dir,
                             const Foundation &foundation, const BoundaryConditions &bcs,
                             double &A, double &Alt, double &bVal)
{
  double Tair = bcs.indoorTemp;
  double hc = foundation.getConvectionCoeff(*told_ptr,
                                            Tair,0.0,0.00208,false,surfacePtr->tilt);
  double hr = getSimpleInteriorIRCoeff(surfacePtr->emissivity,
                                       *told_ptr,Tair);
  int sign = (dir==0)? -1 : 1;

  A = kcoeff[sdim][dir]/dist[sdim][dir] + (hc + hr);
  if (dim == sdim) {
    Alt = -kcoeff[sdim][dir]/dist[sdim][dir];
    bVal = (hc + hr)*Tair + heatGain;
  } else {
    Alt = 0.0;
    bVal = *(told_ptr + sign*stepsize[sdim])*kcoeff[sdim][1]/dist[sdim][1] + (hc + hr)*Tair + heatGain;
  }
}

void BoundaryCell::efCellADI(const int &dim, const int &sdim, const int &dir,
                             const Foundation &foundation, const BoundaryConditions &bcs,
                             double &A, double &Alt, double &bVal)
{
  double Tair = bcs.outdoorTemp;
  double v = bcs.localWindSpeed;
  double eSky = bcs.skyEmissivity;
  double tilt = surfacePtr->tilt;
  double F = getEffectiveExteriorViewFactor(eSky,tilt);
  double hc = foundation.getConvectionCoeff(*told_ptr,Tair,v,foundation.surfaceRoughness,true,tilt);
  double hr = getExteriorIRCoeff(surfacePtr->emissivity,*told_ptr,Tair,eSky,tilt);
  int sign = (dir==0)? -1 : 1;

  A = kcoeff[sdim][dir]/dist[sdim][dir] + (hc + hr);
  if (dim == sdim) {
    Alt = -kcoeff[sdim][dir]/dist[sdim][dir];
    bVal = (hc + hr*pow(F,0.25))*Tair + heatGain;
  } else {
    Alt = 0.0;
    bVal = *(told_ptr + sign*stepsize[sdim])*kcoeff[sdim][dir]/dist[sdim][dir] + (hc + hr*pow(F,0.25))*Tair + heatGain;
  }
}

void BoundaryCell::zfCellMatrix(double &A, double &Alt, double &bVal)
{
  A = 1.0;
  Alt = -1.0;
  bVal = 0.0;
}

void BoundaryCell::ifCellMatrix(const int &dim, const int &dir,
                                const Foundation &foundation, const BoundaryConditions &bcs,
                                double &A, double &Alt, double &bVal)
{
  double Tair = bcs.indoorTemp;
  double hc = foundation.getConvectionCoeff(*told_ptr,
                                            Tair,0.0,0.00208,false,surfacePtr->tilt);
  double hr = getSimpleInteriorIRCoeff(surfacePtr->emissivity,
                                       *told_ptr,Tair);

  A = kcoeff[dim][dir] / dist[dim][dir] + (hc + hr);
  Alt = -kcoeff[dim][dir] / dist[dim][dir];
  bVal = (hc + hr) * Tair + heatGain;
}

void BoundaryCell::efCellMatrix(const int &dim, const int &dir,
                                const Foundation &foundation, const BoundaryConditions &bcs,
                                double &A, double &Alt, double &bVal)
{
  double Tair = bcs.outdoorTemp;
  double v = bcs.localWindSpeed;
  double eSky = bcs.skyEmissivity;
  double tilt = surfacePtr->tilt;
  double F = getEffectiveExteriorViewFactor(eSky,tilt);
  double hc = foundation.getConvectionCoeff(*told_ptr,Tair,v,foundation.surfaceRoughness,true,tilt);
  double hr = getExteriorIRCoeff(surfacePtr->emissivity,*told_ptr,Tair,eSky,tilt);

  A = kcoeff[dim][dir] / dist[dim][dir] + (hc + hr);
  Alt = -kcoeff[dim][dir] / dist[dim][dir];
  bVal = (hc + hr * pow(F, 0.25)) * Tair + heatGain;
}

ZeroThicknessCell::ZeroThicknessCell(const std::size_t &index, const CellType cellType,
                                 const std::size_t &i, const std::size_t &j, const std::size_t &k,
                                 std::size_t *stepsize,
                                 const Foundation &foundation, Surface *surfacePtr, Block *blockPtr,
                                 Mesher *meshPtr):
        Cell(index, cellType, i, j, k, stepsize, foundation, surfacePtr, blockPtr, meshPtr)
{}

std::vector<double> ZeroThicknessCell::calculateHeatFlux(int ndims, double &TNew,
                                                         std::size_t nX, std::size_t nY, std::size_t nZ,
                                                         const std::vector< std::shared_ptr<Cell> > &cell_v)
{
  std::vector<double> Qflux;
  double Qx = 0;
  double Qy = 0;
  double Qz = 0;

  std::vector<double> Qm;
  std::vector<double> Qp;

  if (isEqual(meshPtr[0].deltas[coords[0]], 0.0))
  {
    Qm = cell_v[index - stepsize[0]]->calculateHeatFlux(ndims, *(&TNew - stepsize[0]), nX, nY, nZ, cell_v);
    Qp = cell_v[index + stepsize[0]]->calculateHeatFlux(ndims, *(&TNew + stepsize[0]), nX, nY, nZ, cell_v);
  }
  if (isEqual(meshPtr[1].deltas[coords[1]], 0.0))
  {
    Qm = cell_v[index - stepsize[1]]->calculateHeatFlux(ndims, *(&TNew - stepsize[1]), nX, nY, nZ, cell_v);
    Qp = cell_v[index + stepsize[1]]->calculateHeatFlux(ndims, *(&TNew + stepsize[1]), nX, nY, nZ, cell_v);
  }
  if (isEqual(meshPtr[2].deltas[coords[2]], 0.0))
  {
    Qm = cell_v[index - stepsize[2]]->calculateHeatFlux(ndims, *(&TNew - stepsize[2]), nX, nY, nZ, cell_v);
    Qp = cell_v[index + stepsize[2]]->calculateHeatFlux(ndims, *(&TNew + stepsize[2]), nX, nY, nZ, cell_v);
  }

  Qx = (Qm[0] + Qp[0])*0.5;
  Qy = (Qm[1] + Qp[1])*0.5;
  Qz = (Qm[2] + Qp[2])*0.5;

  Qflux.push_back(Qx);
  Qflux.push_back(Qy);
  Qflux.push_back(Qz);

  return Qflux;
}


}

#endif
