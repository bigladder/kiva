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

void Cell::setConductivities(const std::vector< std::shared_ptr<Cell> > &cell_v) {
  getKXP(cell_v);
  getKXM(cell_v);
  getKYP(cell_v);
  getKYM(cell_v);
  getKZP(cell_v);
  getKZM(cell_v);
}

void Cell::getKXP(const std::vector< std::shared_ptr<Cell> > &cell_v)
{
  if (coords[0] == meshPtr[0].centers.size() - 1)
  {
    // For boundary cells assume that the cell on the other side of the
    // boundary is the same as the current cell
    kxp = conductivity;
  }
  else
  {
    kxp = 1/(meshPtr[0].deltas[coords[0]]/(2*dist[0][1]*conductivity) +
             meshPtr[0].deltas[coords[0] + 1]/(2*dist[0][1]*cell_v[index + stepsize[0]]->conductivity));
  }
}

void Cell::getKXM(const std::vector< std::shared_ptr<Cell> > &cell_v)
{
  if (coords[0] == 0)
  {
    // For boundary cells assume that the cell on the other side of the
    // boundary is the same as the current cell
    kxm = conductivity;
  }
  else
  {
    kxm = 1/(meshPtr[0].deltas[coords[0]]/(2*dist[0][0]*conductivity) +
             meshPtr[0].deltas[coords[0] - 1]/(2*dist[0][0]*cell_v[index - stepsize[0]]->conductivity));
  }
}

void Cell::getKYP(const std::vector< std::shared_ptr<Cell> > &cell_v)
{
  if (coords[1] == meshPtr[1].centers.size() - 1)
  {
    // For boundary cells assume that the cell on the other side of the
    // boundary is the same as the current cell
    kyp = conductivity;
  }
  else
  {
    kyp = 1/(meshPtr[1].deltas[coords[1]]/(2*dist[1][1]*conductivity) +
             meshPtr[1].deltas[coords[1] + 1]/(2*dist[1][1]*cell_v[index + stepsize[1]]->conductivity));
  }
}

void Cell::getKYM(const std::vector< std::shared_ptr<Cell> > &cell_v)
{
  if (coords[1] == 0)
  {
    // For boundary cells assume that the cell on the other side of the
    // boundary is the same as the current cell
    kym = conductivity;
  }
  else
  {
    kym = 1/(meshPtr[1].deltas[coords[1]]/(2*dist[1][0]*conductivity) +
             meshPtr[1].deltas[coords[1] - 1]/(2*dist[1][0]*cell_v[index - stepsize[1]]->conductivity));
  }
}

void Cell::getKZP(const std::vector< std::shared_ptr<Cell> > &cell_v)
{
  if (coords[2] == meshPtr[2].centers.size() - 1)
  {
    // For boundary cells assume that the cell on the other side of the
    // boundary is the same as the current cell
    kzp = conductivity;
  }
  else
  {
    kzp = 1/(meshPtr[2].deltas[coords[2]]/(2*dist[2][1]*conductivity) +
             meshPtr[2].deltas[coords[2] + 1]/(2*dist[2][1]*cell_v[index + stepsize[2]]->conductivity));
  }
}

void Cell::getKZM(const std::vector< std::shared_ptr<Cell> > &cell_v)
{
  if (coords[2] == 0)
  {
    // For boundary cells assume that the cell on the other side of the
    // boundary is the same as the current cell
    kzm = conductivity;
  }
  else
  {
    kzm = 1/(meshPtr[2].deltas[coords[2]]/(2*dist[2][0]*conductivity) +
             meshPtr[2].deltas[coords[2] - 1]/(2*dist[2][0]*cell_v[index - stepsize[2]]->conductivity));
  }
}

void Cell::setPDEcoefficients(int ndims, bool cylindrical) {
  if (ndims > 1) {
    // Radial X terms
    if (cylindrical)
    {
      cxp_c = (dist[0][0]*kxp) / ((dist[0][0] + dist[0][1])*dist[0][1]);
      cxm_c = (dist[0][1]*kxm) / ((dist[0][0] + dist[0][1])*dist[0][0]);
    }
    else
    {
      cxp_c = 0.0;
      cxm_c = 0.0;
    }

    // Cartesian X terms
    cxp = (2*kxp) /  ((dist[0][0] + dist[0][1])*dist[0][1]);
    cxm = -1*(2*kxm) / ((dist[0][0] + dist[0][1])*dist[0][0]);
  }

  // Cartesian Z terms
  czp = (2*kzp) / ((dist[2][0] + dist[2][1])*dist[2][1]);
  czm = -1*(2*kzm) / ((dist[2][0] + dist[2][1])*dist[2][0]);

  // Cartesian Y terms
  if (ndims == 3)
  {
    cyp = (2*kyp) / ((dist[1][0] + dist[1][1])*dist[1][1]);
    cym = -1*(2*kym) / ((dist[1][0] + dist[1][1])*dist[1][0]);
  }
  else
  {
    cyp = 0.0;
    cym = 0.0;
  }
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

  double CXP = cxp*theta;
  double CXM = cxm*theta;
  double CZP = czp*theta;
  double CZM = czm*theta;
  double CYP = cyp*theta;
  double CYM = cym*theta;
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
      CXPC = cxp_c*theta/r;
      CXMC = cxm_c*theta/r;
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

  double CXP = cxp*theta;
  double CXM = cxm*theta;
  double CZP = czp*theta;
  double CZM = czm*theta;
  double Q = heatGain*theta;

  if (foundation.numberOfDimensions == 3) {
    double CYP = cyp * theta;
    double CYM = cym * theta;
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
      CXPC = cxp_c*theta/r;
      CXMC = cxm_c*theta/r;
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

  double CXP = cxp*theta;
  double CXM = cxm*theta;
  double CZP = czp*theta;
  double CZM = czm*theta;
  double Q = heatGain*theta;

  if (foundation.numberOfDimensions == 3) {
    double CYP = cyp*theta;
    double CYM = cym*theta;
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
      CXPC = cxp_c*theta/r;
      CXMC = cxm_c*theta/r;
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

void Cell::calcCellMatrix(Foundation::NumericalScheme scheme, double timestep, const Foundation &foundation,
                          const BoundaryConditions &/*bcs*/,
                          double &A, double &Aip, double &Aim, double &Ajp, double &Ajm,
                          double &Akp, double &Akm, double &bVal)
{
  if (scheme == Foundation::NS_STEADY_STATE)
  {
    calcCellSteadyState(foundation, A, Aip, Aim, Ajp, Ajm, Akp, Akm, bVal);
  }
  else
  {
    double theta = timestep/
                   (density*specificHeat);

    double f;
    if (scheme == Foundation::NS_IMPLICIT)
      f = 1.0;
    else
      f = 0.5;

    double CXP = cxp*theta;
    double CXM = cxm*theta;
    double CZP = czp*theta;
    double CZM = czm*theta;
    double CYP = cyp*theta;
    double CYM = cym*theta;
    double Q = heatGain*theta;

    if (foundation.numberOfDimensions == 3)
    {
      A = (1.0 + f*(CXP + CZP + CYP - CXM - CZM - CYM));
      Aim = f*CXM;
      Aip = f*(-CXP);
      Akm = f*CZM;
      Akp = f*(-CZP);
      Ajm = f*CYM;
      Ajp = f*(-CYP);

      bVal = *told_ptr*(1.0 + (1-f)*(CXM + CZM + CYM - CXP - CZP - CYP))
             - *(told_ptr - stepsize[0])*(1-f)*CXM
             + *(told_ptr + stepsize[0])*(1-f)*CXP
             - *(told_ptr - stepsize[1])*(1-f)*CZM
             + *(told_ptr + stepsize[1])*(1-f)*CZP
             - *(told_ptr - stepsize[2])*(1-f)*CYM
             + *(told_ptr + stepsize[2])*(1-f)*CYP
             + Q;
    }
    else if (foundation.numberOfDimensions == 2)
    {
      double CXPC = 0;
      double CXMC = 0;

      if (coords[0] != 0)
      {
        CXPC = cxp_c*theta/r;
        CXMC = cxm_c*theta/r;
      }
      A = (1.0 + f*(CXPC + CXP + CZP - CXMC - CXM - CZM));
      Aim = f*(CXMC + CXM);
      Aip = f*(-CXPC - CXP);
      Akm = f*CZM;
      Akp = f*(-CZP);

      bVal = *told_ptr*(1.0 + (1-f)*(CXMC + CXM + CZM - CXPC - CXP - CZP))
             - *(told_ptr - stepsize[0])*(1-f)*(CXMC + CXM)
             + *(told_ptr + stepsize[0])*(1-f)*(CXPC + CXP)
             - *(told_ptr - stepsize[2])*(1-f)*CZM
             + *(told_ptr + stepsize[2])*(1-f)*CZP
             + Q;
    }
    else
    {
      A = (1.0 + f*(CZP - CZM));
      Akm = f*CZM;
      Akp = f*(-CZP);

      bVal = *told_ptr*(1.0 + (1-f)*(CZM - CZP))
             - *(told_ptr - stepsize[2])*(1-f)*CZM
             + *(told_ptr + stepsize[2])*(1-f)*CZP
             + Q;
    }
  }

}

void Cell::calcCellSteadyState(const Foundation &foundation,
                               double &A, double &Aip, double &Aim, double &Ajp, double &Ajm,
                               double &Akp, double &Akm, double &bVal)
{
  double CXP = cxp;
  double CXM = cxm;
  double CZP = czp;
  double CZM = czm;
  double CYP = cyp;
  double CYM = cym;
  double Q = heatGain;

  if (foundation.numberOfDimensions == 3)
  {
    A = (CXM + CZM + CYM - CXP - CZP - CYP);
    Aim = -CXM;
    Aip = CXP;
    Akm = -CZM;
    Akp = CZP;
    Ajm = -CYM;
    Ajp = CYP;

    bVal = -Q;
  }
  else if (foundation.numberOfDimensions == 2)
  {
    double CXPC = 0;
    double CXMC = 0;

    if (coords[0] != 0)
    {
      CXPC = cxp_c/r;
      CXMC = cxm_c/r;
    }
    A = (CXMC + CXM + CZM - CXPC - CXP - CZP);
    Aim = (-CXMC - CXM);
    Aip = (CXPC + CXP);
    Akm = -CZM;
    Akp = CZP;

    bVal = -Q;
  }
  else
  {
    A = (CZM - CZP);
    Akm = -CZM;
    Akp = CZP;

    bVal = -Q;
  }
}

void Cell::calcCellADI(int dim, const Foundation &foundation, double timestep,
                       const BoundaryConditions &/*bcs*/,
                       double &Am, double &A, double &Ap, double &bVal) {
  double theta = timestep / (foundation.numberOfDimensions
                             *density*specificHeat);

  double CXP = cxp*theta;
  double CXM = cxm*theta;
  double CZP = czp*theta;
  double CZM = czm*theta;
  double CYP = cyp*theta;
  double CYM = cym*theta;
  double Q = heatGain*theta;

  double f = foundation.fADI;

  if (foundation.numberOfDimensions == 3)
  {
    if (dim == 1) // x
    {
      A = 1.0 + (3 - 2*f)*(CXP - CXM);
      Am = (3 - 2*f)*CXM;
      Ap = (3 - 2*f)*(-CXP);

      bVal = *told_ptr*(1.0 + f*(CZM + CYM - CZP - CYP))
             - *(told_ptr-stepsize[2])*f*CZM
             + *(told_ptr+stepsize[2])*f*CZP
             - *(told_ptr-stepsize[1])*f*CYM
             + *(told_ptr+stepsize[1])*f*CYP
             + Q;
    }
    else if (dim == 2) // y
    {
      A = (1.0 + (3 - 2*f)*(CYP - CYM));
      Am = (3 - 2*f)*CYM;
      Ap = (3 - 2*f)*(-CYP);

      bVal = *told_ptr*(1.0 + f*(CXM + CZM - CXP - CZP))
             - *(told_ptr-stepsize[0])*f*CXM
             + *(told_ptr+stepsize[0])*f*CXP
             - *(told_ptr-stepsize[2])*f*CZM
             + *(told_ptr+stepsize[2])*f*CZP
             + Q;
    }
    else //if (dim == 3) // z
    {
      A = (1.0 + (3 - 2*f)*(CZP - CZM));
      Am = (3 - 2*f)*CZM;
      Ap = (3 - 2*f)*(-CZP);

      bVal = *told_ptr*(1.0 + f*(CXM + CYM - CXP - CYP))
             - *(told_ptr-stepsize[0])*f*CXM
             + *(told_ptr+stepsize[0])*f*CXP
             - *(told_ptr-stepsize[1])*f*CYM
             + *(told_ptr+stepsize[1])*f*CYP
             + Q;
    }

  }
  else if (foundation.numberOfDimensions == 2)
  {
    double CXPC = 0;
    double CXMC = 0;
    if (coords[0] != 0)
    {
      CXPC = cxp_c*theta/r;
      CXMC = cxm_c*theta/r;
    }
    if (dim == 1) // x
    {
      A = 1.0 + (2 - f)*(CXPC + CXP - CXMC - CXM);
      Am = (2 - f)*(CXMC + CXM);
      Ap = (2 - f)*(-CXPC - CXP);

      bVal = *told_ptr*(1.0 + f*(CZM - CZP))
             - *(told_ptr-stepsize[2])*f*CZM
             + *(told_ptr+stepsize[2])*f*CZP
             + Q;
    }
    else //if (dim == 3) // z
    {
      A = 1.0 + (2 - f)*(CZP - CZM);
      Am = (2 - f)*CZM;
      Ap = (2 - f)*(-CZP);

      bVal = *told_ptr*(1.0 + f*(CXMC + CXM - CXPC - CXP))
             - *(told_ptr-stepsize[0])*f*(CXMC + CXM)
             + *(told_ptr+stepsize[0])*f*(CXPC + CXP)
             + Q;
    }
  }
  else
  {
    A = 1.0 + CZP - CZM;
    Am = CZM;
    Ap = -CZP;

    bVal = *told_ptr + Q;
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
  double CZP = -kzp*dist[2][0]/(dist[2][1]+dist[2][0])/dist[2][1];
  double CZM = -kzm*dist[2][1]/(dist[2][1]+dist[2][0])/dist[2][0];

  if (ndims > 1)
  {
    CXP = -kxp*dist[0][0]/(dist[0][1]+dist[0][0])/dist[0][1];
    CXM = -kxm*dist[0][1]/(dist[0][1]+dist[0][0])/dist[0][0];
  }


  if (ndims == 3)
  {
    CYP = -kyp*dist[1][0]/(dist[1][1]+dist[1][0])/dist[1][1];
    CYM = -kym*dist[1][1]/(dist[1][1]+dist[1][0])/dist[1][0];
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

void ExteriorAirCell::calcCellADI(int /*dim*/, const Foundation &/*foundation*/, double /*timestep*/,
                                  const BoundaryConditions &bcs,
                                  double &/*Am*/, double &A, double &/*Ap*/, double &bVal)
{
  doOutdoorTemp(bcs, A, bVal);
};


void ExteriorAirCell::calcCellMatrix(Foundation::NumericalScheme /*scheme*/, double /*timestep*/, const Foundation &/*foundation*/,
                                     const BoundaryConditions &bcs,
                                     double &A, double &/*Aip*/, double &/*Aim*/, double &/*Ajp*/, double &/*Ajm*/,
                                     double &/*Akp*/, double &/*Akm*/, double &bVal)
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

void InteriorAirCell::calcCellMatrix(Foundation::NumericalScheme /*scheme*/, double /*timestep*/, const Foundation &/*foundation*/,
                                     const BoundaryConditions &bcs,
                                     double &A, double &/*Aip*/, double &/*Aim*/, double &/*Ajp*/, double &/*Ajm*/,
                                     double &/*Akp*/, double &/*Akm*/, double &bVal)
{
  doIndoorTemp(bcs, A, bVal);
}

void InteriorAirCell::calcCellADI(int /*dim*/, const Foundation &/*foundation*/, double /*timestep*/,
                                  const BoundaryConditions &bcs,
                                  double &/*Am*/, double &A, double &/*Ap*/, double &bVal)
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
          U = (kxp * *(told_ptr + stepsize[0]) / dist[0][1] +
                      (hc + hr)*Tair + q)/(kxp/dist[0][1] + (hc + hr));
          break;
        case Surface::X_POS:
          U = (kxm * *(&U - stepsize[0])/dist[0][0] +
                      (hc + hr)*Tair + q)/(kxm/dist[0][0] + (hc + hr));
          break;
        case Surface::Y_NEG:
          U = (kyp * *(told_ptr + stepsize[1]) / dist[1][1] +
                      (hc + hr)*Tair + q)/(kyp/dist[1][1] + (hc + hr));
          break;
        case Surface::Y_POS:
          U = (kym * *(&U - stepsize[1])/dist[1][0] +
                      (hc + hr)*Tair + q)/(kym/dist[1][0] + (hc + hr));
          break;
        case Surface::Z_NEG:
          U = (kzp * *(told_ptr + stepsize[2]) / dist[2][1] +
                      (hc + hr)*Tair + q)/(kzp/dist[2][1] + (hc + hr));
          break;
        case Surface::Z_POS:
          U = (kzm * *(&U - stepsize[2])/dist[2][0] +
                      (hc + hr)*Tair + q)/(kzm/dist[2][0] + (hc + hr));
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
          U = (kxp * *(told_ptr + stepsize[0]) / dist[0][1] +
                      (hc + hr*pow(F,0.25))*Tair + q)/(kxp/dist[0][1] + (hc + hr));
          break;
        case Surface::X_POS:
          U = (kxm * *(&U - stepsize[0])/dist[0][0] +
                      (hc + hr*pow(F,0.25))*Tair + q)/(kxm/dist[0][0] + (hc + hr));
          break;
        case Surface::Y_NEG:
          U = (kyp * *(told_ptr + stepsize[1]) / dist[1][1] +
                      (hc + hr*pow(F,0.25))*Tair + q)/(kyp/dist[1][1] + (hc + hr));
          break;
        case Surface::Y_POS:
          U = (kym * *(&U - stepsize[1])/dist[1][0] +
                      (hc + hr*pow(F,0.25))*Tair + q)/(kym/dist[1][0] + (hc + hr));
          break;
        case Surface::Z_NEG:
          U = (kzp * *(told_ptr + stepsize[2]) / dist[2][1] +
                      (hc + hr*pow(F,0.25))*Tair + q)/(kzp/dist[2][1] + (hc + hr));
          break;
        case Surface::Z_POS:
          U = (kzm * *(&U - stepsize[2])/dist[2][0] +
                      (hc + hr*pow(F,0.25))*Tair + q)/(kzm/dist[2][0] + (hc + hr));
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
          V = (kxp * *(&V + stepsize[0])/dist[0][1] +
                      (hc + hr)*Tair + q)/(kxp/dist[0][1] + (hc + hr));
          break;
        case Surface::X_POS:
          V = (kxm * *(told_ptr - stepsize[0]) / dist[0][0] +
                      (hc + hr)*Tair + q)/(kxm/dist[0][0] + (hc + hr));
          break;
        case Surface::Y_NEG:
          V = (kyp * *(&V + stepsize[1])/dist[1][1] +
                      (hc + hr)*Tair + q)/(kyp/dist[1][1] + (hc + hr));
          break;
        case Surface::Y_POS:
          V = (kym * *(told_ptr - stepsize[1]) / dist[1][0] +
                      (hc + hr)*Tair + q)/(kym/dist[1][0] + (hc + hr));
          break;
        case Surface::Z_NEG:
          V = (kzp * *(&V + stepsize[2])/dist[2][1] +
                      (hc + hr)*Tair + q)/(kzp/dist[2][1] + (hc + hr));
          break;
        case Surface::Z_POS:
          V = (kzm * *(told_ptr - stepsize[2]) / dist[2][0] +
                      (hc + hr)*Tair + q)/(kzm/dist[2][0] + (hc + hr));
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
          V = (kxp * *(&V + stepsize[0])/dist[0][1] +
                      (hc + hr*pow(F,0.25))*Tair + q)/(kxp/dist[0][1] + (hc + hr));
          break;
        case Surface::X_POS:
          V = (kxm * *(told_ptr - stepsize[0]) / dist[0][0] +
                      (hc + hr*pow(F,0.25))*Tair + q)/(kxm/dist[0][0] + (hc + hr));
          break;
        case Surface::Y_NEG:
          V = (kyp * *(&V + stepsize[1])/dist[1][1] +
                      (hc + hr*pow(F,0.25))*Tair + q)/(kyp/dist[1][1] + (hc + hr));
          break;
        case Surface::Y_POS:
          V = (kym * *(told_ptr - stepsize[1]) / dist[1][0] +
                      (hc + hr*pow(F,0.25))*Tair + q)/(kym/dist[1][0] + (hc + hr));
          break;
        case Surface::Z_NEG:
          V = (kzp * *(&V + stepsize[2])/dist[2][1] +
                      (hc + hr*pow(F,0.25))*Tair + q)/(kzp/dist[2][1] + (hc + hr));
          break;
        case Surface::Z_POS:
          V = (kzm * *(told_ptr - stepsize[2]) / dist[2][0] +
                      (hc + hr*pow(F,0.25))*Tair + q)/(kzm/dist[2][0] + (hc + hr));
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
          return (kxp * *(told_ptr + stepsize[0])/dist[0][1] +
                         (hc + hr)*Tair + q)/(kxp/dist[0][1] + (hc + hr));
        case Surface::X_POS:
          return (kxm * *(told_ptr - stepsize[0])/dist[0][0] +
                         (hc + hr)*Tair + q)/(kxm/dist[0][0] + (hc + hr));
        case Surface::Y_NEG:
          return (kyp * *(told_ptr + stepsize[1])/dist[1][1] +
                         (hc + hr)*Tair + q)/(kyp/dist[1][1] + (hc + hr));
        case Surface::Y_POS:
          return (kym * *(told_ptr - stepsize[1])/dist[1][0] +
                         (hc + hr)*Tair + q)/(kym/dist[1][0] + (hc + hr));
        case Surface::Z_NEG:
          return (kzp * *(told_ptr + stepsize[2])/dist[2][1] +
                         (hc + hr)*Tair + q)/(kzp/dist[2][1] + (hc + hr));
        case Surface::Z_POS:
          return (kzm * *(told_ptr - stepsize[2])/dist[2][0] +
                         (hc + hr)*Tair + q)/(kzm/dist[2][0] + (hc + hr));
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
          return (kxp * *(told_ptr + stepsize[0])/dist[0][1] +
                         (hc + hr*pow(F,0.25))*Tair + q)/(kxp/dist[0][1] + (hc + hr));
        case Surface::X_POS:
          return (kxm * *(told_ptr - stepsize[0])/dist[0][0] +
                         (hc + hr*pow(F,0.25))*Tair + q)/(kxm/dist[0][0] + (hc + hr));
        case Surface::Y_NEG:
          return (kyp * *(told_ptr + stepsize[1])/dist[1][1] +
                         (hc + hr*pow(F,0.25))*Tair + q)/(kyp/dist[1][1] + (hc + hr));
        case Surface::Y_POS:
          return (kym * *(told_ptr - stepsize[1])/dist[1][0] +
                         (hc + hr*pow(F,0.25))*Tair + q)/(kym/dist[1][0] + (hc + hr));
        case Surface::Z_NEG:
          return (kzp * *(told_ptr + stepsize[2])/dist[2][1] +
                         (hc + hr*pow(F,0.25))*Tair + q)/(kzp/dist[2][1] + (hc + hr));
        case Surface::Z_POS:
          return (kzm * *(told_ptr - stepsize[2])/dist[2][0] +
                         (hc + hr*pow(F,0.25))*Tair + q)/(kzm/dist[2][0] + (hc + hr));
      }
    }
    break;
  }
}

void BoundaryCell::calcCellADI(int dim, const Foundation &foundation, double /*timestep*/,
                               const BoundaryConditions &bcs,
                               double &Am, double &A, double &Ap, double &bVal)
{
  switch (surfacePtr->boundaryConditionType)
  {
    case Surface::ZERO_FLUX:
    {
      switch (surfacePtr->orientation)
      {
        case Surface::X_NEG:
          A = 1.0;

          if (dim == 1)
          {
            Ap = -1.0;
            bVal = 0;
          }
          else
          {
            Ap = 0.0;
            bVal = *(told_ptr+stepsize[0]);
          }
          break;
        case Surface::X_POS:
          A = 1.0;
          if (dim == 1)
          {
            Am = -1.0;
            bVal = 0;
          }
          else
          {
            Am = 0.0;
            bVal = *(told_ptr-stepsize[0]);
          }
          break;
        case Surface::Y_NEG:
          A = 1.0;
          if (dim == 2)
          {
            Ap = -1.0;
            bVal = 0;
          }
          else
          {
            Ap = 0.0;
            bVal = *(told_ptr+stepsize[1]);
          }
          break;
        case Surface::Y_POS:
          A = 1.0;
          if (dim == 2)
          {
            Am = -1.0;
            bVal = 0;
          }
          else
          {
            Am = 0.0;
            bVal = *(told_ptr-stepsize[1]);
          }
          break;
        case Surface::Z_NEG:
          A = 1.0;
          if (dim == 3)
          {
            Ap = -1.0;
            bVal = 0;
          }
          else
          {
            Ap = 0.0;
            bVal = *(told_ptr+stepsize[2]);
          }
          break;
        case Surface::Z_POS:
          A = 1.0;
          if (dim == 3)
          {
            Am = -1.0;
            bVal = 0;
          }
          else
          {
            Am = 0.0;
            bVal = *(told_ptr-stepsize[2]);
          }
          break;
      }
    }
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
          A = kxp/dist[0][1] + (hc + hr);
          if (dim == 1)
          {
            Ap = -kxp/dist[0][1];
            bVal = (hc + hr)*Tair + q;
          }
          else
          {
            Ap = 0.0;
            bVal = *(told_ptr+stepsize[0])*kxp/dist[0][1] + (hc + hr)*Tair + q;
          }
          break;
        case Surface::X_POS:
          A = kxm/dist[0][0] + (hc + hr);
          if (dim == 1)
          {
            Am = -kxm/dist[0][0];
            bVal = (hc + hr)*Tair + q;
          }
          else
          {
            Am = 0.0;
            bVal = *(told_ptr-stepsize[0])*kxm/dist[0][0] + (hc + hr)*Tair + q;
          }
          break;
        case Surface::Y_NEG:
          A = kyp/dist[1][1] + (hc + hr);
          if (dim == 2)
          {
            Ap = -kyp/dist[1][1];
            bVal = (hc + hr)*Tair + q;
          }
          else
          {
            Ap = 0.0;
            bVal = *(told_ptr+stepsize[1])*kyp/dist[1][1] + (hc + hr)*Tair + q;
          }
          break;
        case Surface::Y_POS:
          A = kym/dist[1][0] + (hc + hr);
          if (dim == 2)
          {
            Am = -kym/dist[1][0];
            bVal = (hc + hr)*Tair + q;
          }
          else
          {
            Am = 0.0;
            bVal = *(told_ptr-stepsize[1])*kym/dist[1][0] + (hc + hr)*Tair + q;
          }
          break;
        case Surface::Z_NEG:
          A = kzp/dist[2][1] + (hc + hr);
          if (dim == 3)
          {
            Ap = -kzp/dist[2][1];
            bVal = (hc + hr)*Tair + q;
          }
          else
          {
            Ap = 0.0;
            bVal = *(told_ptr+stepsize[2])*kzp/dist[2][1] + (hc + hr)*Tair + q;
          }
          break;
        case Surface::Z_POS:
          A = kzm/dist[2][0] + (hc + hr);
          if (dim == 3)
          {
            Am = -kzm/dist[2][0];
            bVal = (hc + hr)*Tair + q;
          }
          else
          {
            Am = 0.0;
            bVal = *(told_ptr-stepsize[2])*kzm/dist[2][0] + (hc + hr)*Tair + q;
          }
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
          A = kxp/dist[0][1] + (hc + hr);
          if (dim == 1)
          {
            Ap = -kxp/dist[0][1];
            bVal = (hc + hr*pow(F,0.25))*Tair + q;
          }
          else
          {
            Ap = 0.0;
            bVal = *(told_ptr+stepsize[0])*kxp/dist[0][1] + (hc + hr*pow(F,0.25))*Tair + q;
          }
          break;
        case Surface::X_POS:
          A = kxm/dist[0][0] + (hc + hr);
          if (dim == 1)
          {
            Am = -kxm/dist[0][0];
            bVal = (hc + hr*pow(F,0.25))*Tair + q;
          }
          else
          {
            Am = 0.0;
            bVal = *(told_ptr-stepsize[0])*kxm/dist[0][0] + (hc + hr*pow(F,0.25))*Tair + q;
          }
          break;
        case Surface::Y_NEG:
          A = kyp/dist[1][1] + (hc + hr);
          if (dim == 2)
          {
            Ap = -kyp/dist[1][1];
            bVal = (hc + hr*pow(F,0.25))*Tair + q;
          }
          else
          {
            Ap = 0.0;
            bVal = *(told_ptr+stepsize[1])*kyp/dist[1][1] + (hc + hr*pow(F,0.25))*Tair + q;
          }
          break;
        case Surface::Y_POS:
          A = kym/dist[1][0] + (hc + hr);
          if (dim == 2)
          {
            Am = -kym/dist[1][0];
            bVal = (hc + hr*pow(F,0.25))*Tair + q;
          }
          else
          {
            Am = 0.0;
            bVal = *(told_ptr-stepsize[1])*kym/dist[1][0] + (hc + hr*pow(F,0.25))*Tair + q;
          }
          break;
        case Surface::Z_NEG:
          A = kzp/dist[2][1] + (hc + hr);
          if (dim == 3)
          {
            Ap = -kzp/dist[2][1];
            bVal = (hc + hr*pow(F,0.25))*Tair + q;
          }
          else
          {
            Ap = 0.0;
            bVal = *(told_ptr+stepsize[2])*kzp/dist[2][1] + (hc + hr*pow(F,0.25))*Tair + q;
          }
          break;
        case Surface::Z_POS:
          A = kzm/dist[2][0] + (hc + hr);
          if (dim == 3)
          {
            Am = -kzm/dist[2][0];
            bVal = (hc + hr*pow(F,0.25))*Tair + q;
          }
          else
          {
            Am = 0.0;
            bVal = *(told_ptr-stepsize[2])*kzm/dist[2][0] + (hc + hr*pow(F,0.25))*Tair + q;
          }
          break;
      }
    }
      break;
  }

}

void BoundaryCell::calcCellMatrix(Kiva::Foundation::NumericalScheme /*scheme*/, double /*timestep*/,
                                  const Kiva::Foundation &foundation, const Kiva::BoundaryConditions &bcs,
                                  double &A, double &Aip, double &Aim, double &Ajp, double &Ajm,
                                  double &Akp, double &Akm, double &bVal)
{
  switch (surfacePtr->boundaryConditionType)
  {
    case Surface::ZERO_FLUX:
    {
      switch (surfacePtr->orientation)
      {
        case Surface::X_NEG: {
          A = 1.0;
          Aip = -1.0;
          bVal = 0.0;
          break;
        }
        case Surface::X_POS: {
          A = 1.0;
          Aim = -1.0;
          bVal = 0.0;
          break;
        }
        case Surface::Y_NEG: {
          A = 1.0;
          Ajp = -1.0;
          bVal = 0.0;
          break;
        }
        case Surface::Y_POS: {
          A = 1.0;
          Ajm = -1.0;
          bVal = 0.0;
          break;
        }
        case Surface::Z_NEG: {
          A = 1.0;
          Akp = -1.0;
          bVal = 0.0;
          break;
        }
        case Surface::Z_POS: {
          A = 1.0;
          Akm = -1.0;
          bVal = 0.0;
          break;
        }
      }
    }
      break;
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
        case Surface::X_NEG: {
          A = kxp / dist[0][1] + (hc + hr);
          Aip = -kxp / dist[0][1];
          bVal = (hc + hr) * Tair + q;
          break;
        }
        case Surface::X_POS: {
          A = kxm / dist[0][0] + (hc + hr);
          Aim = -kxm / dist[0][0];
          bVal = (hc + hr) * Tair + q;
          break;
        }
        case Surface::Y_NEG: {
          A = kyp / dist[1][1] + (hc + hr);
          Ajp = -kyp / dist[1][1];
          bVal = (hc + hr) * Tair + q;
          break;
        }
        case Surface::Y_POS: {
          A = kym / dist[1][0] + (hc + hr);
          Ajm = -kym / dist[1][0];
          bVal = (hc + hr) * Tair + q;
          break;
        }
        case Surface::Z_NEG: {
          A = kzp / dist[2][1] + (hc + hr);
          Akp = -kzp / dist[2][1];
          bVal = (hc + hr) * Tair + q;
          break;
        }
        case Surface::Z_POS: {
          A = kzm / dist[2][0] + (hc + hr);
          Akm = -kzm / dist[2][0];
          bVal = (hc + hr) * Tair + q;
          break;
        }
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
        case Surface::X_NEG: {
          A = kxp / dist[0][1] + (hc + hr);
          Aip = -kxp / dist[0][1];
          bVal = (hc + hr * pow(F, 0.25)) * Tair + q;
          break;
        }
        case Surface::X_POS: {
          A = kxm / dist[0][0] + (hc + hr);
          Aim = -kxm / dist[0][0];
          bVal = (hc + hr * pow(F, 0.25)) * Tair + q;
          break;
        }
        case Surface::Y_NEG: {
          A = kyp / dist[1][1] + (hc + hr);
          Ajp = -kyp / dist[1][1];
          bVal = (hc + hr * pow(F, 0.25)) * Tair + q;
          break;
        }
        case Surface::Y_POS: {
          A = kym / dist[1][0] + (hc + hr);
          Ajm = -kym / dist[1][0];
          bVal = (hc + hr * pow(F, 0.25)) * Tair + q;
          break;
        }
        case Surface::Z_NEG: {
          A = kzp / dist[2][1] + (hc + hr);
          Akp = -kzp / dist[2][1];
          bVal = (hc + hr * pow(F, 0.25)) * Tair + q;
          break;
        }
        case Surface::Z_POS: {
          A = kzm / dist[2][0] + (hc + hr);
          Akm = -kzm / dist[2][0];
          bVal = (hc + hr * pow(F, 0.25)) * Tair + q;
          break;
        }
      }
    }
      break;
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
  double CZP = -kzp*dist[2][0]/(dist[2][1]+dist[2][0])/dist[2][1];
  double CZM = -kzm*dist[2][1]/(dist[2][1]+dist[2][0])/dist[2][0];

  if (ndims > 1)
  {
    CXP = -kxp*dist[0][0]/(dist[0][1]+dist[0][0])/dist[0][1];
    CXM = -kxm*dist[0][1]/(dist[0][1]+dist[0][0])/dist[0][0];
  }


  if (ndims == 3)
  {
    CYP = -kyp*dist[1][0]/(dist[1][1]+dist[1][0])/dist[1][1];
    CYM = -kym*dist[1][1]/(dist[1][1]+dist[1][0])/dist[1][0];
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
      CXP = -kxp/dist[0][1];
      CXM = 0;
    }
      break;
    case Surface::X_POS:
    {
      CXP = 0;
      CXM = -kxm/dist[0][0];
    }
      break;
    case Surface::Y_NEG:
    {
      CYP = -kyp/dist[1][1];
      CYM = 0;
    }
      break;
    case Surface::Y_POS:
    {
      CYP = 0;
      CYM = -kym/dist[1][0];
    }
      break;
    case Surface::Z_NEG:
    {
      CZP = -kzp/dist[2][1];
      CZM = 0;
    }
      break;
    case Surface::Z_POS:
    {
      CZP = 0;
      CZM = -kzm/dist[2][0];
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
