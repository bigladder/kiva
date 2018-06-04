/* Copyright (c) 2012-2017 Big Ladder Software LLC. All rights reserved.
* See the LICENSE file for additional terms and conditions. */

#ifndef Ground_CPP
#define Ground_CPP

#undef PRNTSURF

#include "Ground.hpp"
#include "Errors.hpp"
//#include <unsupported/Eigen/SparseExtra>

namespace Kiva {

static const double PI = 4.0*atan(1.0);

static const bool TDMA = true;

Ground::Ground(Foundation &foundation) : foundation(foundation)
{
  pSolver = std::make_shared<Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>>>();
}

Ground::Ground(Foundation &foundation, GroundOutput::OutputMap &outputMap)
  : foundation(foundation), groundOutput(outputMap)
{
  pSolver = std::make_shared<Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>>>();
}

Ground::~Ground()
{
}

void Ground::buildDomain()
{
  // Create mesh
  foundation.createMeshData();

  // Build matrices for PDE term coefficients
  domain.setDomain(foundation);

  nX = domain.meshX.centers.size();
  nY = domain.meshY.centers.size();
  nZ = domain.meshZ.centers.size();
  num_cells = nX*nY*nZ;

  // Initialize matices
  if (foundation.numericalScheme == Foundation::NS_ADE)
  {
    U.resize(num_cells);
    UOld.resize(num_cells);

    V.resize(num_cells);
    VOld.resize(num_cells);
  }

  if ((foundation.numericalScheme == Foundation::NS_ADI ||
    foundation.numberOfDimensions == 1) && TDMA)
  {
    a1.resize(num_cells, 0.0);
    a2.resize(num_cells, 0.0);
    a3.resize(num_cells, 0.0);
    b_.resize(num_cells, 0.0);
    x_.resize(num_cells);
  }

  pSolver->setMaxIterations(foundation.maxIterations);
  pSolver->setTolerance(foundation.tolerance);
  tripletList.reserve(num_cells*(1+2*foundation.numberOfDimensions));
  Amat.resize(num_cells, num_cells);
  b.resize(num_cells);
  x.resize(num_cells);
  x.fill(283.15);

  TNew.resize(num_cells);
  TOld.resize(num_cells);

  link_cells_to_temp();
}

void Ground::calculateADE()
{
  // Set Old values
  for (size_t index = 0; index < num_cells; ++index)
  {
        UOld[index] = VOld[index] = TOld[index];
  }

  // Solve for new values (Main loop)
  #pragma omp parallel sections num_threads(2)
  {
    #pragma omp section
      calculateADEUpwardSweep();
    #pragma omp section
      calculateADEDownwardSweep();
  }
  for (size_t index = 0; index < num_cells; ++index)
  {
        TNew[index] = 0.5*(U[index] + V[index]);

        // Update old values for next timestep
        TOld[index] = TNew[index];
  }
}

void Ground::calculateADEUpwardSweep()
{
  // Upward sweep (Solve U Matrix starting from 1, 1)
  for (size_t index = 0; index < num_cells; index++)
  {
        Cell* this_cell = &domain.cell[index];
        switch (this_cell->cellType)
        {
        case Cell::BOUNDARY:
          {

          switch (this_cell->surfacePtr->boundaryConditionType)
          {
          case Surface::ZERO_FLUX:
            {
            switch (this_cell->surfacePtr->orientation)
            {
            case Surface::X_NEG:
              U[index] = UOld[index+domain.stepsize_i];
              break;
            case Surface::X_POS:
              U[index] = U[index-domain.stepsize_i];
              break;
            case Surface::Y_NEG:
              U[index] = UOld[index+domain.stepsize_j];
              break;
            case Surface::Y_POS:
              U[index] = U[index-domain.stepsize_j];
              break;
            case Surface::Z_NEG:
              U[index] = UOld[index+domain.stepsize_k];
              break;
            case Surface::Z_POS:
              U[index] = U[index-domain.stepsize_k];
              break;
            }
            }
            break;

          case Surface::CONSTANT_TEMPERATURE:

            U[index] = this_cell->surfacePtr->temperature;
            break;

          case Surface::INTERIOR_TEMPERATURE:

            U[index] = bcs.indoorTemp;
            break;

          case Surface::EXTERIOR_TEMPERATURE:

            U[index] = bcs.outdoorTemp;
            break;

          case Surface::INTERIOR_FLUX:
            {
            double& Tair = bcs.indoorTemp;
            double& q = this_cell->heatGain;

            double hc = getConvectionCoeff(*this_cell->told_ptr,
                    Tair,0.0,0.00208,false,this_cell->surfacePtr->tilt);  // TODO Make roughness a property of the interior surfaces
            double hr = getSimpleInteriorIRCoeff(this_cell->surfacePtr->emissivity,
                               *this_cell->told_ptr,Tair);

            switch (this_cell->surfacePtr->orientation)
            {
            case Surface::X_NEG:
              U[index] = (this_cell->kxp*UOld[index+domain.stepsize_i]/this_cell->dxp +
                   (hc + hr)*Tair + q)/(this_cell->kxp/this_cell->dxp + (hc + hr));
              break;
            case Surface::X_POS:
              U[index] = (this_cell->kxm*U[index-domain.stepsize_i]/this_cell->dxm +
                   (hc + hr)*Tair + q)/(this_cell->kxm/this_cell->dxm + (hc + hr));
              break;
            case Surface::Y_NEG:
              U[index] = (this_cell->kyp*UOld[index+domain.stepsize_j]/this_cell->dyp +
                   (hc + hr)*Tair + q)/(this_cell->kyp/this_cell->dyp + (hc + hr));
              break;
            case Surface::Y_POS:
              U[index] = (this_cell->kym*U[index-domain.stepsize_j]/this_cell->dym +
                   (hc + hr)*Tair + q)/(this_cell->kym/this_cell->dym + (hc + hr));
              break;
            case Surface::Z_NEG:
              U[index] = (this_cell->kzp*UOld[index+domain.stepsize_k]/this_cell->dzp +
                   (hc + hr)*Tair + q)/(this_cell->kzp/this_cell->dzp + (hc + hr));
              break;
            case Surface::Z_POS:
              U[index] = (this_cell->kzm*U[index-domain.stepsize_k]/this_cell->dzm +
                   (hc + hr)*Tair + q)/(this_cell->kzm/this_cell->dzm + (hc + hr));
              break;
            }
            }
            break;

          case Surface::EXTERIOR_FLUX:
            {
            double Tair = bcs.outdoorTemp;
            double v = bcs.localWindSpeed;
            double eSky = bcs.skyEmissivity;
            double tilt = this_cell->surfacePtr->tilt;
            double F = getEffectiveExteriorViewFactor(eSky,tilt);
            double hc = getConvectionCoeff(*this_cell->told_ptr,Tair,v,foundation.surfaceRoughness,true,tilt);
            double hr = getExteriorIRCoeff(this_cell->surfacePtr->emissivity,*this_cell->told_ptr,Tair,eSky,tilt);
            double q = this_cell->heatGain;

            switch (this_cell->surfacePtr->orientation)
            {
            case Surface::X_NEG:
              U[index] = (this_cell->kxp*UOld[index+domain.stepsize_i]/this_cell->dxp +
                    (hc + hr*pow(F,0.25))*Tair + q)/(this_cell->kxp/this_cell->dxp + (hc + hr));
              break;
            case Surface::X_POS:
              U[index] = (this_cell->kxm*U[index-domain.stepsize_i]/this_cell->dxm +
                    (hc + hr*pow(F,0.25))*Tair + q)/(this_cell->kxm/this_cell->dxm + (hc + hr));
              break;
            case Surface::Y_NEG:
              U[index] = (this_cell->kyp*UOld[index+domain.stepsize_j]/this_cell->dyp +
                    (hc + hr*pow(F,0.25))*Tair + q)/(this_cell->kyp/this_cell->dyp + (hc + hr));
              break;
            case Surface::Y_POS:
              U[index] = (this_cell->kym*U[index-domain.stepsize_j]/this_cell->dym +
                    (hc + hr*pow(F,0.25))*Tair + q)/(this_cell->kym/this_cell->dym + (hc + hr));
              break;
            case Surface::Z_NEG:
              U[index] = (this_cell->kzp*UOld[index+domain.stepsize_k]/this_cell->dzp +
                    (hc + hr*pow(F,0.25))*Tair + q)/(this_cell->kzp/this_cell->dzp + (hc + hr));
              break;
            case Surface::Z_POS:
              U[index] = (this_cell->kzm*U[index-domain.stepsize_k]/this_cell->dzm +
                    (hc + hr*pow(F,0.25))*Tair + q)/(this_cell->kzm/this_cell->dzm + (hc + hr));
              break;
            }
            }
            break;
          }
          }
          break;
        case Cell::INTERIOR_AIR:
          U[index] = bcs.indoorTemp;
          break;
        case Cell::EXTERIOR_AIR:
          U[index] = bcs.outdoorTemp;
          break;
        default:
          {
          double theta = timestep/
            (this_cell->density*this_cell->specificHeat);

          double CXP = this_cell->cxp*theta;
          double CXM = this_cell->cxm*theta;
          double CZP = this_cell->czp*theta;
          double CZM = this_cell->czm*theta;
          double CYP = this_cell->cyp*theta;
          double CYM = this_cell->cym*theta;
          double Q = this_cell->heatGain*theta;

          if (foundation.numberOfDimensions == 3)
            U[index] = (UOld[index]*(1.0 - CXP - CZP - CYP)
                - U[index-domain.stepsize_i]*CXM
                + UOld[index+domain.stepsize_i]*CXP
                - U[index-domain.stepsize_k]*CZM
                + UOld[index+domain.stepsize_k]*CZP
                - U[index-domain.stepsize_j]*CYM
                + UOld[index+domain.stepsize_j]*CYP
                + Q) /
                (1.0 - CXM - CZM - CYM);
          else if (foundation.numberOfDimensions == 2)
          {
            double CXPC = 0;
            double CXMC = 0;

            if (this_cell->i != 0)
            {
              CXPC = this_cell->cxp_c*theta/this_cell->r;
              CXMC = this_cell->cxm_c*theta/this_cell->r;
            }
            U[index] = (UOld[index]*(1.0 - CXPC - CXP - CZP)
                - U[index-domain.stepsize_i]*(CXMC + CXM)
                + UOld[index+domain.stepsize_i]*(CXPC + CXP)
                - U[index-domain.stepsize_k]*CZM
                + UOld[index+domain.stepsize_k]*CZP
                + Q) /
                (1.0 - CXMC - CXM - CZM);
          }
          else
          {
            U[index] = (UOld[index]*(1.0 - CZP)
                - U[index-domain.stepsize_k]*CZM
                + UOld[index+domain.stepsize_k]*CZP
                + Q) /
                (1.0 - CZM);
          }
          }
          break;
        }
  }
}

void Ground::calculateADEDownwardSweep()
{
  // Downward sweep (Solve V Matrix starting from I, K)
  for (size_t index = num_cells - 1; /* i >= 0 && */ index < num_cells; index--)
  {
        Cell* this_cell = &domain.cell[index];
        switch (this_cell->cellType)
        {
        case Cell::BOUNDARY:
          {

          switch (this_cell->surfacePtr->boundaryConditionType)
          {
          case Surface::ZERO_FLUX:
            {
            switch (this_cell->surfacePtr->orientation)
            {
            case Surface::X_NEG:
              V[index] = V[index+domain.stepsize_i];
              break;
            case Surface::X_POS:
              V[index] = VOld[index-domain.stepsize_i];
              break;
            case Surface::Y_NEG:
              V[index] = V[index+domain.stepsize_j];
              break;
            case Surface::Y_POS:
              V[index] = VOld[index-domain.stepsize_j];
              break;
            case Surface::Z_NEG:
              V[index] = V[index+domain.stepsize_k];
              break;
            case Surface::Z_POS:
              V[index] = VOld[index-domain.stepsize_k];
              break;
            }
            }
            break;

          case Surface::CONSTANT_TEMPERATURE:

            V[index] = this_cell->surfacePtr->temperature;
            break;

          case Surface::INTERIOR_TEMPERATURE:

            V[index] = bcs.indoorTemp;
            break;

          case Surface::EXTERIOR_TEMPERATURE:

            V[index] = bcs.outdoorTemp;
            break;

          case Surface::INTERIOR_FLUX:
            {
            double& Tair = bcs.indoorTemp;
            double& q = this_cell->heatGain;

            double hc = getConvectionCoeff(*this_cell->told_ptr,
                    Tair,0.0,0.00208,false,this_cell->surfacePtr->tilt);
            double hr = getSimpleInteriorIRCoeff(this_cell->surfacePtr->emissivity,
                               *this_cell->told_ptr,Tair);

            switch (this_cell->surfacePtr->orientation)
            {
            case Surface::X_NEG:
              V[index] = (this_cell->kxp*V[index+domain.stepsize_i]/this_cell->dxp +
                    (hc + hr)*Tair + q)/(this_cell->kxp/this_cell->dxp + (hc + hr));
              break;
            case Surface::X_POS:
              V[index] = (this_cell->kxm*VOld[index-domain.stepsize_i]/this_cell->dxm +
                    (hc + hr)*Tair + q)/(this_cell->kxm/this_cell->dxm + (hc + hr));
              break;
            case Surface::Y_NEG:
              V[index] = (this_cell->kyp*V[index+domain.stepsize_j]/this_cell->dyp +
                    (hc + hr)*Tair + q)/(this_cell->kyp/this_cell->dyp + (hc + hr));
              break;
            case Surface::Y_POS:
              V[index] = (this_cell->kym*VOld[index-domain.stepsize_j]/this_cell->dym +
                    (hc + hr)*Tair + q)/(this_cell->kym/this_cell->dym + (hc + hr));
              break;
            case Surface::Z_NEG:
              V[index] = (this_cell->kzp*V[index+domain.stepsize_k]/this_cell->dzp +
                    (hc + hr)*Tair + q)/(this_cell->kzp/this_cell->dzp + (hc + hr));
              break;
            case Surface::Z_POS:
              V[index] = (this_cell->kzm*VOld[index-domain.stepsize_k]/this_cell->dzm +
                    (hc + hr)*Tair + q)/(this_cell->kzm/this_cell->dzm + (hc + hr));
              break;
            }
            }
            break;

          case Surface::EXTERIOR_FLUX:
            {
            double& Tair = bcs.outdoorTemp;
            double& v = bcs.localWindSpeed;
            double& eSky = bcs.skyEmissivity;
            double tilt = this_cell->surfacePtr->tilt;
            double F = getEffectiveExteriorViewFactor(eSky,tilt);
            double hc = getConvectionCoeff(*this_cell->told_ptr,Tair,v,foundation.surfaceRoughness,true,tilt);
            double hr = getExteriorIRCoeff(this_cell->surfacePtr->emissivity,*this_cell->told_ptr,Tair,eSky,tilt);
            double q = this_cell->heatGain;

            switch (this_cell->surfacePtr->orientation)
            {
            case Surface::X_NEG:
              V[index] = (this_cell->kxp*V[index+domain.stepsize_i]/this_cell->dxp +
                  (hc + hr*pow(F,0.25))*Tair + q)/(this_cell->kxp/this_cell->dxp + (hc + hr));
              break;
            case Surface::X_POS:
              V[index] = (this_cell->kxm*VOld[index-domain.stepsize_i]/this_cell->dxm +
                  (hc + hr*pow(F,0.25))*Tair + q)/(this_cell->kxm/this_cell->dxm + (hc + hr));
              break;
            case Surface::Y_NEG:
              V[index] = (this_cell->kyp*V[index+domain.stepsize_j]/this_cell->dyp +
                  (hc + hr*pow(F,0.25))*Tair + q)/(this_cell->kyp/this_cell->dyp + (hc + hr));
              break;
            case Surface::Y_POS:
              V[index] = (this_cell->kym*VOld[index-domain.stepsize_j]/this_cell->dym +
                  (hc + hr*pow(F,0.25))*Tair + q)/(this_cell->kym/this_cell->dym + (hc + hr));
              break;
            case Surface::Z_NEG:
              V[index] = (this_cell->kzp*VOld[index+domain.stepsize_k]/this_cell->dzp +
                  (hc + hr*pow(F,0.25))*Tair + q)/(this_cell->kzp/this_cell->dzp + (hc + hr));
              break;
            case Surface::Z_POS:
              V[index] = (this_cell->kzm*VOld[index-domain.stepsize_k]/this_cell->dzm +
                  (hc + hr*pow(F,0.25))*Tair + q)/(this_cell->kzm/this_cell->dzm + (hc + hr));
              break;
            }
            }
            break;
          }
          }
          break;

        case Cell::INTERIOR_AIR:
          V[index] = bcs.indoorTemp;
          break;
        case Cell::EXTERIOR_AIR:
          V[index] = bcs.outdoorTemp;
          break;
        default:
          {
          double theta = timestep/
            (this_cell->density*this_cell->specificHeat);

          double CXP = this_cell->cxp*theta;
          double CXM = this_cell->cxm*theta;
          double CZP = this_cell->czp*theta;
          double CZM = this_cell->czm*theta;
          double CYP = this_cell->cyp*theta;
          double CYM = this_cell->cym*theta;
          double Q = this_cell->heatGain*theta;

          if (foundation.numberOfDimensions == 3)
            V[index] = (VOld[index]*(1.0 + CXM + CZM + CYM)
                - VOld[index-domain.stepsize_i]*CXM
                + V[index+domain.stepsize_i]*CXP
                - VOld[index-domain.stepsize_k]*CZM
                + V[index+domain.stepsize_k]*CZP
                - VOld[index-domain.stepsize_j]*CYM
                + V[index+domain.stepsize_j]*CYP
                + Q) /
                (1.0 + CXP + CZP + CYP);
          else if (foundation.numberOfDimensions == 2)
          {
            double CXPC = 0;
            double CXMC = 0;

            if (this_cell->i != 0)
            {
              CXPC = this_cell->cxp_c*theta/this_cell->r;
              CXMC = this_cell->cxm_c*theta/this_cell->r;
            }
            V[index] = (VOld[index]*(1.0 + CXMC + CXM + CZM)
                - VOld[index-domain.stepsize_i]*(CXMC + CXM)
                + V[index+domain.stepsize_i]*(CXPC + CXP)
                - VOld[index-domain.stepsize_k]*CZM
                + V[index+domain.stepsize_k]*CZP
                + Q) /
                (1.0 + CXPC + CXP + CZP);
          }
          else
          {
            V[index] = (VOld[index]*(1.0 + CZM)
                - VOld[index-domain.stepsize_k]*CZM
                + V[index+domain.stepsize_k]*CZP
                + Q) /
                (1.0 + CZP);
          }
          }
          break;
        }
  }
}

void Ground::calculateExplicit()
{
  for (size_t index = 0; index < num_cells; index++)
  {
        Cell* this_cell = &domain.cell[index];
        switch (this_cell->cellType)
        {
        case Cell::BOUNDARY:
          {

          switch (this_cell->surfacePtr->boundaryConditionType)
          {
          case Surface::ZERO_FLUX:
            {
            switch (this_cell->surfacePtr->orientation)
            {
            case Surface::X_NEG:
              TNew[index] = *this_cell->i_up_Ptr->told_ptr;
              break;
            case Surface::X_POS:
              TNew[index] = *this_cell->i_down_Ptr->told_ptr;
              break;
            case Surface::Y_NEG:
              TNew[index] = *this_cell->j_up_Ptr->told_ptr;
              break;
            case Surface::Y_POS:
              TNew[index] = *this_cell->j_down_Ptr->told_ptr;
              break;
            case Surface::Z_NEG:
              TNew[index] = *this_cell->k_up_Ptr->told_ptr;
              break;
            case Surface::Z_POS:
              TNew[index] = *this_cell->k_down_Ptr->told_ptr;
              break;
            }
            }
            break;

          case Surface::CONSTANT_TEMPERATURE:

            TNew[index] = this_cell->surfacePtr->temperature;
            break;

          case Surface::INTERIOR_TEMPERATURE:

            TNew[index] = bcs.indoorTemp;
            break;

          case Surface::EXTERIOR_TEMPERATURE:

            TNew[index] = bcs.outdoorTemp;
            break;

          case Surface::INTERIOR_FLUX:
            {
            double& Tair = bcs.indoorTemp;
            double& q = this_cell->heatGain;

            double hc = getConvectionCoeff(*this_cell->told_ptr,
                    Tair,0.0,0.00208,false,this_cell->surfacePtr->tilt);
            double hr = getSimpleInteriorIRCoeff(this_cell->surfacePtr->emissivity,
                               *this_cell->told_ptr,Tair);

            switch (this_cell->surfacePtr->orientation)
            {
            case Surface::X_NEG:
              TNew[index] = (this_cell->kxp * *this_cell->i_up_Ptr->told_ptr/this_cell->dxp +
                    (hc + hr)*Tair + q)/(this_cell->kxp/this_cell->dxp + (hc + hr));
              break;
            case Surface::X_POS:
              TNew[index] = (this_cell->kxm * *this_cell->i_down_Ptr->told_ptr/this_cell->dxm +
                    (hc + hr)*Tair + q)/(this_cell->kxm/this_cell->dxm + (hc + hr));
              break;
            case Surface::Y_NEG:
              TNew[index] = (this_cell->kyp * *this_cell->j_up_Ptr->told_ptr/this_cell->dyp +
                    (hc + hr)*Tair + q)/(this_cell->kyp/this_cell->dyp + (hc + hr));
              break;
            case Surface::Y_POS:
              TNew[index] = (this_cell->kym * *this_cell->j_down_Ptr->told_ptr/this_cell->dym +
                    (hc + hr)*Tair + q)/(this_cell->kym/this_cell->dym + (hc + hr));
              break;
            case Surface::Z_NEG:
              TNew[index] = (this_cell->kzp * *this_cell->k_up_Ptr->told_ptr/this_cell->dzp +
                    (hc + hr)*Tair + q)/(this_cell->kzp/this_cell->dzp + (hc + hr));
              break;
            case Surface::Z_POS:
              TNew[index] = (this_cell->kzm * *this_cell->k_down_Ptr->told_ptr/this_cell->dzm +
                    (hc + hr)*Tair + q)/(this_cell->kzm/this_cell->dzm + (hc + hr));
              break;
            }
            }
            break;

          case Surface::EXTERIOR_FLUX:
            {
            double& Tair = bcs.outdoorTemp;
            double& v = bcs.localWindSpeed;
            double& eSky = bcs.skyEmissivity;
            double tilt = this_cell->surfacePtr->tilt;
            double F = getEffectiveExteriorViewFactor(eSky,tilt);
            double hc = getConvectionCoeff(*this_cell->told_ptr,Tair,v,foundation.surfaceRoughness,true,tilt);
            double hr = getExteriorIRCoeff(this_cell->surfacePtr->emissivity,*this_cell->told_ptr,Tair,eSky,tilt);
            double q = this_cell->heatGain;

            switch (this_cell->surfacePtr->orientation)
            {
            case Surface::X_NEG:
              TNew[index] = (this_cell->kxp * *this_cell->i_up_Ptr->told_ptr/this_cell->dxp +
                  (hc + hr*pow(F,0.25))*Tair + q)/(this_cell->kxp/this_cell->dxp + (hc + hr));
              break;
            case Surface::X_POS:
              TNew[index] = (this_cell->kxm * *this_cell->i_down_Ptr->told_ptr/this_cell->dxm +
                  (hc + hr*pow(F,0.25))*Tair + q)/(this_cell->kxm/this_cell->dxm + (hc + hr));
              break;
            case Surface::Y_NEG:
              TNew[index] = (this_cell->kyp * *this_cell->j_up_Ptr->told_ptr/this_cell->dyp +
                  (hc + hr*pow(F,0.25))*Tair + q)/(this_cell->kyp/this_cell->dyp + (hc + hr));
              break;
            case Surface::Y_POS:
              TNew[index] = (this_cell->kym * *this_cell->j_down_Ptr->told_ptr/this_cell->dym +
                  (hc + hr*pow(F,0.25))*Tair + q)/(this_cell->kym/this_cell->dym + (hc + hr));
              break;
            case Surface::Z_NEG:
              TNew[index] = (this_cell->kzp * *this_cell->k_up_Ptr->told_ptr/this_cell->dzp +
                  (hc + hr*pow(F,0.25))*Tair + q)/(this_cell->kzp/this_cell->dzp + (hc + hr));
              break;
            case Surface::Z_POS:
              TNew[index] = (this_cell->kzm * *this_cell->k_down_Ptr->told_ptr/this_cell->dzm +
                  (hc + hr*pow(F,0.25))*Tair + q)/(this_cell->kzm/this_cell->dzm + (hc + hr));
              break;
            }
            }
            break;
          }
          }
          break;
        case Cell::INTERIOR_AIR:
          TNew[index] = bcs.indoorTemp;
          break;

        case Cell::EXTERIOR_AIR:
          TNew[index] = bcs.outdoorTemp;
          break;
        default:
          {
          double theta = timestep/
            (this_cell->density*this_cell->specificHeat);

          double CXP = this_cell->cxp*theta;
          double CXM = this_cell->cxm*theta;
          double CZP = this_cell->czp*theta;
          double CZM = this_cell->czm*theta;
          double CYP = this_cell->cyp*theta;
          double CYM = this_cell->cym*theta;
          double Q = this_cell->heatGain*theta;

          if (foundation.numberOfDimensions == 3)
            TNew[index] = *this_cell->told_ptr*(1.0 + CXM + CZM + CYM - CXP - CZP - CYP)
                - *this_cell->i_down_Ptr->told_ptr*CXM
                + *this_cell->i_up_Ptr->told_ptr*CXP
                - *this_cell->k_down_Ptr->told_ptr*CZM
                + *this_cell->k_up_Ptr->told_ptr*CZP
                - *this_cell->j_down_Ptr->told_ptr*CYM
                + *this_cell->j_up_Ptr->told_ptr*CYP
                + Q;
          else if (foundation.numberOfDimensions == 2)
          {
            double CXPC = 0;
            double CXMC = 0;

            if (this_cell->i != 0)
            {
              CXPC = this_cell->cxp_c*theta/this_cell->r;
              CXMC = this_cell->cxm_c*theta/this_cell->r;
            }

            TNew[index] = *this_cell->told_ptr*(1.0 + CXMC + CXM + CZM - CXPC - CXP - CZP)
                - *this_cell->i_down_Ptr->told_ptr*(CXMC + CXM)
                + *this_cell->i_up_Ptr->told_ptr*(CXPC + CXP)
                - *this_cell->k_down_Ptr->told_ptr*CZM
                + *this_cell->k_up_Ptr->told_ptr*CZP
                + Q;
          }
          else
          {
            TNew[index] = *this_cell->told_ptr*(1.0 + CZM - CZP)
                - *this_cell->k_down_Ptr->told_ptr*CZM
                + *this_cell->k_up_Ptr->told_ptr*CZP
                + Q;
          }
          }
          break;
        }
  }
  TOld.assign(TNew.begin(), TNew.end());
}

void Ground::calculateMatrix(Foundation::NumericalScheme scheme)
{
  for (int index = 0; index < num_cells; index++)
  {
        Cell* this_cell = &domain.cell[index];
        int index_ip = index+domain.stepsize_i;
        int index_im = index-domain.stepsize_i;
        int index_jp = index+domain.stepsize_j;
        int index_jm = index-domain.stepsize_j;
        int index_kp = index+domain.stepsize_k;
        int index_km = index-domain.stepsize_k;

        double A, Aip, Aim, Ajp, Ajm, Akp, Akm, bVal = 0.0;

        switch (this_cell->cellType)
        {
        case Cell::BOUNDARY:
          {

          switch (this_cell->surfacePtr->boundaryConditionType)
          {
          case Surface::ZERO_FLUX:
            {
            switch (this_cell->surfacePtr->orientation)
            {
            case Surface::X_NEG:
              A = 1.0;
              Aip = -1.0;
              bVal = 0.0;

              setAmatValue(index,index,A);
              setAmatValue(index,index_ip,Aip);
              setbValue(index,bVal);
              break;
            case Surface::X_POS:
              A = 1.0;
              Aim = -1.0;
              bVal = 0.0;

              setAmatValue(index,index,A);
              setAmatValue(index,index_im,Aim);
              setbValue(index,bVal);
              break;
            case Surface::Y_NEG:
              A = 1.0;
              Ajp = -1.0;
              bVal = 0.0;

              setAmatValue(index,index,A);
              setAmatValue(index,index_jp,Ajp);
              setbValue(index,bVal);
              break;
            case Surface::Y_POS:
              A = 1.0;
              Ajm = -1.0;
              bVal = 0.0;

              setAmatValue(index,index,A);
              setAmatValue(index,index_jm,Ajm);
              setbValue(index,bVal);
              break;
            case Surface::Z_NEG:
              A = 1.0;
              Akp = -1.0;
              bVal = 0.0;

              setAmatValue(index,index,A);
              setAmatValue(index,index_kp,Akp);
              setbValue(index,bVal);
              break;
            case Surface::Z_POS:
              A = 1.0;
              Akm = -1.0;
              bVal = 0.0;

              setAmatValue(index,index,A);
              setAmatValue(index,index_km,Akm);
              setbValue(index,bVal);
              break;
            }
            }
            break;
          case Surface::CONSTANT_TEMPERATURE:
            A = 1.0;
            bVal = this_cell->surfacePtr->temperature;

            setAmatValue(index,index,A);
            setbValue(index,bVal);
            break;
          case Surface::INTERIOR_TEMPERATURE:
            A = 1.0;
            bVal = bcs.indoorTemp;

            setAmatValue(index,index,A);
            setbValue(index,bVal);
            break;
          case Surface::EXTERIOR_TEMPERATURE:
            A = 1.0;
            bVal = bcs.outdoorTemp;

            setAmatValue(index,index,A);
            setbValue(index,bVal);
            break;
          case Surface::INTERIOR_FLUX:
            {
            double& Tair = bcs.indoorTemp;
            double& q = this_cell->heatGain;

            double hc = getConvectionCoeff(*this_cell->told_ptr,
                    Tair,0.0,0.00208,false,this_cell->surfacePtr->tilt);
            double hr = getSimpleInteriorIRCoeff(this_cell->surfacePtr->emissivity,
                               *this_cell->told_ptr,Tair);

            switch (this_cell->surfacePtr->orientation)
            {
            case Surface::X_NEG:
              A = this_cell->kxp/this_cell->dxp + (hc + hr);
              Aip = -this_cell->kxp/this_cell->dxp;
              bVal = (hc + hr)*Tair + q;

              setAmatValue(index,index,A);
              setAmatValue(index,index_ip,Aip);
              setbValue(index,bVal);
              break;
            case Surface::X_POS:
              A = this_cell->kxm/this_cell->dxm + (hc + hr);
              Aim = -this_cell->kxm/this_cell->dxm;
              bVal = (hc + hr)*Tair + q;

              setAmatValue(index,index,A);
              setAmatValue(index,index_im,Aim);
              setbValue(index,bVal);
              break;
            case Surface::Y_NEG:
              A = this_cell->kyp/this_cell->dyp + (hc + hr);
              Ajp = -this_cell->kyp/this_cell->dyp;
              bVal = (hc + hr)*Tair + q;

              setAmatValue(index,index,A);
              setAmatValue(index,index_jp,Ajp);
              setbValue(index,bVal);
              break;
            case Surface::Y_POS:
              A = this_cell->kym/this_cell->dym + (hc + hr);
              Ajm = -this_cell->kym/this_cell->dym;
              bVal = (hc + hr)*Tair + q;

              setAmatValue(index,index,A);
              setAmatValue(index,index_jm,Ajm);
              setbValue(index,bVal);
              break;
            case Surface::Z_NEG:
              A = this_cell->kzp/this_cell->dzp + (hc + hr);
              Akp = -this_cell->kzp/this_cell->dzp;
              bVal = (hc + hr)*Tair + q;

              setAmatValue(index,index,A);
              setAmatValue(index,index_kp,Akp);
              setbValue(index,bVal);
              break;
            case Surface::Z_POS:
              A = this_cell->kzm/this_cell->dzm + (hc + hr);
              Akm = -this_cell->kzm/this_cell->dzm;
              bVal = (hc + hr)*Tair + q;

              setAmatValue(index,index,A);
              setAmatValue(index,index_km,Akm);
              setbValue(index,bVal);
              break;
            }
            }
            break;

          case Surface::EXTERIOR_FLUX:
            {
            double& Tair = bcs.outdoorTemp;
            double& v = bcs.localWindSpeed;
            double& eSky = bcs.skyEmissivity;
            double tilt = this_cell->surfacePtr->tilt;
            double F = getEffectiveExteriorViewFactor(eSky,tilt);
            double hc = getConvectionCoeff(*this_cell->told_ptr,Tair,v,foundation.surfaceRoughness,true,tilt);
            double hr = getExteriorIRCoeff(this_cell->surfacePtr->emissivity,*this_cell->told_ptr,Tair,eSky,tilt);
            double q = this_cell->heatGain;

            switch (this_cell->surfacePtr->orientation)
            {
            case Surface::X_NEG:
              A = this_cell->kxp/this_cell->dxp + (hc + hr);
              Aip = -this_cell->kxp/this_cell->dxp;
              bVal = (hc + hr*pow(F,0.25))*Tair + q;

              setAmatValue(index,index,A);
              setAmatValue(index,index_ip,Aip);
              setbValue(index,bVal);
              break;
            case Surface::X_POS:
              A = this_cell->kxm/this_cell->dxm + (hc + hr);
              Aim = -this_cell->kxm/this_cell->dxm;
              bVal = (hc + hr*pow(F,0.25))*Tair + q;

              setAmatValue(index,index,A);
              setAmatValue(index,index_im,Aim);
              setbValue(index,bVal);
              break;
            case Surface::Y_NEG:
              A = this_cell->kyp/this_cell->dyp + (hc + hr);
              Ajp = -this_cell->kyp/this_cell->dyp;
              bVal = (hc + hr*pow(F,0.25))*Tair + q;

              setAmatValue(index,index,A);
              setAmatValue(index,index_jp,Ajp);
              setbValue(index,bVal);
              break;
            case Surface::Y_POS:
              A = this_cell->kym/this_cell->dym + (hc + hr);
              Ajm = -this_cell->kym/this_cell->dym;
              bVal = (hc + hr*pow(F,0.25))*Tair + q;

              setAmatValue(index,index,A);
              setAmatValue(index,index_jm,Ajm);
              setbValue(index,bVal);
              break;
            case Surface::Z_NEG:
              A = this_cell->kzp/this_cell->dzp + (hc + hr);
              Akp = -this_cell->kzp/this_cell->dzp;
              bVal = (hc + hr*pow(F,0.25))*Tair + q;

              setAmatValue(index,index,A);
              setAmatValue(index,index_kp,Akp);
              setbValue(index,bVal);
              break;
            case Surface::Z_POS:
              A = this_cell->kzm/this_cell->dzm + (hc + hr);
              Akm = -this_cell->kzm/this_cell->dzm;
              bVal = (hc + hr*pow(F,0.25))*Tair + q;

              setAmatValue(index,index,A);
              setAmatValue(index,index_km,Akm);
              setbValue(index,bVal);
              break;
            }
            }
            break;
          }
          }
          break;
        case Cell::INTERIOR_AIR:
          A = 1.0;
          bVal = bcs.indoorTemp;

          setAmatValue(index,index,A);
          setbValue(index,bVal);
          break;
        case Cell::EXTERIOR_AIR:
          A = 1.0;
          bVal = bcs.outdoorTemp;

          setAmatValue(index,index,A);
          setbValue(index,bVal);
          break;
        default:
          {
          if (scheme == Foundation::NS_STEADY_STATE)
          {
            double CXP = this_cell->cxp;
            double CXM = this_cell->cxm;
            double CZP = this_cell->czp;
            double CZM = this_cell->czm;
            double CYP = this_cell->cyp;
            double CYM = this_cell->cym;
            double Q = this_cell->heatGain;

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

              setAmatValue(index,index,A);
              setAmatValue(index,index_ip,Aip);
              setAmatValue(index,index_im,Aim);
              setAmatValue(index,index_jp,Ajp);
              setAmatValue(index,index_jm,Ajm);
              setAmatValue(index,index_kp,Akp);
              setAmatValue(index,index_km,Akm);
              setbValue(index,bVal);
            }
            else if (foundation.numberOfDimensions == 2)
            {
              double CXPC = 0;
              double CXMC = 0;

              if (this_cell->i != 0)
              {
                CXPC = this_cell->cxp_c/this_cell->r;
                CXMC = this_cell->cxm_c/this_cell->r;
              }
              A = (CXMC + CXM + CZM - CXPC - CXP - CZP);
              Aim = (-CXMC - CXM);
              Aip = (CXPC + CXP);
              Akm = -CZM;
              Akp = CZP;

              bVal = -Q;

              setAmatValue(index,index,A);
              setAmatValue(index,index_ip,Aip);
              setAmatValue(index,index_im,Aim);
              setAmatValue(index,index_kp,Akp);
              setAmatValue(index,index_km,Akm);
              setbValue(index,bVal);
            }
            else
            {
              A = (CZM - CZP);
              Akm = -CZM;
              Akp = CZP;

              bVal = -Q;

              setAmatValue(index,index,A);
              setAmatValue(index,index_kp,Akp);
              setAmatValue(index,index_km,Akm);
              setbValue(index,bVal);
            }
          }
          else
          {
            double theta = timestep/
              (this_cell->density*this_cell->specificHeat);

            double f;
            if (scheme == Foundation::NS_IMPLICIT)
              f = 1.0;
            else
              f = 0.5;

            double CXP = this_cell->cxp*theta;
            double CXM = this_cell->cxm*theta;
            double CZP = this_cell->czp*theta;
            double CZM = this_cell->czm*theta;
            double CYP = this_cell->cyp*theta;
            double CYM = this_cell->cym*theta;
            double Q = this_cell->heatGain*theta;

            if (foundation.numberOfDimensions == 3)
            {
              A = (1.0 + f*(CXP + CZP + CYP - CXM - CZM - CYM));
              Aim = f*CXM;
              Aip = f*(-CXP);
              Akm = f*CZM;
              Akp = f*(-CZP);
              Ajm = f*CYM;
              Ajp = f*(-CYP);

              bVal = *this_cell->told_ptr*(1.0 + (1-f)*(CXM + CZM + CYM - CXP - CZP - CYP))
                 - *this_cell->i_down_Ptr->told_ptr*(1-f)*CXM
                 + *this_cell->i_up_Ptr->told_ptr*(1-f)*CXP
                 - *this_cell->j_down_Ptr->told_ptr*(1-f)*CZM
                 + *this_cell->j_up_Ptr->told_ptr*(1-f)*CZP
                 - *this_cell->k_down_Ptr->told_ptr*(1-f)*CYM
                 + *this_cell->k_up_Ptr->told_ptr*(1-f)*CYP
                 + Q;

              setAmatValue(index,index,A);
              setAmatValue(index,index_ip,Aip);
              setAmatValue(index,index_im,Aim);
              setAmatValue(index,index_jp,Ajp);
              setAmatValue(index,index_jm,Ajm);
              setAmatValue(index,index_kp,Akp);
              setAmatValue(index,index_km,Akm);
              setbValue(index,bVal);
            }
            else if (foundation.numberOfDimensions == 2)
            {
              double CXPC = 0;
              double CXMC = 0;

              if (this_cell->i != 0)
              {
                CXPC = this_cell->cxp_c*theta/this_cell->r;
                CXMC = this_cell->cxm_c*theta/this_cell->r;
              }
              A = (1.0 + f*(CXPC + CXP + CZP - CXMC - CXM - CZM));
              Aim = f*(CXMC + CXM);
              Aip = f*(-CXPC - CXP);
              Akm = f*CZM;
              Akp = f*(-CZP);

              bVal = *this_cell->told_ptr*(1.0 + (1-f)*(CXMC + CXM + CZM - CXPC - CXP - CZP))
                 - *this_cell->i_down_Ptr->told_ptr*(1-f)*(CXMC + CXM)
                 + *this_cell->i_up_Ptr->told_ptr*(1-f)*(CXPC + CXP)
                 - *this_cell->k_down_Ptr->told_ptr*(1-f)*CZM
                 + *this_cell->k_up_Ptr->told_ptr*(1-f)*CZP
                 + Q;

              setAmatValue(index,index,A);
              setAmatValue(index,index_ip,Aip);
              setAmatValue(index,index_im,Aim);
              setAmatValue(index,index_kp,Akp);
              setAmatValue(index,index_km,Akm);
              setbValue(index,bVal);
            }
            else
            {
              A = (1.0 + f*(CZP - CZM));
              Akm = f*CZM;
              Akp = f*(-CZP);

              bVal = *this_cell->told_ptr*(1.0 + (1-f)*(CZM - CZP))
                 - *this_cell->k_down_Ptr->told_ptr*(1-f)*CZM
                 + *this_cell->k_up_Ptr->told_ptr*(1-f)*CZP
                 + Q;

              setAmatValue(index,index,A);
              setAmatValue(index,index_kp,Akp);
              setAmatValue(index,index_km,Akm);
              setbValue(index,bVal);
            }
          }
          }
          break;
        }
  }

  solveLinearSystem();

  // Read solution into temperature matrix
  TNew = getXvalues();
  // Update old values for next timestep
  TOld.assign(TNew.begin(), TNew.end());
  clearAmat();
}

void Ground::calculateADI(int dim)
{
  std::size_t dest_index;
  for (size_t index = 0; index < num_cells; index++)
  {
        Cell* this_cell = &domain.cell[index];
        index = this_cell->index;
        dest_index = domain.dest_index_vector[index][dim-1];

        double A{0.0}, Ap{0.0}, Am{0.0}, bVal{0.0};


        switch (this_cell->cellType)
        {
        case Cell::BOUNDARY:
          {
          switch (this_cell->surfacePtr->boundaryConditionType)
          {
          case Surface::ZERO_FLUX:
            {
            switch (this_cell->surfacePtr->orientation)
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
                bVal = *this_cell->i_up_Ptr->told_ptr;
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
                bVal = *this_cell->i_down_Ptr->told_ptr;
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
                bVal = *this_cell->j_up_Ptr->told_ptr;
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
                bVal = *this_cell->j_down_Ptr->told_ptr;
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
                bVal = *this_cell->k_up_Ptr->told_ptr;
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
                bVal = *this_cell->k_down_Ptr->told_ptr;
              }
              break;
            }
            }
            break;
          case Surface::CONSTANT_TEMPERATURE:
            A = 1.0;
            bVal = this_cell->surfacePtr->temperature;
            break;
          case Surface::INTERIOR_TEMPERATURE:
            A = 1.0;
            bVal = bcs.indoorTemp;
            break;
          case Surface::EXTERIOR_TEMPERATURE:
            A = 1.0;
            bVal = bcs.outdoorTemp;
            break;
          case Surface::INTERIOR_FLUX:
            {
            double& Tair = bcs.indoorTemp;
            double& q = this_cell->heatGain;

            double hc = getConvectionCoeff(*this_cell->told_ptr,
                    Tair,0.0,0.00208,false,this_cell->surfacePtr->tilt);
            double hr = getSimpleInteriorIRCoeff(this_cell->surfacePtr->emissivity,
                               *this_cell->told_ptr,Tair);

            switch (this_cell->surfacePtr->orientation)
            {
            case Surface::X_NEG:
              A = this_cell->kxp/this_cell->dxp + (hc + hr);
              if (dim == 1)
              {
                Ap = -this_cell->kxp/this_cell->dxp;
                bVal = (hc + hr)*Tair + q;
              }
              else
              {
                Ap = 0.0;
                bVal = *this_cell->i_up_Ptr->told_ptr*this_cell->kxp/this_cell->dxp + (hc + hr)*Tair + q;
              }
              break;
            case Surface::X_POS:
              A = this_cell->kxm/this_cell->dxm + (hc + hr);
              if (dim == 1)
              {
                Am = -this_cell->kxm/this_cell->dxm;
                bVal = (hc + hr)*Tair + q;
              }
              else
              {
                Am = 0.0;
                bVal = *this_cell->i_down_Ptr->told_ptr*this_cell->kxm/this_cell->dxm + (hc + hr)*Tair + q;
              }
              break;
            case Surface::Y_NEG:
              A = this_cell->kyp/this_cell->dyp + (hc + hr);
              if (dim == 2)
              {
                Ap = -this_cell->kyp/this_cell->dyp;
                bVal = (hc + hr)*Tair + q;
              }
              else
              {
                Ap = 0.0;
                bVal = *this_cell->j_up_Ptr->told_ptr*this_cell->kyp/this_cell->dyp + (hc + hr)*Tair + q;
              }
              break;
            case Surface::Y_POS:
              A = this_cell->kym/this_cell->dym + (hc + hr);
              if (dim == 2)
              {
                Am = -this_cell->kym/this_cell->dym;
                bVal = (hc + hr)*Tair + q;
              }
              else
              {
                Am = 0.0;
                bVal = *this_cell->j_down_Ptr->told_ptr*this_cell->kym/this_cell->dym + (hc + hr)*Tair + q;
              }
              break;
            case Surface::Z_NEG:
              A = this_cell->kzp/this_cell->dzp + (hc + hr);
              if (dim == 3)
              {
                Ap = -this_cell->kzp/this_cell->dzp;
                bVal = (hc + hr)*Tair + q;
              }
              else
              {
                Ap = 0.0;
                bVal = *this_cell->k_up_Ptr->told_ptr*this_cell->kzp/this_cell->dzp + (hc + hr)*Tair + q;
              }
              break;
            case Surface::Z_POS:
              A = this_cell->kzm/this_cell->dzm + (hc + hr);
              if (dim == 3)
              {
                Am = -this_cell->kzm/this_cell->dzm;
                bVal = (hc + hr)*Tair + q;
              }
              else
              {
                Am = 0.0;
                bVal = *this_cell->k_down_Ptr->told_ptr*this_cell->kzm/this_cell->dzm + (hc + hr)*Tair + q;
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
            double tilt = this_cell->surfacePtr->tilt;
            double F = getEffectiveExteriorViewFactor(eSky,tilt);
            double hc = getConvectionCoeff(*this_cell->told_ptr,Tair,v,foundation.surfaceRoughness,true,tilt);
            double hr = getExteriorIRCoeff(this_cell->surfacePtr->emissivity,*this_cell->told_ptr,Tair,eSky,tilt);
            double q = this_cell->heatGain;

            switch (this_cell->surfacePtr->orientation)
            {
            case Surface::X_NEG:
              A = this_cell->kxp/this_cell->dxp + (hc + hr);
              if (dim == 1)
              {
                Ap = -this_cell->kxp/this_cell->dxp;
                bVal = (hc + hr*pow(F,0.25))*Tair + q;
              }
              else
              {
                Ap = 0.0;
                bVal = *this_cell->i_up_Ptr->told_ptr*this_cell->kxp/this_cell->dxp + (hc + hr*pow(F,0.25))*Tair + q;
              }
              break;
            case Surface::X_POS:
              A = this_cell->kxm/this_cell->dxm + (hc + hr);
              if (dim == 1)
              {
                Am = -this_cell->kxm/this_cell->dxm;
                bVal = (hc + hr*pow(F,0.25))*Tair + q;
              }
              else
              {
                Am = 0.0;
                bVal = *this_cell->i_down_Ptr->told_ptr*this_cell->kxm/this_cell->dxm + (hc + hr*pow(F,0.25))*Tair + q;
              }
              break;
            case Surface::Y_NEG:
              A = this_cell->kyp/this_cell->dyp + (hc + hr);
              if (dim == 2)
              {
                Ap = -this_cell->kyp/this_cell->dyp;
                bVal = (hc + hr*pow(F,0.25))*Tair + q;
              }
              else
              {
                Ap = 0.0;
                bVal = *this_cell->j_up_Ptr->told_ptr*this_cell->kyp/this_cell->dyp + (hc + hr*pow(F,0.25))*Tair + q;
              }
              break;
            case Surface::Y_POS:
              A = this_cell->kym/this_cell->dym + (hc + hr);
              if (dim == 2)
              {
                Am = -this_cell->kym/this_cell->dym;
                bVal = (hc + hr*pow(F,0.25))*Tair + q;
              }
              else
              {
                Am = 0.0;
                bVal = *this_cell->j_down_Ptr->told_ptr*this_cell->kym/this_cell->dym + (hc + hr*pow(F,0.25))*Tair + q;
              }
              break;
            case Surface::Z_NEG:
              A = this_cell->kzp/this_cell->dzp + (hc + hr);
              if (dim == 3)
              {
                Ap = -this_cell->kzp/this_cell->dzp;
                bVal = (hc + hr*pow(F,0.25))*Tair + q;
              }
              else
              {
                Ap = 0.0;
                bVal = *this_cell->k_up_Ptr->told_ptr*this_cell->kzp/this_cell->dzp + (hc + hr*pow(F,0.25))*Tair + q;
              }
              break;
            case Surface::Z_POS:
              A = this_cell->kzm/this_cell->dzm + (hc + hr);
              if (dim == 3)
              {
                Am = -this_cell->kzm/this_cell->dzm;
                bVal = (hc + hr*pow(F,0.25))*Tair + q;
              }
              else
              {
                Am = 0.0;
                bVal = *this_cell->k_down_Ptr->told_ptr*this_cell->kzm/this_cell->dzm + (hc + hr*pow(F,0.25))*Tair + q;
              }
              break;
            }
            }
            break;
          }
          }
          break;
        case Cell::INTERIOR_AIR:
          A = 1.0;
          bVal = bcs.indoorTemp;
          break;
        case Cell::EXTERIOR_AIR:
          A = 1.0;
          bVal = bcs.outdoorTemp;
          break;
        default:
          {
          double theta = timestep / (foundation.numberOfDimensions
                                      *this_cell->density*this_cell->specificHeat);

          double CXP = this_cell->cxp*theta;
          double CXM = this_cell->cxm*theta;
          double CZP = this_cell->czp*theta;
          double CZM = this_cell->czm*theta;
          double CYP = this_cell->cyp*theta;
          double CYM = this_cell->cym*theta;
          double Q = this_cell->heatGain*theta;

          double f = foundation.fADI;

          if (foundation.numberOfDimensions == 3)
          {
            if (dim == 1) // x
            {
              A = 1.0 + (3 - 2*f)*(CXP - CXM);
              Am = (3 - 2*f)*CXM;
              Ap = (3 - 2*f)*(-CXP);

              bVal = *this_cell->told_ptr*(1.0 + f*(CZM + CYM - CZP - CYP))
                   - *this_cell->k_down_Ptr->told_ptr*f*CZM
                   + *this_cell->k_up_Ptr->told_ptr*f*CZP
                   - *this_cell->j_down_Ptr->told_ptr*f*CYM
                   + *this_cell->j_up_Ptr->told_ptr*f*CYP
                   + Q;
            }
            else if (dim == 2) // y
            {
              A = (1.0 + (3 - 2*f)*(CYP - CYM));
              Am = (3 - 2*f)*CYM;
              Ap = (3 - 2*f)*(-CYP);

              bVal = *this_cell->told_ptr*(1.0 + f*(CXM + CZM - CXP - CZP))
                   - *this_cell->i_down_Ptr->told_ptr*f*CXM
                   + *this_cell->i_up_Ptr->told_ptr*f*CXP
                   - *this_cell->k_down_Ptr->told_ptr*f*CZM
                   + *this_cell->k_up_Ptr->told_ptr*f*CZP
                   + Q;
            }
            else //if (dim == 3) // z
            {
              A = (1.0 + (3 - 2*f)*(CZP - CZM));
              Am = (3 - 2*f)*CZM;
              Ap = (3 - 2*f)*(-CZP);

              bVal = *this_cell->told_ptr*(1.0 + f*(CXM + CYM - CXP - CYP))
                   - *this_cell->i_down_Ptr->told_ptr*f*CXM
                   + *this_cell->i_up_Ptr->told_ptr*f*CXP
                   - *this_cell->j_down_Ptr->told_ptr*f*CYM
                   + *this_cell->j_up_Ptr->told_ptr*f*CYP
                   + Q;
            }

          }
          else if (foundation.numberOfDimensions == 2)
          {
            double CXPC = 0;
            double CXMC = 0;
            if (this_cell->i != 0)
            {
              CXPC = this_cell->cxp_c*theta/this_cell->r;
              CXMC = this_cell->cxm_c*theta/this_cell->r;
            }
            if (dim == 1) // x
            {
              A = 1.0 + (2 - f)*(CXPC + CXP - CXMC - CXM);
              Am = (2 - f)*(CXMC + CXM);
              Ap = (2 - f)*(-CXPC - CXP);

              bVal = *this_cell->told_ptr*(1.0 + f*(CZM - CZP))
                   - *this_cell->k_down_Ptr->told_ptr*f*CZM
                   + *this_cell->k_up_Ptr->told_ptr*f*CZP
                   + Q;
            }
            else //if (dim == 3) // z
            {
              A = 1.0 + (2 - f)*(CZP - CZM);
              Am = (2 - f)*CZM;
              Ap = (2 - f)*(-CZP);

              bVal = *this_cell->told_ptr*(1.0 + f*(CXMC + CXM - CXPC - CXP))
                   - *this_cell->i_down_Ptr->told_ptr*f*(CXMC + CXM)
                   + *this_cell->i_up_Ptr->told_ptr*f*(CXPC + CXP)
                   + Q;
            }
          }
          else
          {
            A = 1.0 + CZP - CZM;
            Am = CZM;
            Ap = -CZP;

            bVal = *this_cell->told_ptr + Q;
          }
          }
          break;
        }
        setValuesADI(dest_index, Am, A, Ap, bVal);
  }

  solveLinearSystem();

  for (size_t index = 0; index < num_cells; ++index)
  {
    dest_index = domain.dest_index_vector[index][dim-1];

      // Read solution into temperature matrix
    TNew[index] = getxValue(dest_index);
  }
  // Update old values for next timestep
  TOld.assign(TNew.begin(), TNew.end());

  clearAmat();
}

void Ground::calculate(BoundaryConditions& boundaryConidtions, double ts)
{
  bcs = boundaryConidtions;
  timestep = ts;
  // update boundary conditions
  setSolarBoundaryConditions();
  setInteriorRadiationBoundaryConditions();

  // Calculate Temperatures
  switch(foundation.numericalScheme)
  {
  case Foundation::NS_ADE:
    calculateADE();
    break;
  case Foundation::NS_EXPLICIT:
    calculateExplicit();
    break;
  case Foundation::NS_ADI:
    {
    if (foundation.numberOfDimensions > 1)
      calculateADI(1);
    if (foundation.numberOfDimensions == 3)
      calculateADI(2);
    calculateADI(3);
    }
    break;
  case Foundation::NS_IMPLICIT:
    calculateMatrix(Foundation::NS_IMPLICIT);
    break;
  case Foundation::NS_CRANK_NICOLSON:
    calculateMatrix(Foundation::NS_CRANK_NICOLSON);
    break;
  case Foundation::NS_STEADY_STATE:
    calculateMatrix(Foundation::NS_STEADY_STATE);
    break;
  }

}

void Ground::setAmatValue(const int i,const int j,const double val)
{
  if ((foundation.numericalScheme == Foundation::NS_ADI ||
    foundation.numberOfDimensions == 1) && TDMA)
  {
    if (j < i)
      a1[i] = val;
    else if (j == i)
      a2[i] = val;
    else
      a3[i] = val;
  }
  else
  {
    tripletList.emplace_back(i,j,val);
  }
}

void Ground::setbValue(const int i,const double val)
{
  if ((foundation.numericalScheme == Foundation::NS_ADI ||
    foundation.numberOfDimensions == 1) && TDMA)
  {
    b_[i] = val;
  }
  else
  {
    b(i) = val;
  }
}

void Ground::setValuesADI(const std::size_t & index, const double & Am, const double & A,
                          const double & Ap, const double & bVal) {
    a1[index] = Am;
    a2[index] = A;
    a3[index] = Ap;
    b_[index] = bVal;
};

void Ground::solveLinearSystem()
{
  if ((foundation.numericalScheme == Foundation::NS_ADI ||
    foundation.numberOfDimensions == 1) && TDMA)
  {
    solveTDM(a1,a2,a3,b_,x_);
  }
  else
  {
    int iters;
    double residual;

    bool success;

    Amat.setFromTriplets(tripletList.begin(), tripletList.end());
    pSolver->compute(Amat);
    x = pSolver->solveWithGuess(b,x);
    int status = pSolver->info();

//    Eigen::saveMarket(Amat, "Amat.mtx");
//    Eigen::saveMarketVector(b, "b.mtx");
    success = status == Eigen::Success;
    if (!success) {
      iters = pSolver->iterations();
      residual = pSolver->error();

      std::stringstream ss;
      ss << "Solution did not converge after " << iters << " iterations. The final residual was: (" << residual << ").";
      showMessage(MSG_ERR, ss.str());
    }
  }
}

void Ground::clearAmat()
{
  if ((foundation.numericalScheme == Foundation::NS_ADI ||
    foundation.numberOfDimensions == 1) && TDMA)
  {
    std::fill(a1.begin(), a1.end(), 0.0);
    std::fill(a2.begin(), a2.end(), 0.0);
    std::fill(a3.begin(), a3.end(), 0.0);
    std::fill(b_.begin(), b_.end(), 0.0);

  }
  else
  {
    tripletList.clear();
    tripletList.reserve(nX*nY*nZ*(1+2*foundation.numberOfDimensions));
  }
}

double Ground::getxValue(const int i)
{
  if ((foundation.numericalScheme == Foundation::NS_ADI ||
    foundation.numberOfDimensions == 1) && TDMA)
  {
    return x_[i];
  }
  else
  {
    return x(i);
  }
}

std::vector<double> Ground::getXvalues()
{
  if ((foundation.numericalScheme == Foundation::NS_ADI ||
    foundation.numberOfDimensions == 1) && TDMA)
  {
    return x_;
  }
  else
  {
    std::vector<double> v(x.data(), x.data() + x.rows());
    return v;
  }
}

double Ground::getConvectionCoeff(double Tsurf,
                  double Tamb,
                  double Vair,
                    double roughness,
                  bool isExterior,
                  double tilt)
{
  if (foundation.convectionCalculationMethod == Foundation::CCM_AUTO)
    return getDOE2ConvectionCoeff(tilt,0.0,0.0,Tsurf,Tamb,Vair,roughness);
  else //if (foundation.convectionCalculationMethod == Foundation::CCM_CONSTANT_COEFFICIENT)
  {
    if (isExterior)
      return foundation.exteriorConvectiveCoefficient;
    else
      return foundation.interiorConvectiveCoefficient;
  }
}

double Ground::getSurfaceArea(Surface::SurfaceType surfaceType)
{
  double totalArea = 0;

  // Find surface(s)
  for (size_t s = 0; s < foundation.surfaces.size(); s++)
  {
    if (foundation.surfaces[s].type == surfaceType)
    {
      Surface surface;
      surface = foundation.surfaces[s];

      totalArea += surface.area;
    }
  }

  return totalArea;
}

void Ground::calculateSurfaceAverages(){
  for (auto output : groundOutput.outputMap) {
    Surface::SurfaceType surface = output.first;
    std::vector<GroundOutput::OutputType> outTypes = output.second;

    double constructionRValue = 0.0;
    double surfaceArea = foundation.surfaceAreas[surface];

    if (surface == Surface::ST_SLAB_CORE) {
      constructionRValue = foundation.slab.totalResistance();
    }
    else if (surface == Surface::ST_SLAB_PERIM) {
      constructionRValue = foundation.slab.totalResistance();
    }
    else if (surface == Surface::ST_WALL_INT) {
      constructionRValue = foundation.wall.totalResistance();
    }

    double totalHeatTransferRate = 0.0;
    //double TA = 0;
    double HA = 0.0;
    double totalArea = 0.0;

    double& Tair = bcs.indoorTemp;

    if (foundation.hasSurface[surface]) {
      // Find surface(s)
      for (size_t s = 0; s < foundation.surfaces.size(); s++)
      {
        double tilt = foundation.surfaces[s].tilt;
        if (foundation.surfaces[s].type == surface)
        {

          #ifdef PRNTSURF
            std::ofstream output;
            output.open("surface.csv");
            output  << "x, T, h, q, dx\n";
          #endif

          for (auto index : foundation.surfaces[s].indices)
          {
            Cell* this_cell = &domain.cell[index];
            double h = getConvectionCoeff(TNew[index],Tair,0.0,0.00208,false,tilt)
                 + getSimpleInteriorIRCoeff(this_cell->surfacePtr->emissivity,
                     TNew[index],Tair);

            double& A = this_cell->area;

            totalArea += A;
            totalHeatTransferRate += h*A*(Tair - TNew[index]);
            //TA += TNew[index]*A;
            HA += h*A;

            #ifdef PRNTSURF
              output <<
                domain.meshX.centers[i] << ", " <<
                TNew[index] << ", " <<
                h << ", " <<
                h*(Tair - TNew[index]) << ", " <<
                domain.meshX.deltas[i] << "\n";
            #endif

          }

          #ifdef PRNTSURF
            output.close();
          #endif

        }
      }
    }

    if (totalArea > 0.0) {
      double Tavg = Tair - totalHeatTransferRate/HA;
      double hAvg = HA/totalArea;

      groundOutput.outputValues[{surface,GroundOutput::OT_TEMP}] = Tavg;
      groundOutput.outputValues[{surface,GroundOutput::OT_FLUX}] = totalHeatTransferRate/totalArea;
      groundOutput.outputValues[{surface,GroundOutput::OT_RATE}] = totalHeatTransferRate/totalArea*surfaceArea;
      groundOutput.outputValues[{surface,GroundOutput::OT_CONV}] = hAvg;

      groundOutput.outputValues[{surface,GroundOutput::OT_EFF_TEMP}] = Tair - (totalHeatTransferRate/totalArea)*(constructionRValue+1/hAvg) - 273.15;
    }
    else {
      groundOutput.outputValues[{surface,GroundOutput::OT_TEMP}] = Tair;
      groundOutput.outputValues[{surface,GroundOutput::OT_FLUX}] = 0.0;
      groundOutput.outputValues[{surface,GroundOutput::OT_RATE}] = 0.0;
      groundOutput.outputValues[{surface,GroundOutput::OT_CONV}] = 0.0;

      groundOutput.outputValues[{surface,GroundOutput::OT_EFF_TEMP}] = Tair - 273.15;
    }
  }
}

double Ground::getSurfaceAverageValue(std::pair<Surface::SurfaceType, GroundOutput::OutputType> output)
{
  return groundOutput.outputValues[output];
}

std::vector<double> Ground::calculateHeatFlux(Cell* this_cell)
{
  std::vector<double> Qflux;
  double Qx = 0;
  double Qy = 0;
  double Qz = 0;

  double CXP = 0;
  double CXM = 0;
  double CYP = 0;
  double CYM = 0;
  double CZP = -this_cell->kzp*this_cell->dzm/(this_cell->dzp+this_cell->dzm)/this_cell->dzp;
  double CZM = -this_cell->kzm*this_cell->dzp/(this_cell->dzp+this_cell->dzm)/this_cell->dzm;

  if (foundation.numberOfDimensions > 1)
  {
    CXP = -this_cell->kxp*this_cell->dxm/(this_cell->dxp+this_cell->dxm)/this_cell->dxp;
    CXM = -this_cell->kxm*this_cell->dxp/(this_cell->dxp+this_cell->dxm)/this_cell->dxm;
  }


  if (foundation.numberOfDimensions == 3)
  {
    CYP = -this_cell->kyp*this_cell->dym/(this_cell->dyp+this_cell->dym)/this_cell->dyp;
    CYM = -this_cell->kym*this_cell->dyp/(this_cell->dyp+this_cell->dym)/this_cell->dym;
  }

  double DTXP = 0;
  double DTXM = 0;
  double DTYP = 0;
  double DTYM = 0;
  double DTZP = 0;
  double DTZM = 0;

  std::size_t index = this_cell->index;

  if (this_cell->i != nX - 1)
    DTXP = TNew[index+domain.stepsize_i]-TNew[index];

  if (this_cell->i != 0)
    DTXM = TNew[index]-TNew[index-domain.stepsize_i];

  if (this_cell->j != nY - 1)
    DTYP = TNew[index+domain.stepsize_j]-TNew[index];

  if (this_cell->j != 0)
    DTYM = TNew[index]-TNew[index-domain.stepsize_j];

  if (this_cell->k != nZ - 1)
    DTZP = TNew[index+domain.stepsize_k]-TNew[index];

  if (this_cell->k != 0)
    DTZM = TNew[index]-TNew[index-domain.stepsize_k];

  switch (this_cell->cellType)
  {
    case Cell::BOUNDARY:
      {
        switch (this_cell->surfacePtr->orientation)
        {
          case Surface::X_NEG:
            {
              CXP = -this_cell->kxp/this_cell->dxp;
              CXM = 0;
            }
          break;
          case Surface::X_POS:
            {
              CXP = 0;
              CXM = -this_cell->kxm/this_cell->dxm;
            }
          break;
          case Surface::Y_NEG:
            {
              CYP = -this_cell->kyp/this_cell->dyp;
              CYM = 0;
            }
          break;
          case Surface::Y_POS:
            {
              CYP = 0;
              CYM = -this_cell->kym/this_cell->dym;
            }
          break;
          case Surface::Z_NEG:
            {
              CZP = -this_cell->kzp/this_cell->dzp;
              CZM = 0;
            }
          break;
          case Surface::Z_POS:
            {
              CZP = 0;
              CZM = -this_cell->kzm/this_cell->dzm;
            }
          break;
        }
        Qx = CXP*DTXP + CXM*DTXM;
        Qy = CYP*DTYP + CYM*DTYM;
        Qz = CZP*DTZP + CZM*DTZM;
      }
      break;
    case Cell::INTERIOR_AIR:
      break;
    case Cell::EXTERIOR_AIR:
      break;
    case Cell::ZERO_THICKNESS:
      {
        //int numZeroDims = domain.getNumZeroDims(i,j,k);

        std::vector<double> Qm;
        std::vector<double> Qp;

        if (isEqual(domain.meshX.deltas[this_cell->i], 0.0))
        {
          Qm = calculateHeatFlux(this_cell->i_down_Ptr);
          Qp = calculateHeatFlux(this_cell->i_up_Ptr);
        }
        if (isEqual(domain.meshY.deltas[this_cell->j], 0.0))
        {
          Qm = calculateHeatFlux(this_cell->j_down_Ptr);
          Qp = calculateHeatFlux(this_cell->j_up_Ptr);
        }
        if (isEqual(domain.meshZ.deltas[this_cell->k], 0.0))
        {
          Qm = calculateHeatFlux(this_cell->k_down_Ptr);
          Qp = calculateHeatFlux(this_cell->k_up_Ptr);
        }

        Qx = (Qm[0] + Qp[0])*0.5;
        Qy = (Qm[1] + Qp[1])*0.5;
        Qz = (Qm[2] + Qp[2])*0.5;
      }
      break;
    default:
      {
        Qx = CXP*DTXP + CXM*DTXM;
        Qy = CYP*DTYP + CYM*DTYM;
        Qz = CZP*DTZP + CZM*DTZM;
      }
    break;
  }

  Qflux.push_back(Qx);
  Qflux.push_back(Qy);
  Qflux.push_back(Qz);

  return Qflux;
}

void Ground::calculateBoundaryLayer()
{
  Foundation fd = foundation;

  BoundaryConditions preBCs;
  preBCs.localWindSpeed = 0;
  preBCs.outdoorTemp = 273.15;
  preBCs.indoorTemp = 293.15;
  fd.coordinateSystem = Foundation::CS_CARTESIAN;
  fd.numberOfDimensions = 2;
  fd.reductionStrategy = Foundation::RS_AP;
  fd.numericalScheme = Foundation::NS_STEADY_STATE;
  fd.farFieldWidth = 100;

  Ground pre(fd);
  pre.buildDomain();
  pre.calculate(preBCs);

  std::vector<double> x2s;
  std::vector<double> fluxSums;

  double fluxSum = 0.0;

  double x1_0 = 0.0;

  bool firstIndex = true;

  size_t i_min = pre.domain.meshX.getNearestIndex(boost::geometry::area(foundation.polygon)/
      boost::geometry::perimeter(foundation.polygon));

  size_t k = pre.domain.meshZ.getNearestIndex(0.0);

  size_t j = pre.nY/2;

  for (size_t i = i_min; i < pre.nX; i++)
  {
    std::size_t index = i + pre.nX*j + pre.nX*pre.nY*k;
    double Qz = pre.calculateHeatFlux(&pre.domain.cell[index])[2];
    double x1 = pre.domain.meshX.dividers[i];
    double x2 = pre.domain.meshX.dividers[i+1];

    if (Qz > 0.0)
    {
      fluxSum += std::max(Qz,0.0)*(x2-x1);

      if (firstIndex)
        x1_0 = x1;
      x2s.push_back(x2);
      fluxSums.push_back(fluxSum);

      firstIndex = false;
    }

  }

  //std::ofstream output;
  //output.open("Boundary.csv");

  //output << 0.0 << ", " << 0.0 << "\n";

  boundaryLayer.push_back(std::make_pair(0,0));

  for (std::size_t i = 0; i < fluxSums.size() - 1; i++) // last cell is a zero-thickness cell, so don't include it.
  {
    //output << x2s[i] - x1_0 << ", " << fluxSums[i]/fluxSum << "\n";
    boundaryLayer.push_back(std::make_pair(x2s[i] - x1_0,fluxSums[i]/fluxSum));
  }

}

double Ground::getBoundaryValue(double dist)
{
  double val = 0.0;
  if (dist > boundaryLayer[boundaryLayer.size()-1].first)
    val = 1.0;
  else
  {
    for (std::size_t i = 0; i < boundaryLayer.size()-1; i++)
    {
      if (dist >= boundaryLayer[i].first && dist < boundaryLayer[i+1].first)
      {
        double m = (boundaryLayer[i+1].first - boundaryLayer[i].first)/
            (boundaryLayer[i+1].second - boundaryLayer[i].second);
        val = (dist - boundaryLayer[i].first)/m + boundaryLayer[i].second;
        continue;
      }
    }
  }
  return val;
}

double Ground::getBoundaryDistance(double val)
{
  double dist = 0.0;
  if (val > 1.0 || val < 0.0)
  {
    showMessage(MSG_ERR, "Boundary value passed not between 0.0 and 1.0.");
  }
  else
  {
    for (std::size_t i = 0; i < boundaryLayer.size()-1; i++)
    {
      if (val >= boundaryLayer[i].second && val < boundaryLayer[i+1].second)
      {
        double m = (boundaryLayer[i+1].second - boundaryLayer[i].second)/
            (boundaryLayer[i+1].first - boundaryLayer[i].first);
        dist = (val - boundaryLayer[i].second)/m + boundaryLayer[i].first;
        continue;
      }
    }
  }
  return dist;
}

void Ground::setNewBoundaryGeometry()
{
  double area = boost::geometry::area(foundation.polygon);
  double perimeter = boost::geometry::perimeter(foundation.polygon);
  double interiorPerimeter = 0.0;

  std::size_t nV = foundation.polygon.outer().size();
  for (std::size_t v = 0; v < nV; v++)
  {
    std::size_t vPrev, vNext, vNext2;

    if (v == 0)
      vPrev = nV - 1;
    else
      vPrev = v - 1;

    if (v == nV -1)
      vNext = 0;
    else
      vNext = v + 1;

    if (v == nV - 2)
      vNext2 = 0;
    else if (v == nV -1)
      vNext2 = 1;
    else
      vNext2 = v + 2;

    Point a = foundation.polygon.outer()[vPrev];
    Point b = foundation.polygon.outer()[v];
    Point c = foundation.polygon.outer()[vNext];
    Point d = foundation.polygon.outer()[vNext2];

    // Correct U-turns
    if (foundation.isExposedPerimeter[vPrev] && foundation.isExposedPerimeter[v] && foundation.isExposedPerimeter[vNext])
    {
      if (isEqual(getAngle(a,b,c) + getAngle(b,c,d),PI))
      {
        double AB = getDistance(a,b);
        double BC = getDistance(b,c);
        double CD = getDistance(c,d);
        double edgeDistance = BC;
        double reductionDistance = std::min(AB,CD);
        double reductionValue = 1 - getBoundaryValue(edgeDistance);
        perimeter -= 2*reductionDistance*reductionValue;
      }
    }

    if (foundation.isExposedPerimeter[vPrev] && foundation.isExposedPerimeter[v])
    {
      double alpha = getAngle(a,b,c);
      double A = getDistance(a,b);
      double B = getDistance(b,c);


      if (sin(alpha) > 0)
      {
        double f = getBoundaryDistance(1-sin(alpha/2)/(1+cos(alpha/2)))/sin(alpha/2);

        // Chamfer
        double d = f/cos(alpha/2);
        if (A < d || B < d)
        {
          A = std::min(A,B);
          B = std::min(A,B);
        }
        else
        {
          A = d;
          B = d;
        }
        double C = sqrt(A*A + B*B - 2*A*B*cos(alpha));

        perimeter += C - (A + B);

      }
    }

    if (!foundation.isExposedPerimeter[v])
    {
      interiorPerimeter += getDistance(b,c);
    }

  }

  foundation.reductionStrategy = Foundation::RS_CUSTOM;
  foundation.twoParameters = false;
  foundation.reductionLength2 = area/(perimeter - interiorPerimeter);

}

void Ground::setSolarBoundaryConditions()
{
  if (foundation.numberOfDimensions == 1) {
    return;
  }
  for (std::size_t s = 0; s < foundation.surfaces.size() ; s++)
  {
    if (foundation.surfaces[s].type == Surface::ST_GRADE
        || foundation.surfaces[s].type == Surface::ST_WALL_EXT)
    {

      double& azi = bcs.solarAzimuth;
      double& alt = bcs.solarAltitude;
      double& qDN = bcs.directNormalFlux;
      double& qDH = bcs.diffuseHorizontalFlux;
      double qGH = cos(PI/2 - alt)*qDN + qDH;
      double pssf;
      double q;

      double incidence = 0.0;
      double aziYPos = foundation.orientation;
      double aziXPos = PI/2 + foundation.orientation;
      double aziYNeg = PI + foundation.orientation;
      double aziXNeg = 3*PI/2 + foundation.orientation;

      double tilt = foundation.surfaces[s].tilt;
      if (foundation.surfaces[s].orientation == Surface::Z_POS)
      {
        incidence = cos(PI/2 - alt);
      }
      else if (foundation.surfaces[s].orientation == Surface::Z_NEG)
      {
        incidence = cos(PI/2 - alt - PI);
      }
      else
      {
        if (foundation.numberOfDimensions == 2)
        {
          // incidence is the average incidence on the exterior of a vertical cylinder
          // 2*(int(cos(alt)*cos(x),x,0,PI/2))/(2*PI)
          // 2*(integral of incidence over a quarter of the cylinder) = lit portion
          // divide by the total radians in the circle (2*PI)
          // = 2*(cos(alt))/(2*PI)
          // = cos(alt)/PI
          incidence = cos(alt)/PI;
        }
        else
        {
          double aziSurf;
          if (foundation.surfaces[s].orientation == Surface::Y_POS)
          {
            aziSurf = aziYPos;
          }
          else if (foundation.surfaces[s].orientation == Surface::X_POS)
          {
            aziSurf = aziXPos;
          }
          else if (foundation.surfaces[s].orientation == Surface::Y_NEG)
          {
            aziSurf = aziYNeg;
          }
          else //if (foundation.surfaces[s].orientation == Surface::X_NEG)
          {
            aziSurf = aziXNeg;
          }

          if (foundation.numberOfDimensions == 3 && !foundation.useSymmetry)
          {
            // incidence = cos(alt)*cos(azi-aziSurf)*sin(tilt)+sin(alt)*cos(tilt)
            // simplifies for tilt = PI/2 to = cos(alt)*cos(azi-aziSurf)
            incidence = cos(alt)*cos(azi-aziSurf);
          }
          else // if (foundation.coordinateSystem == Foundation::CS_3D_SYMMETRY)
          {
            // if symmetric, use average incidence (one side will be facing the sun,
            // the other won't).
            if (foundation.surfaces[s].orientation == Surface::Y_POS ||
                foundation.surfaces[s].orientation == Surface::Y_NEG)
            {
              if (foundation.isXSymm)
              {
                double incidenceYPos = cos(alt)*cos(azi-aziYPos);
                if (incidenceYPos < 0)
                  incidenceYPos = 0;

                double incidenceYNeg = cos(alt)*cos(azi-aziYNeg);
                if (incidenceYNeg < 0)
                  incidenceYNeg = 0;

                incidence = (incidenceYPos + incidenceYNeg)/2.0;

              }
              else
              {
                incidence = cos(alt)*cos(azi-aziSurf);
              }
            }

            if (foundation.surfaces[s].orientation == Surface::X_POS ||
                foundation.surfaces[s].orientation == Surface::X_NEG)
            {
              if (foundation.isYSymm)
              {
                double incidenceXPos = cos(alt)*cos(azi-aziXPos);
                if (incidenceXPos < 0)
                  incidenceXPos = 0;

                double incidenceXNeg = cos(alt)*cos(azi-aziXNeg);
                if (incidenceXNeg < 0)
                  incidenceXNeg = 0;

                incidence = (incidenceXPos + incidenceXNeg)/2.0;

              }
              else
              {
                incidence = cos(alt)*cos(azi-aziSurf);
              }
            }
          }
        }
      }

      // if sun is below horizon, incidence is zero
      if (sin(alt) < 0)
        incidence = 0;
      if (incidence < 0)
        incidence = 0;

      double Fsky = (1.0 + cos(tilt))/2.0;
      double Fg = 1.0 - Fsky;
      double rho_g = 1.0 - foundation.soilAbsorptivity;

      for (auto index: foundation.surfaces[s].indices)
      {
        Cell* this_cell = &domain.cell[index];
        double alpha = this_cell->surfacePtr->absorptivity;

        if (qGH > 0.0)
        {
          pssf = incidence;
          q = alpha*(qDN*pssf + qDH*Fsky + qGH*Fg*rho_g);
        }
        else
        {
          q = 0;
        }
        this_cell->heatGain = q;

      }
    }
  }
}

void Ground::setInteriorRadiationBoundaryConditions()
{
  for (std::size_t s = 0; s < foundation.surfaces.size() ; s++)
  {
    if (foundation.surfaces[s].type == Surface::ST_SLAB_CORE
        || foundation.surfaces[s].type == Surface::ST_SLAB_PERIM
        || foundation.surfaces[s].type == Surface::ST_WALL_INT)
    {
      for (auto index: foundation.surfaces[s].indices)
      {
        Cell* this_cell = &domain.cell[index];
        if (foundation.surfaces[s].type == Surface::ST_WALL_INT) {
          this_cell->heatGain = bcs.wallAbsRadiation;
        }
        else {
          this_cell->heatGain = bcs.slabAbsRadiation;
        }
      }
    }
  }
}

void Ground::link_cells_to_temp()
{
  for (Cell & this_cell : domain.cell) {
    this_cell.told_ptr = &TOld[this_cell.index];
  }
}


double getArrayValue(std::vector<std::vector<std::vector<double>>> Mat, std::size_t i, std::size_t j, std::size_t k)
{
  return Mat[i][j][k];
}

}

#endif
