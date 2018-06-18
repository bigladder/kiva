/* Copyright (c) 2012-2017 Big Ladder Software LLC. All rights reserved.
* See the LICENSE file for additional terms and conditions. */

#include "fixtures/bestest-fixture.hpp"
#include "fixtures/typical-fixture.hpp"

#include "Errors.hpp"

using namespace Kiva;


class CellFixture : public BESTESTFixture {
protected:
  void SetUp() {
    specifySystem();
    ground = std::make_shared<Ground>(fnd,outputMap);
    ground->foundation.createMeshData();
    Domain domain = ground->domain;
    domain.setDomain(ground->foundation);
    cell_vector = domain.cell;
  };

  std::vector< std::shared_ptr<Cell> > cell_vector;
};

TEST_F( CellFixture, cell_basics)
{
  EXPECT_EQ(cell_vector[0]->i, 0);
  EXPECT_EQ(cell_vector[0]->j, 0);
  EXPECT_EQ(cell_vector[0]->k, 0);
  EXPECT_EQ(cell_vector[0]->index, 0);
  EXPECT_EQ(cell_vector[0]->cellType, CellType::BOUNDARY);
  EXPECT_EQ(cell_vector[0]->surfacePtr->type, Surface::SurfaceType::ST_DEEP_GROUND);

  EXPECT_EQ(cell_vector[120]->i, 38);
  EXPECT_EQ(cell_vector[120]->j, 0);
  EXPECT_EQ(cell_vector[120]->k, 2);
  EXPECT_EQ(cell_vector[120]->index, 120);
  EXPECT_EQ(cell_vector[120]->cellType, CellType::NORMAL);
//  TODO: replace surfacePtr and blockPtr NULL with nullptr. This EXPECT_EQ errors.
//  EXPECT_EQ(cell_vector[120]->surfacePtr, NULL);
}

TEST_F( GC10aADIFixture, calcCellADI)
{
  double A{0.0}, Ap{0.0}, Am{0.0}, bVal{0.0};

  auto this_cell = ground->domain.cell[0];
  this_cell->calcCellADI(1, fnd, 3600.0, bcs, Am, A, Ap, bVal);
  EXPECT_DOUBLE_EQ(A, 1);
  EXPECT_DOUBLE_EQ(Ap, 0);
  EXPECT_DOUBLE_EQ(Am, 0);
  EXPECT_DOUBLE_EQ(bVal, this_cell->surfacePtr->temperature);

  this_cell = ground->domain.cell[120];
  this_cell->calcCellADI(1, fnd, 3600.0, bcs, Am, A, Ap, bVal);
  double theta = 3600.0 / (fnd.numberOfDimensions
                             *this_cell->density*this_cell->specificHeat);
  double f = fnd.fADI;
  EXPECT_DOUBLE_EQ(A, 1.0 + (2 - f)*(this_cell->cxp - this_cell->cxm)*theta);
  EXPECT_DOUBLE_EQ(Ap, (2 - f)*(-this_cell->cxp*theta));
  EXPECT_DOUBLE_EQ(Am, (2 - f)*(this_cell->cxm*theta));
  EXPECT_DOUBLE_EQ(bVal, *this_cell->told_ptr*(1.0 + f*(this_cell->czm - this_cell->czp)*theta)
                       - *(this_cell->told_ptr - ground->domain.stepsize[2])*f*this_cell->czm*theta
                       + *(this_cell->told_ptr + ground->domain.stepsize[2])*f*this_cell->czp*theta
                       + this_cell->heatGain*theta);
}

TEST_F( GC10aImplicitFixture, calcCellMatrix)
{
  double A{0}, Aip{0}, Aim{0}, Ajp{0}, Ajm{0}, Akp{0}, Akm{0}, bVal{0};

  auto this_cell = ground->domain.cell[0];
  this_cell->calcCellMatrix(fnd.numericalScheme, 3600.0, fnd, bcs, A, Aip, Aim, Ajp, Ajm, Akp, Akm, bVal);
  EXPECT_DOUBLE_EQ(A, 1);
  EXPECT_DOUBLE_EQ(Aip, 0);
  EXPECT_DOUBLE_EQ(Aim, 0);
  EXPECT_DOUBLE_EQ(bVal, this_cell->surfacePtr->temperature);

  this_cell = ground->domain.cell[120];
  this_cell->calcCellMatrix(fnd.numericalScheme, 3600.0, fnd, bcs, A, Aip, Aim, Ajp, Ajm, Akp, Akm, bVal);
//  TODO: upgrade from regression tests to physical-ish unittests.
  double theta = 3600.0 / (this_cell->density*this_cell->specificHeat);
  EXPECT_DOUBLE_EQ(A, (1.0 + (this_cell->cxp + this_cell->czp -
                              this_cell->cxm - this_cell->czm)*theta));
  EXPECT_DOUBLE_EQ(Aip, -this_cell->cxp*theta);
  EXPECT_DOUBLE_EQ(Aim, this_cell->cxm*theta);
  EXPECT_DOUBLE_EQ(bVal, *this_cell->told_ptr + this_cell->heatGain*theta);
}

TEST_F( GC10aSteadyStateFixture, calcCellMatrixSS)
{
  double A{0}, Aip{0}, Aim{0}, Ajp{0}, Ajm{0}, Akp{0}, Akm{0}, bVal{0};
  auto this_cell = ground->domain.cell[0];
  this_cell->calcCellMatrix(fnd.numericalScheme, 3600.0, fnd, bcs, A, Aip, Aim, Ajp, Ajm, Akp, Akm, bVal);
  EXPECT_DOUBLE_EQ(A, 1);
  EXPECT_DOUBLE_EQ(Aip, 0);
  EXPECT_DOUBLE_EQ(Aim, 0);
  EXPECT_DOUBLE_EQ(bVal, this_cell->surfacePtr->temperature);

  this_cell = ground->domain.cell[120];
  this_cell->calcCellMatrix(fnd.numericalScheme, 3600.0, fnd, bcs, A, Aip, Aim, Ajp, Ajm, Akp, Akm, bVal);
  EXPECT_DOUBLE_EQ(A, this_cell->cxm + this_cell->czm - this_cell->cxp - this_cell->czp);
  EXPECT_DOUBLE_EQ(Aip, this_cell->cxp);
  EXPECT_DOUBLE_EQ(Aim, -this_cell->cxm);
  EXPECT_DOUBLE_EQ(bVal, 0);
}
