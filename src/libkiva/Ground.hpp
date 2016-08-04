/* Copyright (c) 2012-2016 Big Ladder Software. All rights reserved.
* See the LICENSE file for additional terms and conditions. */

#ifndef GROUND_HPP_
#define GROUND_HPP_

#include "Mesher.hpp"
#include "Domain.hpp"
#include "BoundaryConditions.hpp"
#include "Foundation.hpp"
#include "GroundOutput.hpp"
#include "Algorithms.hpp"

#include <cmath>
#include <vector>
#include <string>

#include <boost/lexical_cast.hpp>
#include <boost/multi_array.hpp>

#include "lis.h"

class Ground
{
public:

  // Constructor
  Ground(Foundation &foundation);
  Ground(Foundation &foundation, GroundOutput::OutputMap &outputMap);

  // Destructor
  virtual ~Ground();

  Domain domain;
  Foundation &foundation;
  GroundOutput groundOutput;

  size_t nX, nY, nZ;

  boost::multi_array<double, 3> TNew; // solution, n+1
  boost::multi_array<double, 3> TOld; // solution, n

  void buildDomain();

  void calculateBoundaryLayer();
  void setNewBoundaryGeometry();
  void calculate(BoundaryConditions& boundaryConidtions, double ts=0.0);

  std::vector<double> calculateHeatFlux(const size_t &i, const size_t &j, const size_t &k);


  void calculateSurfaceAverages();
  double getSurfaceAverageValue(std::pair<Surface::SurfaceType, GroundOutput::OutputType> output);


private:

  double timestep;

  BoundaryConditions bcs;
  // Data structures

  // ADE
  boost::multi_array<double, 3> U; // ADE upper sweep, n+1
  boost::multi_array<double, 3> UOld; // ADE upper sweep, n
  boost::multi_array<double, 3> V; // ADE lower sweep, n+1
  boost::multi_array<double, 3> VOld; // ADE lower sweep, n

  // ADI
  std::vector<double> a1; // lower diagonal
  std::vector<double> a2; // main diagonal
  std::vector<double> a3; // upper diagonal
  std::vector<double> b_; // right-hand side
  std::vector<double> x_; // solution

  // Implicit
  LIS_MATRIX Amat;
  LIS_VECTOR b, x;

  LIS_SOLVER solver;

  std::vector<char> solverOptions;

private:

  // Calculators (Called from main calculator)
  void calculateADE();

  void calculateADEUpwardSweep();

  void calculateADEDownwardSweep();

  void calculateExplicit();

  void calculateMatrix(Foundation::NumericalScheme scheme);

  void calculateADI(int dim);

  // Misc. Functions
  void setAmatValue(const int i, const int j, const double val);
  void setbValue(const int i, const double val);
  void solveLinearSystem();
  void clearAmat();
  double getxValue(const int i);

  double getConvectionCoeff(double Tsurf,
                double Tamb,
                double Vair,
                double roughness,
                bool isExterior,
                double tilt);


  double getSurfaceArea(Surface::SurfaceType surfaceType);

  void setSolarBoundaryConditions();

  std::vector<std::pair<double,double>> boundaryLayer;

  double getBoundaryValue(double dist);

  double getBoundaryDistance(double val);

};



#endif /* GROUND_HPP_ */
