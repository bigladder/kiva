/* Ground.c++ is part of Kiva (Written by Neal Kruis)
 * Copyright (C) 2012-2013 Big Ladder Software <info@bigladdersoftware.com>
 * All rights reserved.
 *
 * Kiva is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Kiva is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Kiva.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef Ground_CPP
#define Ground_CPP

//#define USE_LIS_SOLVER

#include "Ground.h"

static const double PI = 4.0*atan(1.0);

static const bool TDMA = true;

Ground::Ground(WeatherData &weatherData, Foundation &foundation,
               SimulationControl &simulationControl, std::string outputFileName) :
               foundation(foundation), simulationControl(simulationControl),
               weatherData(weatherData)
{
  // set up output file
  outputFile.open(outputFileName.c_str());
  outputFile << "Time Stamp" << printOutputHeaders() << std::endl;

  // Build Domain Object
  buildDomain();

  // Initial Conditions
  initializeConditions();

  initializePlots();
}

Ground::~Ground()
{

  outputFile.close();

#if defined(USE_LIS_SOLVER)
  lis_matrix_destroy(Amat);
  lis_vector_destroy(x);
  lis_vector_destroy(b);
  lis_solver_destroy(solver); // for whatever reason, this causes a crash
#endif

}

void Ground::buildDomain()
{
  timestep = simulationControl.timestep.total_seconds();

  // Create mesh
  if (foundation.deepGroundBoundary == Foundation::DGB_AUTO)
    foundation.deepGroundTemperature = weatherData.dryBulbTemp.getAverage();

  std::cout << "Creating Domain..." << std::endl;
  foundation.createMeshData();

  // Build matrices for PDE term coefficients
  domain.setDomain(foundation);

  nX = domain.meshX.centers.size();
  nY = domain.meshY.centers.size();
  nZ = domain.meshZ.centers.size();

  std::cout << "  X Cells: " << nX << std::endl;
  std::cout << "  Y Cells: " << nY << std::endl;
  std::cout << "  Z Cells: " << nZ << std::endl;
  std::cout << "  Total Cells: " << nX*nY*nZ << std::endl;

  //domain.printCellTypes();
}

void Ground::initializeConditions()
{

  annualAverageDryBulbTemperature = weatherData.dryBulbTemp.getAverage();

  // Initialize matices
  if (foundation.numericalScheme == Foundation::NS_ADE)
  {
    U.resize(boost::extents[nX][nY][nZ]);
    UOld.resize(boost::extents[nX][nY][nZ]);

    V.resize(boost::extents[nX][nY][nZ]);
    VOld.resize(boost::extents[nX][nY][nZ]);
  }

  if (foundation.numericalScheme == Foundation::NS_ADI && TDMA)
  {
    a1.resize(nX*nY*nZ, 0.0);
    a2.resize(nX*nY*nZ, 0.0);
    a3.resize(nX*nY*nZ, 0.0);
    b_.resize(nX*nY*nZ, 0.0);
    x_.resize(nX*nY*nZ);
  }

#if defined(USE_LIS_SOLVER)
    std::string solverOptionsString = "-i ";
    solverOptionsString.append(foundation.solver);
    solverOptionsString.append(" -p ");
    solverOptionsString.append(foundation.preconditioner);
    solverOptionsString.append(" -maxiter ");
    solverOptionsString.append(std::to_string(foundation.maxIterations));
    solverOptionsString.append(" -initx_zeros false -tol ");
    solverOptionsString.append(std::to_string(foundation.tolerance));

    std::vector<char> solverChars(solverOptionsString.begin(),solverOptionsString.end());
    solverChars.push_back('\0');
    solverOptions = solverChars;
#endif

  if (foundation.numericalScheme == Foundation::NS_CRANK_NICOLSON ||
      foundation.numericalScheme == Foundation::NS_IMPLICIT ||
      foundation.numericalScheme == Foundation::NS_STEADY_STATE ||
      foundation.initializationMethod == Foundation::IM_STEADY_STATE ||
      foundation.implicitAccelPeriods > 0 ||
      (foundation.numericalScheme == Foundation::NS_ADI && !(TDMA)))
  {

#if defined(USE_LIS_SOLVER)
    lis_matrix_create(LIS_COMM_WORLD,&Amat);
    lis_matrix_set_size(Amat,nX*nY*nZ,nX*nY*nZ);

    lis_vector_create(LIS_COMM_WORLD,&b);
    lis_vector_set_size(b,0,nX*nY*nZ);

    lis_vector_duplicate(b,&x);

    lis_vector_set_all(283.15,x);
    lis_solver_create(&solver);
    lis_solver_set_option(&solverOptions[0],solver);
#else
        Amat = boost::numeric::ublas::compressed_matrix<double,
        boost::numeric::ublas::column_major, 0,
        boost::numeric::ublas::unbounded_array<int>,
        boost::numeric::ublas::unbounded_array<double> >(nX*nY*nZ,nX*nY*nZ);  // Coefficient Matrix
        b = boost::numeric::ublas::vector<double> (nX*nY*nZ);  // constant and unknown vectors
        x = boost::numeric::ublas::vector<double> (nX*nY*nZ);  // constant and unknown vectors
#endif

  }
  else
  {
    // Create a place holder matrix/vector system of size 1 to avoid problems when destroying the LIS components
#if defined(USE_LIS_SOLVER)
    lis_matrix_create(LIS_COMM_WORLD,&Amat);
    lis_matrix_set_size(Amat,1,1);

    lis_vector_create(LIS_COMM_WORLD,&b);
    lis_vector_set_size(b,0,1);

    lis_vector_duplicate(b,&x);

    lis_vector_set_all(283.15,x);
    lis_solver_create(&solver);
    lis_solver_set_option(&solverOptions[0],solver);
#endif
  }

  TNew.resize(boost::extents[nX][nY][nZ]);
  TOld.resize(boost::extents[nX][nY][nZ]);

  tNow = 0.0;

  prevStatusUpdate = boost::posix_time::second_clock::local_time();

  initPeriod = true;

  if (foundation.numericalScheme != Foundation::NS_STEADY_STATE)
  {
    std::cout << "Initializing Temperatures..." << std::endl;

    // Calculate initial time in seconds (simulation start minus warmup and acceleration periods)
    double simulationTimestep = simulationControl.timestep.total_seconds();
    double accelTimestep = foundation.implicitAccelTimestep*3600;

    double accelDuration = accelTimestep*foundation.implicitAccelPeriods;
    double warmupDuration = foundation.warmupDays*24*3600;

    double tInit = -warmupDuration - simulationTimestep - accelDuration - accelTimestep;

    // Calculate initial conditions
    if (foundation.initializationMethod == Foundation::IM_STEADY_STATE)
    {
      Foundation::NumericalScheme tempNS = foundation.numericalScheme;
      foundation.numericalScheme = Foundation::NS_STEADY_STATE;
      calculate(tInit);
      foundation.numericalScheme = tempNS;
    }
    else
    {
      tNow = tInit;
      for (size_t k = 0; k < nZ; ++k)
      {
        for (size_t j = 0; j < nY; ++j)
        {
          for (size_t i = 0; i < nX; ++i)
          {
            TOld[i][j][k]= getInitialTemperature(tInit,
                domain.meshZ.centers[k]);
          }
        }
      }
    }

    // Calculate implicit acceleration
    if (foundation.implicitAccelPeriods > 0)
    {
      boost::posix_time::time_duration tempTimestep = simulationControl.timestep;
      simulationControl.timestep = boost::posix_time::hours(foundation.implicitAccelTimestep);
      timestep = accelTimestep;

      double tAccelStart = -warmupDuration - simulationTimestep - accelDuration; // [s] Acceleration start time
      double tAccelEnd = -warmupDuration - simulationTimestep; // [s] Acceleration end time

      Foundation::NumericalScheme tempNS = foundation.numericalScheme;
      foundation.numericalScheme = Foundation::NS_IMPLICIT;

      for (double t = tAccelStart; t <= tAccelEnd; t += timestep)
      {
        calculate(t);
      }

      foundation.numericalScheme = tempNS;
      simulationControl.timestep = tempTimestep;
      timestep = simulationControl.timestep.total_seconds();

    }

    // Calculate warmup
    if (foundation.warmupDays > 0)
    {

      double tWarmupStart = -warmupDuration; // [s] Acceleration start time
      double tWarmupEnd = -simulationTimestep; // [s] Simulation end time

      for (double t = tWarmupStart; t <= tWarmupEnd; t += timestep)
      {
        calculate(t);
      }

    }

  }
  initPeriod = false;

}

void Ground::initializePlots()
{
  for (std::size_t p = 0; p < foundation.outputAnimations.size(); p++)
  {
    if (!foundation.outputAnimations[p].startDateSet)
      foundation.outputAnimations[p].startDate = simulationControl.startDate;

    if (!foundation.outputAnimations[p].endDateSet)
      foundation.outputAnimations[p].endDate = simulationControl.endDate;

    if (foundation.coordinateSystem == Foundation::CS_3D ||
      foundation.coordinateSystem == Foundation::CS_3D_SYMMETRY)
    {
      if (!foundation.outputAnimations[p].xRangeSet)
      {
        foundation.outputAnimations[p].xRange.first = domain.meshX.dividers[0];
        foundation.outputAnimations[p].xRange.second = domain.meshX.dividers[nX];
      }

      if (!foundation.outputAnimations[p].yRangeSet)
      {
        foundation.outputAnimations[p].yRange.first = domain.meshY.dividers[0];
        foundation.outputAnimations[p].yRange.second = domain.meshY.dividers[nY];
      }

      if (!foundation.outputAnimations[p].zRangeSet)
      {

        if (!foundation.outputAnimations[p].xRangeSet && !foundation.outputAnimations[p].yRangeSet)
        {
          foundation.outputAnimations[p].zRange.first = 0;
          foundation.outputAnimations[p].zRange.second = 0;
        }
        else
        {
          foundation.outputAnimations[p].zRange.first = domain.meshZ.dividers[0];
          foundation.outputAnimations[p].zRange.second = domain.meshZ.dividers[nZ];
        }
      }


    }
    else
    {
      if (!foundation.outputAnimations[p].xRangeSet)
      {
        foundation.outputAnimations[p].xRange.first = domain.meshX.dividers[0];
        foundation.outputAnimations[p].xRange.second = domain.meshX.dividers[nX];
      }

      if (!foundation.outputAnimations[p].yRangeSet)
      {
        foundation.outputAnimations[p].yRange.first = 0.5;
        foundation.outputAnimations[p].yRange.second = 0.5;
      }

      if (!foundation.outputAnimations[p].zRangeSet)
      {
        foundation.outputAnimations[p].zRange.first = domain.meshZ.dividers[0];
        foundation.outputAnimations[p].zRange.second = domain.meshZ.dividers[nZ];
      }
    }

    plots.push_back(GroundPlot(foundation.outputAnimations[p],domain,foundation.blocks));

    boost::posix_time::ptime simStart = simulationControl.startTime;
    boost::posix_time::ptime startTime(foundation.outputAnimations[p].startDate,boost::posix_time::hours(0));;
    boost::posix_time::ptime endTime(foundation.outputAnimations[p].endDate + boost::gregorian::days(1));
    boost::posix_time::time_duration untilStart =  startTime - simStart;
    boost::posix_time::time_duration untilEnd =  endTime - simStart;

    plots[p].tStart = untilStart.total_seconds();
    plots[p].nextPlotTime = untilStart.total_seconds();
    plots[p].tEnd = untilEnd.total_seconds();
  }
}

void Ground::simulate()
{
  std::cout << "Beginning Simulation..." << std::endl;

  boost::posix_time::ptime simStart = simulationControl.startTime;
  boost::posix_time::ptime simEnd(simulationControl.endDate + boost::gregorian::days(1));
  boost::posix_time::time_duration simDuration =  simEnd - simStart;


  double tstart = 0.0; // [s] Simulation start time
  prevOutputTime = -foundation.outputReport.minFrequency.total_seconds();
  double tend = simDuration.total_seconds(); // [s] Simulation end time

  for (double t = tstart; t < tend; t = t + timestep)
  {

    percentComplete = round(t/tend*1000)/10.0;
    calculate(t);

    if (t - prevOutputTime >= foundation.outputReport.minFrequency.total_seconds())
    {
      outputFile << to_simple_string(getSimTime(t)) << printOutputLine() << std::endl;
      prevOutputTime = t;
    }

  }

  std::cout << "  " << getSimTime(tend - timestep) << " (100%)" << std::endl;

}

void Ground::calculateADE()
{
  // Set Old values
  for (size_t k = 0; k < nZ; ++k)
  {
    for (size_t j = 0; j < nY; ++j)
    {
      for (size_t i = 0; i < nX; ++i)
      {
        UOld[i][j][k] = VOld[i][j][k] = TOld[i][j][k];
      }

    }
  }

  // Solve for new values (Main loop)
  #pragma omp parallel sections num_threads(2)
  {
    #pragma omp section
      calculateADEUpwardSweep();
    #pragma omp section
      calculateADEDownwardSweep();
  }

  for (size_t k = 0; k < nZ; ++k)
  {
    for (size_t j = 0; j < nY; ++j)
    {
      for (size_t i = 0; i < nX; ++i)
      {
        // Calculate average of sweeps
        TNew[i][j][k] = 0.5*(U[i][j][k] + V[i][j][k]);

        // Update old values for next timestep
        TOld[i][j][k] = TNew[i][j][k];
      }
    }
  }
}

void Ground::calculateADEUpwardSweep()
{
  // Upward sweep (Solve U Matrix starting from 1, 1)
  for (size_t k = 0; k < nZ; k++)
  {
    for (size_t j = 0; j < nY; ++j)
    {
      for (size_t i = 0; i < nX; i++)
      {
        switch (domain.cell[i][j][k].cellType)
        {
        case Cell::BOUNDARY:
          {
          double tilt;
          if (domain.cell[i][j][k].surface.orientation == Surface::Z_POS)
            tilt = 0;
          else if (domain.cell[i][j][k].surface.orientation == Surface::Z_NEG)
            tilt = PI;
          else
            tilt = PI/2.0;

          switch (domain.cell[i][j][k].surface.boundaryConditionType)
          {
          case Surface::ZERO_FLUX:
            {
            switch (domain.cell[i][j][k].surface.orientation)
            {
            case Surface::X_NEG:
              U[i][j][k] = UOld[i+1][j][k];
              break;
            case Surface::X_POS:
              U[i][j][k] = U[i-1][j][k];
              break;
            case Surface::Y_NEG:
              U[i][j][k] = UOld[i][j+1][k];
              break;
            case Surface::Y_POS:
              U[i][j][k] = U[i][j-1][k];
              break;
            case Surface::Z_NEG:
              U[i][j][k] = UOld[i][j][k+1];
              break;
            case Surface::Z_POS:
              U[i][j][k] = U[i][j][k-1];
              break;
            }
            }
            break;

          case Surface::CONSTANT_TEMPERATURE:

            U[i][j][k] = domain.cell[i][j][k].surface.temperature;
            break;

          case Surface::INTERIOR_TEMPERATURE:

            U[i][j][k] = foundation.indoorAirTemperature;
            break;

          case Surface::EXTERIOR_TEMPERATURE:

            U[i][j][k] = getOutdoorTemperature();
            break;

          case Surface::INTERIOR_FLUX:
            {
            double Tair = foundation.indoorAirTemperature;
            double q = domain.cell[i][j][k].heatGain;

            double hc = getConvectionCoeff(TOld[i][j][k],
                    Tair,0.0,1.52,false,tilt);
            double hr = getSimpleInteriorIRCoeff(domain.cell[i][j][k].surface.emissivity,
                               TOld[i][j][k],Tair);

            switch (domain.cell[i][j][k].surface.orientation)
            {
            case Surface::X_NEG:
              U[i][j][k] = (domain.getKXP(i,j,k)*UOld[i+1][j][k]/domain.getDXP(i) +
                   (hc + hr)*Tair + q)/(domain.getKXP(i,j,k)/domain.getDXP(i) + (hc + hr));
              break;
            case Surface::X_POS:
              U[i][j][k] = (domain.getKXM(i,j,k)*U[i-1][j][k]/domain.getDXM(i) +
                   (hc + hr)*Tair + q)/(domain.getKXM(i,j,k)/domain.getDXM(i) + (hc + hr));
              break;
            case Surface::Y_NEG:
              U[i][j][k] = (domain.getKYP(i,j,k)*UOld[i][j+1][k]/domain.getDYP(j) +
                   (hc + hr)*Tair + q)/(domain.getKYP(i,j,k)/domain.getDYP(j) + (hc + hr));
              break;
            case Surface::Y_POS:
              U[i][j][k] = (domain.getKYM(i,j,k)*U[i][j-1][k]/domain.getDYM(j) +
                   (hc + hr)*Tair + q)/(domain.getKYM(i,j,k)/domain.getDYM(j) + (hc + hr));
              break;
            case Surface::Z_NEG:
              U[i][j][k] = (domain.getKZP(i,j,k)*UOld[i][j][k+1]/domain.getDZP(k) +
                   (hc + hr)*Tair + q)/(domain.getKZP(i,j,k)/domain.getDZP(k) + (hc + hr));
              break;
            case Surface::Z_POS:
              U[i][j][k] = (domain.getKZM(i,j,k)*U[i][j][k-1]/domain.getDZM(k) +
                   (hc + hr)*Tair + q)/(domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr));
              break;
            }
            }
            break;

          case Surface::EXTERIOR_FLUX:
            {
            double Tair = getOutdoorTemperature();
            double v = getLocalWindSpeed();
            double eSky = weatherData.skyEmissivity.getValue(getSimTime(tNow));
            double F = getEffectiveExteriorViewFactor(eSky,tilt);
            double hc = getConvectionCoeff(TOld[i][j][k],Tair,v,foundation.surfaceRoughness,true,tilt);
            double hr = getExteriorIRCoeff(domain.cell[i][j][k].surface.emissivity,TOld[i][j][k],Tair,eSky,tilt);
            double q = domain.cell[i][j][k].heatGain;

            switch (domain.cell[i][j][k].surface.orientation)
            {
            case Surface::X_NEG:
              U[i][j][k] = (domain.getKXP(i,j,k)*UOld[i+1][j][k]/domain.getDXP(i) +
                    (hc + hr*pow(F,0.25))*Tair + q)/(domain.getKXP(i,j,k)/domain.getDXP(i) + (hc + hr));
              break;
            case Surface::X_POS:
              U[i][j][k] = (domain.getKXM(i,j,k)*U[i-1][j][k]/domain.getDXM(i) +
                    (hc + hr*pow(F,0.25))*Tair + q)/(domain.getKXM(i,j,k)/domain.getDXM(i) + (hc + hr));
              break;
            case Surface::Y_NEG:
              U[i][j][k] = (domain.getKYP(i,j,k)*UOld[i][j+1][k]/domain.getDYP(j) +
                    (hc + hr*pow(F,0.25))*Tair + q)/(domain.getKYP(i,j,k)/domain.getDYP(j) + (hc + hr));
              break;
            case Surface::Y_POS:
              U[i][j][k] = (domain.getKYM(i,j,k)*U[i][j-1][k]/domain.getDYM(j) +
                    (hc + hr*pow(F,0.25))*Tair + q)/(domain.getKYM(i,j,k)/domain.getDYM(j) + (hc + hr));
              break;
            case Surface::Z_NEG:
              U[i][j][k] = (domain.getKZP(i,j,k)*UOld[i][j][k+1]/domain.getDZP(k) +
                    (hc + hr*pow(F,0.25))*Tair + q)/(domain.getKZP(i,j,k)/domain.getDZP(k) + (hc + hr));
              break;
            case Surface::Z_POS:
              U[i][j][k] = (domain.getKZM(i,j,k)*U[i][j][k-1]/domain.getDZM(k) +
                    (hc + hr*pow(F,0.25))*Tair + q)/(domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr));
              break;
            }
            }
            break;
          }
          }
          break;
        case Cell::INTERIOR_AIR:
          U[i][j][k] = foundation.indoorAirTemperature;
          break;
        case Cell::EXTERIOR_AIR:
          U[i][j][k] = getOutdoorTemperature();
          break;
        default:
          {
          double theta = timestep/
            (domain.cell[i][j][k].density*domain.cell[i][j][k].specificHeat);

          double CXP = domain.cell[i][j][k].cxp*theta;
          double CXM = domain.cell[i][j][k].cxm*theta;
          double CZP = domain.cell[i][j][k].czp*theta;
          double CZM = domain.cell[i][j][k].czm*theta;
          double CYP = domain.cell[i][j][k].cyp*theta;
          double CYM = domain.cell[i][j][k].cym*theta;
          double Q = domain.cell[i][j][k].heatGain*theta;

          if (foundation.coordinateSystem == Foundation::CS_3D ||
            foundation.coordinateSystem == Foundation::CS_3D_SYMMETRY)
            U[i][j][k] = (UOld[i][j][k]*(1.0 - CXP - CZP - CYP)
                - U[i-1][j][k]*CXM
                + UOld[i+1][j][k]*CXP
                - U[i][j][k-1]*CZM
                + UOld[i][j][k+1]*CZP
                - U[i][j-1][k]*CYM
                + UOld[i][j+1][k]*CYP
                + Q) /
                (1.0 - CXM - CZM - CYM);
          else
          {
            double CXPC = 0;
            double CXMC = 0;

            if (i != 0)
            {
              double r = domain.meshX.centers[i];
              CXPC = domain.cell[i][j][k].cxp_c*theta/r;
              CXMC = domain.cell[i][j][k].cxm_c*theta/r;
            }
            U[i][j][k] = (UOld[i][j][k]*(1.0 - CXPC - CXP - CZP)
                - U[i-1][j][k]*(CXMC + CXM)
                + UOld[i+1][j][k]*(CXPC + CXP)
                - U[i][j][k-1]*CZM
                + UOld[i][j][k+1]*CZP
                + Q) /
                (1.0 - CXMC - CXM - CZM);
          }
          }
          break;
        }
      }
    }
  }
}

void Ground::calculateADEDownwardSweep()
{
  // Downward sweep (Solve V Matrix starting from I, K)
  for (size_t k = nZ - 1; k >= 0 && k < nZ; k--)
  {
    for (size_t j = nY - 1; j >= 0 && j < nY; j--)
    {
      for (size_t i = nX - 1; i >= 0 && i < nX; i--)
      {
        switch (domain.cell[i][j][k].cellType)
        {
        case Cell::BOUNDARY:
          {
          double tilt;
          if (domain.cell[i][j][k].surface.orientation == Surface::Z_POS)
            tilt = 0;
          else if (domain.cell[i][j][k].surface.orientation == Surface::Z_NEG)
            tilt = PI;
          else
            tilt = PI/2.0;

          switch (domain.cell[i][j][k].surface.boundaryConditionType)
          {
          case Surface::ZERO_FLUX:
            {
            switch (domain.cell[i][j][k].surface.orientation)
            {
            case Surface::X_NEG:
              V[i][j][k] = V[i+1][j][k];
              break;
            case Surface::X_POS:
              V[i][j][k] = VOld[i-1][j][k];
              break;
            case Surface::Y_NEG:
              V[i][j][k] = V[i][j+1][k];
              break;
            case Surface::Y_POS:
              V[i][j][k] = VOld[i][j-1][k];
              break;
            case Surface::Z_NEG:
              V[i][j][k] = V[i][j][k+1];
              break;
            case Surface::Z_POS:
              V[i][j][k] = VOld[i][j][k-1];
              break;
            }
            }
            break;

          case Surface::CONSTANT_TEMPERATURE:

            V[i][j][k] = domain.cell[i][j][k].surface.temperature;
            break;

          case Surface::INTERIOR_TEMPERATURE:

            V[i][j][k] = foundation.indoorAirTemperature;
            break;

          case Surface::EXTERIOR_TEMPERATURE:

            V[i][j][k] = getOutdoorTemperature();
            break;

          case Surface::INTERIOR_FLUX:
            {
            double Tair = foundation.indoorAirTemperature;
            double q = 0;

            double hc = getConvectionCoeff(TOld[i][j][k],
                    Tair,0.0,1.52,false,tilt);
            double hr = getSimpleInteriorIRCoeff(domain.cell[i][j][k].surface.emissivity,
                               TOld[i][j][k],Tair);

            switch (domain.cell[i][j][k].surface.orientation)
            {
            case Surface::X_NEG:
              V[i][j][k] = (domain.getKXP(i,j,k)*V[i+1][j][k]/domain.getDXP(i) +
                    (hc + hr)*Tair + q)/(domain.getKXP(i,j,k)/domain.getDXP(i) + (hc + hr));
              break;
            case Surface::X_POS:
              V[i][j][k] = (domain.getKXM(i,j,k)*VOld[i-1][j][k]/domain.getDXM(i) +
                    (hc + hr)*Tair + q)/(domain.getKXM(i,j,k)/domain.getDXM(i) + (hc + hr));
              break;
            case Surface::Y_NEG:
              V[i][j][k] = (domain.getKYP(i,j,k)*V[i][j+1][k]/domain.getDYP(j) +
                    (hc + hr)*Tair + q)/(domain.getKYP(i,j,k)/domain.getDYP(j) + (hc + hr));
              break;
            case Surface::Y_POS:
              V[i][j][k] = (domain.getKYM(i,j,k)*VOld[i][j-1][k]/domain.getDYM(j) +
                    (hc + hr)*Tair + q)/(domain.getKYM(i,j,k)/domain.getDYM(j) + (hc + hr));
              break;
            case Surface::Z_NEG:
              V[i][j][k] = (domain.getKZP(i,j,k)*V[i][j][k+1]/domain.getDZP(k) +
                    (hc + hr)*Tair + q)/(domain.getKZP(i,j,k)/domain.getDZP(k) + (hc + hr));
              break;
            case Surface::Z_POS:
              V[i][j][k] = (domain.getKZM(i,j,k)*VOld[i][j][k-1]/domain.getDZM(k) +
                    (hc + hr)*Tair + q)/(domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr));
              break;
            }
            }
            break;

          case Surface::EXTERIOR_FLUX:
            {
            double Tair = getOutdoorTemperature();
            double v = getLocalWindSpeed();
            double eSky = weatherData.skyEmissivity.getValue(getSimTime(tNow));
            double F = getEffectiveExteriorViewFactor(eSky,tilt);
            double hc = getConvectionCoeff(TOld[i][j][k],Tair,v,foundation.surfaceRoughness,true,tilt);
            double hr = getExteriorIRCoeff(domain.cell[i][j][k].surface.emissivity,TOld[i][j][k],Tair,eSky,tilt);
            double q = domain.cell[i][j][k].surface.absorptivity*weatherData.globalHorizontalSolar.getValue(getSimTime(tNow));

            switch (domain.cell[i][j][k].surface.orientation)
            {
            case Surface::X_NEG:
              V[i][j][k] = (domain.getKXP(i,j,k)*V[i+1][j][k]/domain.getDXP(i) +
                  (hc + hr*pow(F,0.25))*Tair + q)/(domain.getKXP(i,j,k)/domain.getDXP(i) + (hc + hr));
              break;
            case Surface::X_POS:
              V[i][j][k] = (domain.getKXM(i,j,k)*VOld[i-1][j][k]/domain.getDXM(i) +
                  (hc + hr*pow(F,0.25))*Tair + q)/(domain.getKXM(i,j,k)/domain.getDXM(i) + (hc + hr));
              break;
            case Surface::Y_NEG:
              V[i][j][k] = (domain.getKYP(i,j,k)*V[i][j+1][k]/domain.getDYP(j) +
                  (hc + hr*pow(F,0.25))*Tair + q)/(domain.getKYP(i,j,k)/domain.getDYP(j) + (hc + hr));
              break;
            case Surface::Y_POS:
              V[i][j][k] = (domain.getKYM(i,j,k)*VOld[i][j-1][k]/domain.getDYM(j) +
                  (hc + hr*pow(F,0.25))*Tair + q)/(domain.getKYM(i,j,k)/domain.getDYM(j) + (hc + hr));
              break;
            case Surface::Z_NEG:
              V[i][j][k] = (domain.getKZP(i,j,k)*VOld[i][j][k+1]/domain.getDZP(k) +
                  (hc + hr*pow(F,0.25))*Tair + q)/(domain.getKZP(i,j,k)/domain.getDZP(k) + (hc + hr));
              break;
            case Surface::Z_POS:
              V[i][j][k] = (domain.getKZM(i,j,k)*VOld[i][j][k-1]/domain.getDZM(k) +
                  (hc + hr*pow(F,0.25))*Tair + q)/(domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr));
              break;
            }
            }
            break;
          }
          }
          break;

        case Cell::INTERIOR_AIR:
          V[i][j][k] = foundation.indoorAirTemperature;
          break;
        case Cell::EXTERIOR_AIR:
          V[i][j][k] = getOutdoorTemperature();
          break;
        default:
          {
          double theta = timestep/
            (domain.cell[i][j][k].density*domain.cell[i][j][k].specificHeat);

          double CXP = domain.cell[i][j][k].cxp*theta;
          double CXM = domain.cell[i][j][k].cxm*theta;
          double CZP = domain.cell[i][j][k].czp*theta;
          double CZM = domain.cell[i][j][k].czm*theta;
          double CYP = domain.cell[i][j][k].cyp*theta;
          double CYM = domain.cell[i][j][k].cym*theta;
          double Q = domain.cell[i][j][k].heatGain*theta;

          if (foundation.coordinateSystem == Foundation::CS_3D ||
            foundation.coordinateSystem == Foundation::CS_3D_SYMMETRY)
            V[i][j][k] = (VOld[i][j][k]*(1.0 + CXM + CZM + CYM)
                - VOld[i-1][j][k]*CXM
                + V[i+1][j][k]*CXP
                - VOld[i][j][k-1]*CZM
                + V[i][j][k+1]*CZP
                - VOld[i][j-1][k]*CYM
                + V[i][j+1][k]*CYP
                + Q) /
                (1.0 + CXP + CZP + CYP);
          else
          {
            double CXPC = 0;
            double CXMC = 0;

            if (i != 0)
            {
              double r = domain.meshX.centers[i];
              CXPC = domain.cell[i][j][k].cxp_c*theta/r;
              CXMC = domain.cell[i][j][k].cxm_c*theta/r;
            }
            V[i][j][k] = (VOld[i][j][k]*(1.0 + CXMC + CXM + CZM)
                - VOld[i-1][j][k]*(CXMC + CXM)
                + V[i+1][j][k]*(CXPC + CXP)
                - VOld[i][j][k-1]*CZM
                + V[i][j][k+1]*CZP
                + Q) /
                (1.0 + CXPC + CXP + CZP);
          }
          }
          break;
        }
      }
    }
  }
}

void Ground::calculateExplicit()
{
  for (size_t k = 0; k < nZ; k++)
  {
    for (size_t j = 0; j < nY; ++j)
    {
      for (size_t i = 0; i < nX; i++)
      {
        switch (domain.cell[i][j][k].cellType)
        {
        case Cell::BOUNDARY:
          {
          double tilt;
          if (domain.cell[i][j][k].surface.orientation == Surface::Z_POS)
            tilt = 0;
          else if (domain.cell[i][j][k].surface.orientation == Surface::Z_NEG)
            tilt = PI;
          else
            tilt = PI/2.0;

          switch (domain.cell[i][j][k].surface.boundaryConditionType)
          {
          case Surface::ZERO_FLUX:
            {
            switch (domain.cell[i][j][k].surface.orientation)
            {
            case Surface::X_NEG:
              TNew[i][j][k] = TOld[i+1][j][k];
              break;
            case Surface::X_POS:
              TNew[i][j][k] = TOld[i-1][j][k];
              break;
            case Surface::Y_NEG:
              TNew[i][j][k] = TOld[i][j+1][k];
              break;
            case Surface::Y_POS:
              TNew[i][j][k] = TOld[i][j-1][k];
              break;
            case Surface::Z_NEG:
              TNew[i][j][k] = TOld[i][j][k+1];
              break;
            case Surface::Z_POS:
              TNew[i][j][k] = TOld[i][j][k-1];
              break;
            }
            }
            break;

          case Surface::CONSTANT_TEMPERATURE:

            TNew[i][j][k] = domain.cell[i][j][k].surface.temperature;
            break;

          case Surface::INTERIOR_TEMPERATURE:

            TNew[i][j][k] = foundation.indoorAirTemperature;
            break;

          case Surface::EXTERIOR_TEMPERATURE:

            TNew[i][j][k] = getOutdoorTemperature();
            break;

          case Surface::INTERIOR_FLUX:
            {
            double Tair = foundation.indoorAirTemperature;
            double q = 0;

            double hc = getConvectionCoeff(TOld[i][j][k],
                    Tair,0.0,1.52,false,tilt);
            double hr = getSimpleInteriorIRCoeff(domain.cell[i][j][k].surface.emissivity,
                               TOld[i][j][k],Tair);

            switch (domain.cell[i][j][k].surface.orientation)
            {
            case Surface::X_NEG:
              TNew[i][j][k] = (domain.getKXP(i,j,k)*TOld[i+1][j][k]/domain.getDXP(i) +
                    (hc + hr)*Tair + q)/(domain.getKXP(i,j,k)/domain.getDXP(i) + (hc + hr));
              break;
            case Surface::X_POS:
              TNew[i][j][k] = (domain.getKXM(i,j,k)*TOld[i-1][j][k]/domain.getDXM(i) +
                    (hc + hr)*Tair + q)/(domain.getKXM(i,j,k)/domain.getDXM(i) + (hc + hr));
              break;
            case Surface::Y_NEG:
              TNew[i][j][k] = (domain.getKYP(i,j,k)*TOld[i][j+1][k]/domain.getDYP(j) +
                    (hc + hr)*Tair + q)/(domain.getKYP(i,j,k)/domain.getDYP(j) + (hc + hr));
              break;
            case Surface::Y_POS:
              TNew[i][j][k] = (domain.getKYM(i,j,k)*TOld[i][j-1][k]/domain.getDYM(j) +
                    (hc + hr)*Tair + q)/(domain.getKYM(i,j,k)/domain.getDYM(j) + (hc + hr));
              break;
            case Surface::Z_NEG:
              TNew[i][j][k] = (domain.getKZP(i,j,k)*TOld[i][j][k+1]/domain.getDZP(k) +
                    (hc + hr)*Tair + q)/(domain.getKZP(i,j,k)/domain.getDZP(k) + (hc + hr));
              break;
            case Surface::Z_POS:
              TNew[i][j][k] = (domain.getKZM(i,j,k)*TOld[i][j][k-1]/domain.getDZM(k) +
                    (hc + hr)*Tair + q)/(domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr));
              break;
            }
            }
            break;

          case Surface::EXTERIOR_FLUX:
            {
            double Tair = getOutdoorTemperature();
            double v = getLocalWindSpeed();
            double eSky = weatherData.skyEmissivity.getValue(getSimTime(tNow));
            double F = getEffectiveExteriorViewFactor(eSky,tilt);
            double hc = getConvectionCoeff(TOld[i][j][k],Tair,v,foundation.surfaceRoughness,true,tilt);
            double hr = getExteriorIRCoeff(domain.cell[i][j][k].surface.emissivity,TOld[i][j][k],Tair,eSky,tilt);
            double q = domain.cell[i][j][k].surface.absorptivity*weatherData.globalHorizontalSolar.getValue(getSimTime(tNow));

            switch (domain.cell[i][j][k].surface.orientation)
            {
            case Surface::X_NEG:
              TNew[i][j][k] = (domain.getKXP(i,j,k)*TOld[i+1][j][k]/domain.getDXP(i) +
                  (hc + hr*pow(F,0.25))*Tair + q)/(domain.getKXP(i,j,k)/domain.getDXP(i) + (hc + hr));
              break;
            case Surface::X_POS:
              TNew[i][j][k] = (domain.getKXM(i,j,k)*TOld[i-1][j][k]/domain.getDXM(i) +
                  (hc + hr*pow(F,0.25))*Tair + q)/(domain.getKXM(i,j,k)/domain.getDXM(i) + (hc + hr));
              break;
            case Surface::Y_NEG:
              TNew[i][j][k] = (domain.getKYP(i,j,k)*TOld[i][j+1][k]/domain.getDYP(j) +
                  (hc + hr*pow(F,0.25))*Tair + q)/(domain.getKYP(i,j,k)/domain.getDYP(j) + (hc + hr));
              break;
            case Surface::Y_POS:
              TNew[i][j][k] = (domain.getKYM(i,j,k)*TOld[i][j-1][k]/domain.getDYM(j) +
                  (hc + hr*pow(F,0.25))*Tair + q)/(domain.getKYM(i,j,k)/domain.getDYM(j) + (hc + hr));
              break;
            case Surface::Z_NEG:
              TNew[i][j][k] = (domain.getKZP(i,j,k)*TOld[i][j][k+1]/domain.getDZP(k) +
                  (hc + hr*pow(F,0.25))*Tair + q)/(domain.getKZP(i,j,k)/domain.getDZP(k) + (hc + hr));
              break;
            case Surface::Z_POS:
              TNew[i][j][k] = (domain.getKZM(i,j,k)*TOld[i][j][k-1]/domain.getDZM(k) +
                  (hc + hr*pow(F,0.25))*Tair + q)/(domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr));
              break;
            }
            }
            break;
          }
          }
          break;
        case Cell::INTERIOR_AIR:
          TNew[i][j][k] = foundation.indoorAirTemperature;
          break;

        case Cell::EXTERIOR_AIR:
          TNew[i][j][k] = getOutdoorTemperature();
          break;
        default:
          {
          double theta = timestep/
            (domain.cell[i][j][k].density*domain.cell[i][j][k].specificHeat);

          double CXP = domain.cell[i][j][k].cxp*theta;
          double CXM = domain.cell[i][j][k].cxm*theta;
          double CZP = domain.cell[i][j][k].czp*theta;
          double CZM = domain.cell[i][j][k].czm*theta;
          double CYP = domain.cell[i][j][k].cyp*theta;
          double CYM = domain.cell[i][j][k].cym*theta;
          double Q = domain.cell[i][j][k].heatGain*theta;

          if (foundation.coordinateSystem == Foundation::CS_3D ||
            foundation.coordinateSystem == Foundation::CS_3D_SYMMETRY)
            TNew[i][j][k] = TOld[i][j][k]*(1.0 + CXM + CZM + CYM - CXP - CZP - CYP)
                - TOld[i-1][j][k]*CXM
                + TOld[i+1][j][k]*CXP
                - TOld[i][j][k-1]*CZM
                + TOld[i][j][k+1]*CZP
                - TOld[i][j-1][k]*CYM
                + TOld[i][j+1][k]*CYP
                + Q;
          else
          {
            double CXPC = 0;
            double CXMC = 0;

            if (i != 0)
            {
              double r = domain.meshX.centers[i];
              CXPC = domain.cell[i][j][k].cxp_c*theta/r;
              CXMC = domain.cell[i][j][k].cxm_c*theta/r;
            }

            TNew[i][j][k] = TOld[i][j][k]*(1.0 + CXMC + CXM + CZM - CXPC - CXP - CZP)
                - TOld[i-1][j][k]*(CXMC + CXM)
                + TOld[i+1][j][k]*(CXPC + CXP)
                - TOld[i][j][k-1]*CZM
                + TOld[i][j][k+1]*CZP
                + Q;
          }
          }
          break;
        }
      }
    }
  }
  for (size_t k = 0; k < nZ; ++k)
  {
    for (size_t j = 0; j < nY; ++j)
    {
      for (size_t i = 0; i < nX; ++i)
      {
        // Update old values for next timestep
        TOld[i][j][k] = TNew[i][j][k];
      }
    }
  }
}

void Ground::calculateMatrix(Foundation::NumericalScheme scheme)
{
  for (size_t k = 0; k < nZ; k++)
  {
    for (size_t j = 0; j < nY; j++)
    {
      for (size_t i = 0; i < nX; i++)
      {
        int index = i + nX*j + nX*nY*k;
        int index_ip = (i+1) + nX*j + nX*nY*k;
        int index_im = (i-1) + nX*j + nX*nY*k;
        int index_jp = i + nX*(j+1) + nX*nY*k;
        int index_jm = i + nX*(j-1) + nX*nY*k;
        int index_kp = i + nX*j + nX*nY*(k+1);
        int index_km = i + nX*j + nX*nY*(k-1);

        double A, Aip, Aim, Ajp, Ajm, Akp, Akm, bVal = 0.0;

        switch (domain.cell[i][j][k].cellType)
        {
        case Cell::BOUNDARY:
          {
          double tilt;
          if (domain.cell[i][j][k].surface.orientation == Surface::Z_POS)
            tilt = 0;
          else if (domain.cell[i][j][k].surface.orientation == Surface::Z_NEG)
            tilt = PI;
          else
            tilt = PI/2.0;

          switch (domain.cell[i][j][k].surface.boundaryConditionType)
          {
          case Surface::ZERO_FLUX:
            {
            switch (domain.cell[i][j][k].surface.orientation)
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
            bVal = domain.cell[i][j][k].surface.temperature;

            setAmatValue(index,index,A);
            setbValue(index,bVal);
            break;
          case Surface::INTERIOR_TEMPERATURE:
            A = 1.0;
            bVal = foundation.indoorAirTemperature;

            setAmatValue(index,index,A);
            setbValue(index,bVal);
            break;
          case Surface::EXTERIOR_TEMPERATURE:
            A = 1.0;
            bVal = getOutdoorTemperature();

            setAmatValue(index,index,A);
            setbValue(index,bVal);
            break;
          case Surface::INTERIOR_FLUX:
            {
            double Tair = foundation.indoorAirTemperature;
            double q = 0;

            double hc = getConvectionCoeff(TOld[i][j][k],
                    Tair,0.0,1.52,false,tilt);
            double hr = getSimpleInteriorIRCoeff(domain.cell[i][j][k].surface.emissivity,
                               TOld[i][j][k],Tair);

            switch (domain.cell[i][j][k].surface.orientation)
            {
            case Surface::X_NEG:
              A = domain.getKXP(i,j,k)/domain.getDXP(i) + (hc + hr);
              Aip = -domain.getKXP(i,j,k)/domain.getDXP(i);
              bVal = (hc + hr)*Tair + q;

              setAmatValue(index,index,A);
              setAmatValue(index,index_ip,Aip);
              setbValue(index,bVal);
              break;
            case Surface::X_POS:
              A = domain.getKXM(i,j,k)/domain.getDXM(i) + (hc + hr);
              Aim = -domain.getKXM(i,j,k)/domain.getDXM(i);
              bVal = (hc + hr)*Tair + q;

              setAmatValue(index,index,A);
              setAmatValue(index,index_im,Aim);
              setbValue(index,bVal);
              break;
            case Surface::Y_NEG:
              A = domain.getKYP(i,j,k)/domain.getDYP(j) + (hc + hr);
              Ajp = -domain.getKYP(i,j,k)/domain.getDYP(j);
              bVal = (hc + hr)*Tair + q;

              setAmatValue(index,index,A);
              setAmatValue(index,index_jp,Ajp);
              setbValue(index,bVal);
              break;
            case Surface::Y_POS:
              A = domain.getKYM(i,j,k)/domain.getDYM(j) + (hc + hr);
              Ajm = -domain.getKYM(i,j,k)/domain.getDYM(j);
              bVal = (hc + hr)*Tair + q;

              setAmatValue(index,index,A);
              setAmatValue(index,index_jm,Ajm);
              setbValue(index,bVal);
              break;
            case Surface::Z_NEG:
              A = domain.getKZP(i,j,k)/domain.getDZP(k) + (hc + hr);
              Akp = -domain.getKZP(i,j,k)/domain.getDZP(k);
              bVal = (hc + hr)*Tair + q;

              setAmatValue(index,index,A);
              setAmatValue(index,index_kp,Akp);
              setbValue(index,bVal);
              break;
            case Surface::Z_POS:
              A = domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr);
              Akm = -domain.getKZM(i,j,k)/domain.getDZM(k);
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
            double Tair = getOutdoorTemperature();
            double v = getLocalWindSpeed();
            double eSky = weatherData.skyEmissivity.getValue(getSimTime(tNow));
            double F = getEffectiveExteriorViewFactor(eSky,tilt);
            double hc = getConvectionCoeff(TOld[i][j][k],Tair,v,foundation.surfaceRoughness,true,tilt);
            double hr = getExteriorIRCoeff(domain.cell[i][j][k].surface.emissivity,TOld[i][j][k],Tair,eSky,tilt);
            double q = domain.cell[i][j][k].surface.absorptivity*weatherData.globalHorizontalSolar.getValue(getSimTime(tNow));

            switch (domain.cell[i][j][k].surface.orientation)
            {
            case Surface::X_NEG:
              A = domain.getKXP(i,j,k)/domain.getDXP(i) + (hc + hr);
              Aip = -domain.getKXP(i,j,k)/domain.getDXP(i);
              bVal = (hc + hr*pow(F,0.25))*Tair + q;

              setAmatValue(index,index,A);
              setAmatValue(index,index_ip,Aip);
              setbValue(index,bVal);
              break;
            case Surface::X_POS:
              A = domain.getKXM(i,j,k)/domain.getDXM(i) + (hc + hr);
              Aim = -domain.getKXM(i,j,k)/domain.getDXM(i);
              bVal = (hc + hr*pow(F,0.25))*Tair + q;

              setAmatValue(index,index,A);
              setAmatValue(index,index_im,Aim);
              setbValue(index,bVal);
              break;
            case Surface::Y_NEG:
              A = domain.getKYP(i,j,k)/domain.getDYP(j) + (hc + hr);
              Ajp = -domain.getKYP(i,j,k)/domain.getDYP(j);
              bVal = (hc + hr*pow(F,0.25))*Tair + q;

              setAmatValue(index,index,A);
              setAmatValue(index,index_jp,Ajp);
              setbValue(index,bVal);
              break;
            case Surface::Y_POS:
              A = domain.getKYM(i,j,k)/domain.getDYM(j) + (hc + hr);
              Ajm = -domain.getKYM(i,j,k)/domain.getDYM(j);
              bVal = (hc + hr*pow(F,0.25))*Tair + q;

              setAmatValue(index,index,A);
              setAmatValue(index,index_jm,Ajm);
              setbValue(index,bVal);
              break;
            case Surface::Z_NEG:
              A = domain.getKZP(i,j,k)/domain.getDZP(k) + (hc + hr);
              Akp = -domain.getKZP(i,j,k)/domain.getDZP(k);
              bVal = (hc + hr*pow(F,0.25))*Tair + q;

              setAmatValue(index,index,A);
              setAmatValue(index,index_kp,Akp);
              setbValue(index,bVal);
              break;
            case Surface::Z_POS:
              A = domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr);
              Akm = -domain.getKZM(i,j,k)/domain.getDZM(k);
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
          bVal = foundation.indoorAirTemperature;

          setAmatValue(index,index,A);
          setbValue(index,bVal);
          break;
        case Cell::EXTERIOR_AIR:
          A = 1.0;
          bVal = getOutdoorTemperature();

          setAmatValue(index,index,A);
          setbValue(index,bVal);
          break;
        default:
          {
          if (scheme == Foundation::NS_STEADY_STATE)
          {
            double CXP = domain.cell[i][j][k].cxp;
            double CXM = domain.cell[i][j][k].cxm;
            double CZP = domain.cell[i][j][k].czp;
            double CZM = domain.cell[i][j][k].czm;
            double CYP = domain.cell[i][j][k].cyp;
            double CYM = domain.cell[i][j][k].cym;
            double Q = domain.cell[i][j][k].heatGain;

            if (foundation.coordinateSystem == Foundation::CS_3D ||
              foundation.coordinateSystem == Foundation::CS_3D_SYMMETRY)
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
            else
            {
              double CXPC = 0;
              double CXMC = 0;

              if (i != 0)
              {
                double r = domain.meshX.centers[i];
                CXPC = domain.cell[i][j][k].cxp_c/r;
                CXMC = domain.cell[i][j][k].cxm_c/r;
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
          }
          else
          {
            double theta = timestep/
              (domain.cell[i][j][k].density*domain.cell[i][j][k].specificHeat);

            double f;
            if (scheme == Foundation::NS_IMPLICIT)
              f = 1.0;
            else
              f = 0.5;

            double CXP = domain.cell[i][j][k].cxp*theta;
            double CXM = domain.cell[i][j][k].cxm*theta;
            double CZP = domain.cell[i][j][k].czp*theta;
            double CZM = domain.cell[i][j][k].czm*theta;
            double CYP = domain.cell[i][j][k].cyp*theta;
            double CYM = domain.cell[i][j][k].cym*theta;
            double Q = domain.cell[i][j][k].heatGain*theta;

            if (foundation.coordinateSystem == Foundation::CS_3D ||
              foundation.coordinateSystem == Foundation::CS_3D_SYMMETRY)
            {
              A = (1.0 + f*(CXP + CZP + CYP - CXM - CZM - CYM));
              Aim = f*CXM;
              Aip = f*(-CXP);
              Akm = f*CZM;
              Akp = f*(-CZP);
              Ajm = f*CYM;
              Ajp = f*(-CYP);

              bVal = TOld[i][j][k]*(1.0 + (1-f)*(CXM + CZM + CYM - CXP - CZP - CYP))
                 - TOld[i-1][j][k]*(1-f)*CXM
                 + TOld[i+1][j][k]*(1-f)*CXP
                 - TOld[i][j][k-1]*(1-f)*CZM
                 + TOld[i][j][k+1]*(1-f)*CZP
                 - TOld[i][j-1][k]*(1-f)*CYM
                 + TOld[i][j+1][k]*(1-f)*CYP
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
            else
            {
              double CXPC = 0;
              double CXMC = 0;

              if (i != 0)
              {
                double r = domain.meshX.centers[i];
                CXPC = domain.cell[i][j][k].cxp_c*theta/r;
                CXMC = domain.cell[i][j][k].cxm_c*theta/r;
              }
              A = (1.0 + f*(CXPC + CXP + CZP - CXMC - CXM - CZM));
              Aim = f*(CXMC + CXM);
              Aip = f*(-CXPC - CXP);
              Akm = f*CZM;
              Akp = f*(-CZP);

              bVal = TOld[i][j][k]*(1.0 + (1-f)*(CXMC + CXM + CZM - CXPC - CXP - CZP))
                 - TOld[i-1][j][k]*(1-f)*(CXMC + CXM)
                 + TOld[i+1][j][k]*(1-f)*(CXPC + CXP)
                 - TOld[i][j][k-1]*(1-f)*CZM
                 + TOld[i][j][k+1]*(1-f)*CZP
                 + Q;

              setAmatValue(index,index,A);
              setAmatValue(index,index_ip,Aip);
              setAmatValue(index,index_im,Aim);
              setAmatValue(index,index_kp,Akp);
              setAmatValue(index,index_km,Akm);
              setbValue(index,bVal);
            }
          }
          }
          break;
        }
      }
    }
  }

  solveLinearSystem();

  for (size_t k = 0; k < nZ; ++k)
  {
    for (size_t j = 0; j < nY; ++j)
    {
      for (size_t i = 0; i < nX; ++i)
      {
        int index = i + nX*j + nX*nY*k;
        // Read solution into temperature matrix
        TNew[i][j][k] = getxValue(index);

        // Update old values for next timestep
        TOld[i][j][k] = TNew[i][j][k];
      }
    }
  }

  clearAmat();
}

void Ground::calculateADI(int dim)
{
  for (size_t k = 0; k < nZ; k++)
  {
    for (size_t j = 0; j < nY; j++)
    {
      for (size_t i = 0; i < nX; i++)
      {

        int index;
        if (dim == 1)
          index = i + nX*j + nX*nY*k;
        else if (dim == 2)
          index = j + nY*i + nY*nX*k;
        else if (dim == 3)
          index = k + nZ*i + nZ*nX*j;

        double A, Ap, Am, bVal = 0.0;


        switch (domain.cell[i][j][k].cellType)
        {
        case Cell::BOUNDARY:
          {
          double tilt;
          if (domain.cell[i][j][k].surface.orientation == Surface::Z_POS)
            tilt = 0;
          else if (domain.cell[i][j][k].surface.orientation == Surface::Z_NEG)
            tilt = PI;
          else
            tilt = PI/2.0;

          switch (domain.cell[i][j][k].surface.boundaryConditionType)
          {
          case Surface::ZERO_FLUX:
            {
            switch (domain.cell[i][j][k].surface.orientation)
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
                bVal = TOld[i+1][j][k];
              }

              setAmatValue(index,index,A);
              setAmatValue(index,index+1,Ap);
              setbValue(index,bVal);
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
                bVal = TOld[i-1][j][k];
              }

              setAmatValue(index,index,A);
              setAmatValue(index,index-1,Am);
              setbValue(index,bVal);
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
                bVal = TOld[i][j+1][k];
              }

              setAmatValue(index,index,A);
              setAmatValue(index,index+1,Ap);
              setbValue(index,bVal);
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
                bVal = TOld[i][j-1][k];
              }

              setAmatValue(index,index,A);
              setAmatValue(index,index-1,Am);
              setbValue(index,bVal);
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
                bVal = TOld[i][j][k+1];
              }

              setAmatValue(index,index,A);
              setAmatValue(index,index+1,Ap);
              setbValue(index,bVal);
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
                bVal = TOld[i][j][k-1];
              }

              setAmatValue(index,index,A);
              setAmatValue(index,index-1,Am);
              setbValue(index,bVal);
              break;
            }
            }
            break;
          case Surface::CONSTANT_TEMPERATURE:
            A = 1.0;
            bVal = domain.cell[i][j][k].surface.temperature;

            setAmatValue(index,index,A);
            setbValue(index,bVal);
            break;
          case Surface::INTERIOR_TEMPERATURE:
            A = 1.0;
            bVal = foundation.indoorAirTemperature;

            setAmatValue(index,index,A);
            setbValue(index,bVal);
            break;
          case Surface::EXTERIOR_TEMPERATURE:
            A = 1.0;
            bVal = getOutdoorTemperature();

            setAmatValue(index,index,A);
            setbValue(index,bVal);
            break;
          case Surface::INTERIOR_FLUX:
            {
            double Tair = foundation.indoorAirTemperature;
            double q = 0;

            double hc = getConvectionCoeff(TOld[i][j][k],
                    Tair,0.0,1.52,false,tilt);
            double hr = getSimpleInteriorIRCoeff(domain.cell[i][j][k].surface.emissivity,
                               TOld[i][j][k],Tair);

            switch (domain.cell[i][j][k].surface.orientation)
            {
            case Surface::X_NEG:
              A = domain.getKXP(i,j,k)/domain.getDXP(i) + (hc + hr);
              if (dim == 1)
              {
                Ap = -domain.getKXP(i,j,k)/domain.getDXP(i);
                bVal = (hc + hr)*Tair + q;
              }
              else
              {
                Ap = 0.0;
                bVal = TOld[i+1][j][k]*domain.getKXP(i,j,k)/domain.getDXP(i) + (hc + hr)*Tair + q;
              }

              setAmatValue(index,index,A);
              setAmatValue(index,index+1,Ap);
              setbValue(index,bVal);
              break;
            case Surface::X_POS:
              A = domain.getKXM(i,j,k)/domain.getDXM(i) + (hc + hr);
              if (dim == 1)
              {
                Am = -domain.getKXM(i,j,k)/domain.getDXM(i);
                bVal = (hc + hr)*Tair + q;
              }
              else
              {
                Am = 0.0;
                bVal = TOld[i-1][j][k]*domain.getKXM(i,j,k)/domain.getDXM(i) + (hc + hr)*Tair + q;
              }

              setAmatValue(index,index,A);
              setAmatValue(index,index-1,Am);
              setbValue(index,bVal);
              break;
            case Surface::Y_NEG:
              A = domain.getKYP(i,j,k)/domain.getDYP(j) + (hc + hr);
              if (dim == 2)
              {
                Ap = -domain.getKYP(i,j,k)/domain.getDYP(j);
                bVal = (hc + hr)*Tair + q;
              }
              else
              {
                Ap = 0.0;
                bVal = TOld[i][j+1][k]*domain.getKYP(i,j,k)/domain.getDYP(j) + (hc + hr)*Tair + q;
              }

              setAmatValue(index,index,A);
              setAmatValue(index,index+1,Ap);
              setbValue(index,bVal);
              break;
            case Surface::Y_POS:
              A = domain.getKYM(i,j,k)/domain.getDYM(j) + (hc + hr);
              if (dim == 2)
              {
                Am = -domain.getKYM(i,j,k)/domain.getDYM(j);
                bVal = (hc + hr)*Tair + q;
              }
              else
              {
                Am = 0.0;
                bVal = TOld[i][j-1][k]*domain.getKYM(i,j,k)/domain.getDYM(j) + (hc + hr)*Tair + q;
              }

              setAmatValue(index,index,A);
              setAmatValue(index,index-1,Am);
              setbValue(index,bVal);
              break;
            case Surface::Z_NEG:
              A = domain.getKZP(i,j,k)/domain.getDZP(k) + (hc + hr);
              if (dim == 3)
              {
                Ap = -domain.getKZP(i,j,k)/domain.getDZP(k);
                bVal = (hc + hr)*Tair + q;
              }
              else
              {
                Ap = 0.0;
                bVal = TOld[i][j][k+1]*domain.getKZP(i,j,k)/domain.getDZP(k) + (hc + hr)*Tair + q;
              }

              setAmatValue(index,index,A);
              setAmatValue(index,index+1,Ap);
              setbValue(index,bVal);
              break;
            case Surface::Z_POS:
              A = domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr);
              if (dim == 3)
              {
                Am = -domain.getKZM(i,j,k)/domain.getDZM(k);
                bVal = (hc + hr)*Tair + q;
              }
              else
              {
                Am = 0.0;
                bVal = TOld[i][j][k-1]*domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr)*Tair + q;
              }

              setAmatValue(index,index,A);
              setAmatValue(index,index-1,Am);
              setbValue(index,bVal);
              break;
            }
            }
            break;

          case Surface::EXTERIOR_FLUX:
            {
            double Tair = getOutdoorTemperature();
            double v = getLocalWindSpeed();
            double eSky = weatherData.skyEmissivity.getValue(getSimTime(tNow));
            double F = getEffectiveExteriorViewFactor(eSky,tilt);
            double hc = getConvectionCoeff(TOld[i][j][k],Tair,v,foundation.surfaceRoughness,true,tilt);
            double hr = getExteriorIRCoeff(domain.cell[i][j][k].surface.emissivity,TOld[i][j][k],Tair,eSky,tilt);
            double q = domain.cell[i][j][k].surface.absorptivity*weatherData.globalHorizontalSolar.getValue(getSimTime(tNow));

            switch (domain.cell[i][j][k].surface.orientation)
            {
            case Surface::X_NEG:
              A = domain.getKXP(i,j,k)/domain.getDXP(i) + (hc + hr);
              if (dim == 1)
              {
                Ap = -domain.getKXP(i,j,k)/domain.getDXP(i);
                bVal = (hc + hr*pow(F,0.25))*Tair + q;
              }
              else
              {
                Ap = 0.0;
                bVal = TOld[i+1][j][k]*domain.getKXP(i,j,k)/domain.getDXP(i) + (hc + hr*pow(F,0.25))*Tair + q;
              }

              setAmatValue(index,index,A);
              setAmatValue(index,index+1,Ap);
              setbValue(index,bVal);
              break;
            case Surface::X_POS:
              A = domain.getKXM(i,j,k)/domain.getDXM(i) + (hc + hr);
              if (dim == 1)
              {
                Am = -domain.getKXM(i,j,k)/domain.getDXM(i);
                bVal = (hc + hr*pow(F,0.25))*Tair + q;
              }
              else
              {
                Am = 0.0;
                bVal = TOld[i-1][j][k]*domain.getKXM(i,j,k)/domain.getDXM(i) + (hc + hr*pow(F,0.25))*Tair + q;
              }

              setAmatValue(index,index,A);
              setAmatValue(index,index-1,Am);
              setbValue(index,bVal);
              break;
            case Surface::Y_NEG:
              A = domain.getKYP(i,j,k)/domain.getDYP(j) + (hc + hr);
              if (dim == 2)
              {
                Ap = -domain.getKYP(i,j,k)/domain.getDYP(j);
                bVal = (hc + hr*pow(F,0.25))*Tair + q;
              }
              else
              {
                Ap = 0.0;
                bVal = TOld[i][j+1][k]*domain.getKYP(i,j,k)/domain.getDYP(j) + (hc + hr*pow(F,0.25))*Tair + q;
              }

              setAmatValue(index,index,A);
              setAmatValue(index,index+1,Ap);
              setbValue(index,bVal);
              break;
            case Surface::Y_POS:
              A = domain.getKYM(i,j,k)/domain.getDYM(j) + (hc + hr);
              if (dim == 2)
              {
                Am = -domain.getKYM(i,j,k)/domain.getDYM(j);
                bVal = (hc + hr*pow(F,0.25))*Tair + q;
              }
              else
              {
                Am = 0.0;
                bVal = TOld[i][j-1][k]*domain.getKYM(i,j,k)/domain.getDYM(j) + (hc + hr*pow(F,0.25))*Tair + q;
              }

              setAmatValue(index,index,A);
              setAmatValue(index,index-1,Am);
              setbValue(index,bVal);
              break;
            case Surface::Z_NEG:
              A = domain.getKZP(i,j,k)/domain.getDZP(k) + (hc + hr);
              if (dim == 3)
              {
                Ap = -domain.getKZP(i,j,k)/domain.getDZP(k);
                bVal = (hc + hr*pow(F,0.25))*Tair + q;
              }
              else
              {
                Ap = 0.0;
                bVal = TOld[i][j][k+1]*domain.getKZP(i,j,k)/domain.getDZP(k) + (hc + hr*pow(F,0.25))*Tair + q;
              }

              setAmatValue(index,index,A);
              setAmatValue(index,index+1,Ap);
              setbValue(index,bVal);
              break;
            case Surface::Z_POS:
              A = domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr);
              if (dim == 3)
              {
                Am = -domain.getKZM(i,j,k)/domain.getDZM(k);
                bVal = (hc + hr*pow(F,0.25))*Tair + q;
              }
              else
              {
                Am = 0.0;
                bVal = TOld[i][j][k-1]*domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr*pow(F,0.25))*Tair + q;
              }

              setAmatValue(index,index,A);
              setAmatValue(index,index-1,Am);
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
          bVal = foundation.indoorAirTemperature;

          setAmatValue(index,index,A);
          setbValue(index,bVal);
          break;
        case Cell::EXTERIOR_AIR:
          A = 1.0;
          bVal = getOutdoorTemperature();

          setAmatValue(index,index,A);
          setbValue(index,bVal);
          break;
        default:
          {
          double theta;
          if (foundation.coordinateSystem == Foundation::CS_3D ||
            foundation.coordinateSystem == Foundation::CS_3D_SYMMETRY)
          {
            theta = timestep/
                (3*domain.cell[i][j][k].density*domain.cell[i][j][k].specificHeat);
          }
          else
          {
            theta = timestep/
                (2*domain.cell[i][j][k].density*domain.cell[i][j][k].specificHeat);
          }

          double CXP = domain.cell[i][j][k].cxp*theta;
          double CXM = domain.cell[i][j][k].cxm*theta;
          double CZP = domain.cell[i][j][k].czp*theta;
          double CZM = domain.cell[i][j][k].czm*theta;
          double CYP = domain.cell[i][j][k].cyp*theta;
          double CYM = domain.cell[i][j][k].cym*theta;
          double Q = domain.cell[i][j][k].heatGain*theta;

          double f = foundation.fADI;

          if (foundation.coordinateSystem == Foundation::CS_3D ||
            foundation.coordinateSystem == Foundation::CS_3D_SYMMETRY)
          {
            if (dim == 1) // x
            {
              A = 1.0 + (3 - 2*f)*(CXP - CXM);
              Am = (3 - 2*f)*CXM;
              Ap = (3 - 2*f)*(-CXP);

              bVal = TOld[i][j][k]*(1.0 + f*(CZM + CYM - CZP - CYP))
                   - TOld[i][j][k-1]*f*CZM
                   + TOld[i][j][k+1]*f*CZP
                   - TOld[i][j-1][k]*f*CYM
                   + TOld[i][j+1][k]*f*CYP
                   + Q;
            }
            else if (dim == 2) // y
            {
              A = (1.0 + (3 - 2*f)*(CYP - CYM));
              Am = (3 - 2*f)*CYM;
              Ap = (3 - 2*f)*(-CYP);

              bVal = TOld[i][j][k]*(1.0 + f*(CXM + CZM - CXP - CZP))
                   - TOld[i-1][j][k]*f*CXM
                   + TOld[i+1][j][k]*f*CXP
                   - TOld[i][j][k-1]*f*CZM
                   + TOld[i][j][k+1]*f*CZP
                   + Q;
            }
            else if (dim == 3) // z
            {
              A = (1.0 + (3 - 2*f)*(CZP - CZM));
              Am = (3 - 2*f)*CZM;
              Ap = (3 - 2*f)*(-CZP);

              bVal = TOld[i][j][k]*(1.0 + f*(CXM + CYM - CXP - CYP))
                   - TOld[i-1][j][k]*f*CXM
                   + TOld[i+1][j][k]*f*CXP
                   - TOld[i][j-1][k]*f*CYM
                   + TOld[i][j+1][k]*f*CYP
                   + Q;
            }

          }
          else
          {
            double CXPC = 0;
            double CXMC = 0;
            if (i != 0)
            {
              double r = domain.meshX.centers[i];
              CXPC = domain.cell[i][j][k].cxp_c*theta/r;
              CXMC = domain.cell[i][j][k].cxm_c*theta/r;
            }
            if (dim == 1) // x
            {
              A = 1.0 + (2 - f)*(CXPC + CXP - CXMC - CXM);
              Am = (2 - f)*(CXMC + CXM);
              Ap = (2 - f)*(-CXPC - CXP);

              bVal = TOld[i][j][k]*(1.0 + f*(CZM - CZP))
                   - TOld[i][j][k-1]*f*CZM
                   + TOld[i][j][k+1]*f*CZP
                   + Q;
            }
            else if (dim == 3) // z
            {
              A = 1.0 + (2 - f)*(CZP - CZM);
              Am = (2 - f)*CZM;
              Ap = (2 - f)*(-CZP);

              bVal = TOld[i][j][k]*(1.0 + f*(CXMC + CXM - CXPC - CXP))
                   - TOld[i-1][j][k]*f*(CXMC + CXM)
                   + TOld[i+1][j][k]*f*(CXPC + CXP)
                   + Q;
            }
          }

          setAmatValue(index,index,A);
          setAmatValue(index,index-1,Am);
          setAmatValue(index,index+1,Ap);
          setbValue(index,bVal);

          }
          break;
        }
      }
    }
  }

  solveLinearSystem();

  for (size_t k = 0; k < nZ; ++k)
  {
    for (size_t j = 0; j < nY; ++j)
    {
      for (size_t i = 0; i < nX; ++i)
      {
        int index;
        if (dim == 1)
          index = i + nX*j + nX*nY*k;
        else if (dim == 2)
          index = j + nY*i + nY*nX*k;
        else if (dim == 3)
          index = k + nZ*i + nZ*nX*j;

        // Read solution into temperature matrix
        TNew[i][j][k] = getxValue(index);
        // Update old values for next timestep
        TOld[i][j][k] = TNew[i][j][k];
      }
    }
  }

  clearAmat();
}

void Ground::calculate(double t)
{
  tNow = t;

  // update boundary conditions
  setSolarBoundaryConditions();

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
    calculateADI(1);
    if (foundation.coordinateSystem == Foundation::CS_3D ||
      foundation.coordinateSystem == Foundation::CS_3D_SYMMETRY)
      calculateADI(2);
    calculateADI(3);
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

  plot();

  printStatus(t);
}

void Ground::plot()
{
  for (std::size_t p = 0; p < plots.size(); p++)
  {
    if (plots[p].makeNewFrame(tNow))
    {
      boost::posix_time::ptime time(simulationControl.startDate,boost::posix_time::seconds(tNow));
      std::string timeStamp = to_simple_string(time);
      plots[p].createFrame(TNew, timeStamp);
    }
  }
}

void Ground::printStatus(double t)
{
  boost::posix_time::ptime currentTime = boost::posix_time::second_clock::local_time();
  boost::posix_time::ptime simTime = getSimTime(t);

  if (currentTime - prevStatusUpdate > boost::posix_time::milliseconds(500))
  {
    if (initPeriod)
    {
      std::cout << "  " << simTime << std::endl;
    }
    else
    {
      std::cout << "  " << simTime << " (" << percentComplete << "%)" << std::endl;
    }

    prevStatusUpdate = currentTime;
  }


}

void Ground::setAmatValue(const int i,const int j,const double val)
{
  if (foundation.numericalScheme == Foundation::NS_ADI && TDMA)
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
#if defined(USE_LIS_SOLVER)
    lis_matrix_set_value(LIS_INS_VALUE,i,j,val,Amat);
#else
    Amat(i,j) = val;
#endif
  }
}

void Ground::setbValue(const int i,const double val)
{
  if (foundation.numericalScheme == Foundation::NS_ADI && TDMA)
  {
    b_[i] = val;
  }
  else
  {
#if defined(USE_LIS_SOLVER)
  lis_vector_set_value(LIS_INS_VALUE,i,val,b);
#else
  b(i) = val;
#endif
  }
}

void Ground::solveLinearSystem()
{
  if (foundation.numericalScheme == Foundation::NS_ADI && TDMA)
  {
    solveTDM(a1,a2,a3,b_,x_);
  }
  else
  {
#if defined(USE_LIS_SOLVER)
    lis_matrix_set_type(Amat,LIS_MATRIX_CSR);
    lis_matrix_assemble(Amat);

    lis_solve(Amat,b,x,solver);

    int status;
    lis_solver_get_status(solver, &status);

    if (status != 0) // LIS_MAXITER status
    {
      int iters;
      double residual;

      lis_solver_get_iters(solver, &iters);
      lis_solver_get_residualnorm(solver, &residual);

      std::cout << "Warning: Solution did not converge after ";
      std::cout << iters << " iterations." << "\n";
      std::cout << "  The final residual was: " << residual << "\n";
      std::cout << "  Solver status: " << status << std::endl;

    }
    //lis_output(Amat,b,x,LIS_FMT_MM,"Matrix.mtx");
#else
    umf::umf_solve(Amat,x,b);
#endif
  }
}

void Ground::clearAmat()
{
  if (foundation.numericalScheme == Foundation::NS_ADI && TDMA)
  {
    std::fill(a1.begin(), a1.end(), 0.0);
    std::fill(a2.begin(), a2.end(), 0.0);
    std::fill(a3.begin(), a3.end(), 0.0);
    std::fill(b_.begin(), b_.end(), 0.0);

  }
  else
  {
#if defined(USE_LIS_SOLVER)
  lis_matrix_destroy(Amat);
  lis_matrix_create(LIS_COMM_WORLD,&Amat);
  lis_matrix_set_size(Amat,nX*nY*nZ,nX*nY*nZ);

  lis_vector_destroy(b);
  lis_vector_create(LIS_COMM_WORLD,&b);
  lis_vector_set_size(b,0,nX*nY*nZ);

  lis_solver_destroy(solver);
  lis_solver_create(&solver);
  lis_solver_set_option(&solverOptions[0],solver);

#else
  Amat.clear();
#endif
  }
}

double Ground::getxValue(const int i)
{
  if (foundation.numericalScheme == Foundation::NS_ADI && TDMA)
  {
    return x_[i];
  }
  else
  {
#if defined(USE_LIS_SOLVER)
  double xVal;
  lis_vector_get_value(x,i,&xVal);
  return xVal;
#else
  return x(i);
#endif
  }
}


boost::posix_time::ptime Ground::getSimTime(double t)
{
  boost::posix_time::ptime tp = simulationControl.startTime + boost::posix_time::seconds(t);

  return tp;
}

double Ground::getInitialTemperature(double t, double z)
{
  if (foundation.initializationMethod == Foundation::IM_KUSUDA)
  {
    double minDryBulb = weatherData.dryBulbTemp.getMin();
    double maxDryBulb = weatherData.dryBulbTemp.getMax();

    boost::posix_time::ptime tp = getSimTime(t);
    boost::gregorian::greg_year year = tp.date().year();
    boost::gregorian::date dayBegin(year,boost::gregorian::Jan,1);
    boost::posix_time::ptime tYearStart(dayBegin);
    boost::posix_time::time_duration tSinceYearStart = tp - tYearStart;
    double trel = tSinceYearStart.total_seconds();
    double tshift = 0;  // Assume min temperature occurs at the beginning of the year
    double seconds_in_day = 60.0*60.0*24.0;
    double Tamp = (maxDryBulb - minDryBulb) / 2.0;
    double diff = foundation.soil.conductivity/(foundation.soil.density*foundation.soil.specificHeat);

    return annualAverageDryBulbTemperature
        - Tamp*exp(z*pow(PI/(365*seconds_in_day*diff),0.5))
        *cos(2*PI/(365*seconds_in_day)*(trel - tshift
        - z/2*pow((365*seconds_in_day)/(PI*diff),0.5)));
  }
  else //if (foundation.initializationMethod == Foundation::IM_CONSTANT_TEMPERATURE)
    return foundation.initialTemperature;
}

double Ground::getDeepGroundTemperature()
{
  if (foundation.deepGroundBoundary == Foundation::DGB_AUTO)
    return annualAverageDryBulbTemperature;
  else //if (foundation.deepGroundBoundary == Foundation::DGB_CONSTANT_TEMPERATURE)
    return foundation.deepGroundTemperature;
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

double Ground::getOutdoorTemperature()
{
  if (foundation.outdoorTemperatureMethod == Foundation::OTM_WEATHER_FILE)
    return weatherData.dryBulbTemp.getValue(getSimTime(tNow));
  else // if (foundation.outdoorTemperatureMethod == Foundation::OTM_CONSTANT_TEMPERATURE)
    return foundation.outdoorDryBulbTemperature;
}

double Ground::getLocalWindSpeed()
{
  double vWS = weatherData.windSpeed.getValue(getSimTime(tNow));
  double deltaWS = 270;  // [m]
  double alphaWS = 0.14;
  double deltaLocal = foundation.deltaLocal;  // [m]
  double alphaLocal = foundation.alphaLocal;
  double zWS = 10;  // [m]
  double zLocal = foundation.vegetationHeight;  // [m]

  double vLocal = vWS*pow(deltaWS/zWS,alphaWS)*pow(zLocal/deltaLocal,alphaLocal);

  return vLocal;
}

void Ground::setSolarBoundaryConditions()
{
  for (std::size_t s = 0; s < foundation.surfaces.size() ; s++)
  {
    if (foundation.surfaces[s].name == "Grade"
        || foundation.surfaces[s].name == "Exterior Wall")
    {

      double azi = weatherData.azimuth.getValue(getSimTime(tNow));
      double alt = weatherData.altitude.getValue(getSimTime(tNow));
      double qDN = weatherData.directNormalSolar.getValue(getSimTime(tNow));
      double qGH = weatherData.globalHorizontalSolar.getValue(getSimTime(tNow));
      double qDH = weatherData.diffuseHorizontalSolar.getValue(getSimTime(tNow));
      double pssf;
      double q;

      double incidence;
      double aziYPos = foundation.orientation;
      double aziXPos = PI/2 + foundation.orientation;
      double aziYNeg = PI + foundation.orientation;
      double aziXNeg = 3*PI/2 + foundation.orientation;

      double tilt;
      if (foundation.surfaces[s].orientation == Surface::Z_POS)
      {
        tilt = 0.0;
        incidence = cos(PI/2 - alt);
      }
      else if (foundation.surfaces[s].orientation == Surface::Z_NEG)
      {
        tilt = PI;
        incidence = cos(PI/2 - alt - PI);
      }
      else
      {
        tilt = PI/2.0;

        if (foundation.coordinateSystem == Foundation::CS_2DAXIAL ||
                 foundation.coordinateSystem == Foundation::CS_2DLINEAR)
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
          else if (foundation.surfaces[s].orientation == Surface::X_NEG)
          {
            aziSurf = aziXNeg;
          }

          if (foundation.coordinateSystem == Foundation::CS_3D)
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

#if defined(ENABLE_OPENGL)
      PixelCounter counter(512, 1, false);
#endif

      for (std::size_t index = 0; index < foundation.surfaces[s].indices.size(); index++)
      {
        std::size_t i = boost::get<0>(foundation.surfaces[s].indices[index]);
        std::size_t j = boost::get<1>(foundation.surfaces[s].indices[index]);
        std::size_t k = boost::get<2>(foundation.surfaces[s].indices[index]);

        double alpha = domain.cell[i][j][k].surface.absorptivity;

        if (qGH > 0.0)
        {

#if defined(ENABLE_OPENGL)
          if (isGreaterThan(domain.cell[i][j][k].area, 0.0))
          {

            std::vector<Polygon3> shadedSurface(1);

            double xMin, xMax, yMin, yMax, zMin, zMax;
            xMin = domain.meshX.dividers[i];
            xMax = domain.meshX.dividers[i+1];
            yMin = domain.meshY.dividers[j];
            yMax = domain.meshY.dividers[j+1];
            zMin = domain.meshZ.dividers[k];
            zMax = domain.meshZ.dividers[k+1];


            Polygon3 poly;
            poly.outer().push_back(Point3(xMin,yMin,zMax));
            poly.outer().push_back(Point3(xMin,yMax,zMin));
            poly.outer().push_back(Point3(xMax,yMax,zMax));
            poly.outer().push_back(Point3(xMax,yMin,zMin));
            shadedSurface[0] = poly;

            double areaRatio = counter.getAreaRatio(foundation.orientation,azi,alt,foundation.buildingSurfaces,shadedSurface, 0);

            int pixels = counter.retrievePixelCount(0);
            pssf = areaRatio*pixels;

          }
          else
          {
#else
          pssf = incidence;
#endif
#if defined(ENABLE_OPENGL)
          }
#endif
          q = alpha*(qDN*pssf + qDH*Fsky + qGH*Fg*rho_g);

        }
        else
        {
          q = 0;
        }
        domain.cell[i][j][k].heatGain = q;

      }
    }
  }
}

std::string Ground::printOutputHeaders()
{
  std::string outputHeader = "";

  for (size_t o = 0; o < foundation.outputReport.size(); o++)
  {
    outputHeader += ", " + foundation.outputReport[o].headerText;
  }

  return outputHeader;
}

std::string Ground::printOutputLine()
{
  std::string outputLine = "";

  for (size_t o = 0; o < foundation.outputReport.size(); o++)
  {
    std::string valueString;
    double value;
    if (foundation.outputReport[o].variableID == 0)
    {
      // "Slab Core Average Heat Flux [W/m2]"
      value = getSurfaceAverageHeatFlux("Slab Interior");
      valueString = boost::lexical_cast<std::string>(value);
    }
    else if (foundation.outputReport[o].variableID == 1)
    {
      // "Slab Core Average Temperature [K]"
      value = getSurfaceAverageTemperature("Slab Interior");
      valueString = boost::lexical_cast<std::string>(value);
    }
    else if (foundation.outputReport[o].variableID == 2)
    {
      // "Slab Core Average Effective Temperature [C]"
      value = getSurfaceEffectiveTemperature("Slab Interior",foundation.slab.totalResistance());
      valueString = boost::lexical_cast<std::string>(value);
    }
    else if (foundation.outputReport[o].variableID == 3)
    {
      // "Slab Core Total Heat Transfer Rate [W]"
      value = getSurfaceAverageHeatFlux("Slab Interior")*
          getSurfaceArea("Slab Interior");
      valueString = boost::lexical_cast<std::string>(value);
    }
    else if (foundation.outputReport[o].variableID == 4)
    {
      // "Slab Perimeter Average Heat Flux [W/m2]"
      if (foundation.hasPerimeterSurface)
      {
        value = getSurfaceAverageHeatFlux("Slab Perimeter");
        valueString = boost::lexical_cast<std::string>(value);
      }
      else
        valueString = "null";
    }
    else if (foundation.outputReport[o].variableID == 5)
    {
      // "Slab Perimeter Average Temperature [K]"
      if (foundation.hasPerimeterSurface)
      {
        value = getSurfaceAverageTemperature("Slab Perimeter");
        valueString = boost::lexical_cast<std::string>(value);
      }
      else
        valueString = "null";
    }
      else if (foundation.outputReport[o].variableID == 6)
    {
      // "Slab Perimeter Average Effective Temperature [C]"
      if (foundation.hasPerimeterSurface)
      {
        value = getSurfaceEffectiveTemperature("Slab Perimeter",foundation.slab.totalResistance());
        valueString = boost::lexical_cast<std::string>(value);
      }
      else
        valueString = "null";
    }
      else if (foundation.outputReport[o].variableID == 7)
    {
      // "Slab Perimeter Total Heat Transfer Rate [W]"
      if (foundation.hasPerimeterSurface)
      {
        value = getSurfaceAverageHeatFlux("Slab Perimeter")*
            getSurfaceArea("Slab Perimeter");
        valueString = boost::lexical_cast<std::string>(value);
      }
      else
        valueString = "null";
    }
    else if (foundation.outputReport[o].variableID == 8)
    {
      // "Slab Average Heat Flux [W/m2]"
      double coreArea = getSurfaceArea("Slab Interior");
      double coreTotal = getSurfaceAverageHeatFlux("Slab Interior")*coreArea;
      double perimeterArea = 0.0;
      double perimeterTotal = 0.0;
      if (foundation.hasPerimeterSurface)
      {
        perimeterArea = getSurfaceArea("Slab Perimeter");
        perimeterTotal = getSurfaceAverageHeatFlux("Slab Perimeter")*perimeterArea;
      }

      value = (coreTotal + perimeterTotal)/(coreArea + perimeterArea);
      valueString = boost::lexical_cast<std::string>(value);
    }
    else if (foundation.outputReport[o].variableID == 9)
    {
      // "Slab Average Temperature [K]"
      double coreArea = getSurfaceArea("Slab Interior");
      double coreTotal = getSurfaceAverageTemperature("Slab Interior")*coreArea;
      double perimeterArea = 0.0;
      double perimeterTotal = 0.0;
      if (foundation.hasPerimeterSurface)
      {
        perimeterArea = getSurfaceArea("Slab Perimeter");
        perimeterTotal = getSurfaceAverageTemperature("Slab Perimeter")*perimeterArea;
      }

      value = (coreTotal + perimeterTotal)/(coreArea + perimeterArea);
      valueString = boost::lexical_cast<std::string>(value);
    }
    else if (foundation.outputReport[o].variableID == 10)
    {
      // "Slab Total Heat Transfer Rate [W]"
      double coreTotal = getSurfaceAverageHeatFlux("Slab Interior")*
          getSurfaceArea("Slab Interior");
      double perimeterTotal = 0.0;
      if (foundation.hasPerimeterSurface)
        perimeterTotal = getSurfaceAverageHeatFlux("Slab Perimeter")*
        getSurfaceArea("Slab Perimeter");

      value = coreTotal + perimeterTotal;
      valueString = boost::lexical_cast<std::string>(value);
    }
    else if (foundation.outputReport[o].variableID == 11)
    {
      // "Wall Average Heat Flux [W/m2]"
      if (foundation.foundationDepth > 0.0)
      {
        value = getSurfaceAverageHeatFlux("Interior Wall");
        valueString = boost::lexical_cast<std::string>(value);
      }
      else
        valueString = "null";

    }
    else if (foundation.outputReport[o].variableID == 12)
    {
      // "Wall Average Temperature [K]"
      if (foundation.foundationDepth > 0.0)
      {
        value = getSurfaceAverageTemperature("Interior Wall");
        valueString = boost::lexical_cast<std::string>(value);
      }
      else
        valueString = "null";
    }
    else if (foundation.outputReport[o].variableID == 13)
    {
      // "Wall Average Effective Temperature [C]"
      if (foundation.foundationDepth > 0.0)
      {
        value = getSurfaceEffectiveTemperature("Interior Wall",foundation.wall.totalResistance());
        valueString = boost::lexical_cast<std::string>(value);
      }
      else
        valueString = "null";
    }
    else if (foundation.outputReport[o].variableID == 14)
    {
      // "Wall Total Heat Transfer Rate [W]"
      if (foundation.foundationDepth > 0.0)
      {
        value = getSurfaceAverageHeatFlux("Interior Wall")*
            getSurfaceArea("Interior Wall");
        valueString = boost::lexical_cast<std::string>(value);
      }
      else
        valueString = "null";
    }
    else if (foundation.outputReport[o].variableID == 15)
    {
      // "Foundation Average Heat Flux [W/m2]"
      double coreArea = getSurfaceArea("Slab Interior");
      double coreTotal = getSurfaceAverageHeatFlux("Slab Interior")*coreArea;

      double perimeterArea = 0.0;
      double perimeterTotal = 0.0;
      if (foundation.hasPerimeterSurface)
      {
        perimeterArea = getSurfaceArea("Slab Perimeter");
        perimeterTotal = getSurfaceAverageHeatFlux("Slab Perimeter")*perimeterArea;
      }

      double wallArea = 0.0;
      double wallTotal = 0.0;
      if (foundation.foundationDepth > 0.0)
      {
        wallArea = getSurfaceArea("Interior Wall");
        wallTotal = getSurfaceAverageHeatFlux("Interior Wall")*wallArea;
      }

      value = (coreTotal + perimeterTotal + wallTotal)/(coreArea + perimeterArea + wallArea);
      valueString = boost::lexical_cast<std::string>(value);
    }
    else if (foundation.outputReport[o].variableID == 16)
    {
      // "Foundation Average Temperature [K]"
      double coreArea = getSurfaceArea("Slab Interior");
      double coreTotal = getSurfaceAverageTemperature("Slab Interior")*coreArea;

      double perimeterArea = 0.0;
      double perimeterTotal = 0.0;
      if (foundation.hasPerimeterSurface)
      {
        perimeterArea = getSurfaceArea("Slab Perimeter");
        perimeterTotal = getSurfaceAverageTemperature("Slab Perimeter")*perimeterArea;
      }

      double wallArea = 0.0;
      double wallTotal = 0.0;
      if (foundation.foundationDepth > 0.0)
      {
        wallArea = getSurfaceArea("Interior Wall");
        wallTotal = getSurfaceAverageTemperature("Interior Wall")*wallArea;
      }

      value = (coreTotal + perimeterTotal + wallTotal)/(coreArea + perimeterArea + wallArea);
      valueString = boost::lexical_cast<std::string>(value);
    }
    else if (foundation.outputReport[o].variableID == 17)
    {
      // "Foundation Total Heat Transfer Rate [W]"
      double coreTotal = getSurfaceAverageHeatFlux("Slab Interior")*
          getSurfaceArea("Slab Interior");
      double perimeterTotal = 0.0;
      double wallTotal = 0.0;
      if (foundation.hasPerimeterSurface)
        perimeterTotal = getSurfaceAverageHeatFlux("Slab Perimeter")*
        getSurfaceArea("Slab Perimeter");

      if (foundation.foundationDepth > 0.0)
        wallTotal = getSurfaceAverageHeatFlux("Interior Wall")*
            getSurfaceArea("Interior Wall");

      value = coreTotal + perimeterTotal + wallTotal;
      valueString = boost::lexical_cast<std::string>(value);
    }

    outputLine += ", " + valueString;
  }

  return outputLine;

}


double Ground::getSurfaceAverageHeatFlux(std::string surfaceName)
{
  double totalHeatTransferRate = 0;
  double totalArea = 0;

  // Find surface(s)
  for (size_t s = 0; s < foundation.surfaces.size(); s++)
  {
    if (foundation.surfaces[s].name == surfaceName)
    {
      // Find tilt
      double tilt;
      if (foundation.surfaces[s].orientation == Surface::Z_POS)
        tilt = 0.0;
      else if (foundation.surfaces[s].orientation == Surface::Z_NEG)
        tilt = PI;
      else
        tilt = PI/2.0;

      double Tair = foundation.indoorAirTemperature;

      for (std::size_t index = 0; index < foundation.surfaces[s].indices.size(); index++)
      {
        std::size_t i = boost::get<0>(foundation.surfaces[s].indices[index]);
        std::size_t j = boost::get<1>(foundation.surfaces[s].indices[index]);
        std::size_t k = boost::get<2>(foundation.surfaces[s].indices[index]);

        double h = getConvectionCoeff(TNew[i][j][k],Tair,0.0,1.52,false,tilt)
             + getSimpleInteriorIRCoeff(domain.cell[i][j][k].surface.emissivity,
                 TNew[i][j][k],Tair);

        double A = domain.cell[i][j][k].area;

        totalArea += A;
        totalHeatTransferRate += h*A*(Tair - TNew[i][j][k]);

      }
    }
  }

  return totalHeatTransferRate/totalArea;
}

double Ground::getSurfaceEffectiveTemperature(std::string surfaceName, double constructionRValue)
{
  double totalHeatTransferRate = 0;
  double TA = 0;
  double totalArea = 0;

  double Tair = foundation.indoorAirTemperature;

  // Find surface(s)
  for (size_t s = 0; s < foundation.surfaces.size(); s++)
  {
    if (foundation.surfaces[s].name == surfaceName)
    {
      // Find tilt
      double tilt;
      if (foundation.surfaces[s].orientation == Surface::Z_POS)
        tilt = 0.0;
      else if (foundation.surfaces[s].orientation == Surface::Z_NEG)
        tilt = PI;
      else
        tilt = PI/2.0;


      for (std::size_t index = 0; index < foundation.surfaces[s].indices.size(); index++)
      {
        std::size_t i = boost::get<0>(foundation.surfaces[s].indices[index]);
        std::size_t j = boost::get<1>(foundation.surfaces[s].indices[index]);
        std::size_t k = boost::get<2>(foundation.surfaces[s].indices[index]);

        double h = getConvectionCoeff(TNew[i][j][k],Tair,0.0,1.52,false,tilt)
             + getSimpleInteriorIRCoeff(domain.cell[i][j][k].surface.emissivity,
                 TNew[i][j][k],Tair);

        double A = domain.cell[i][j][k].area;

        totalArea += A;
        totalHeatTransferRate += h*A*(Tair - TNew[i][j][k]);
        TA += TNew[i][j][k]*A;

      }
    }
  }

  double Tavg = TA/totalArea;

  double hAvg = totalHeatTransferRate/(totalArea*(Tair - Tavg));

  return Tair - (totalHeatTransferRate/totalArea)*(constructionRValue+1/hAvg) - 273.15;
}

double Ground::getSurfaceArea(std::string surfaceName)
{
  double totalArea = 0;

  // Find surface(s)
  for (size_t s = 0; s < foundation.surfaces.size(); s++)
  {
    if (foundation.surfaces[s].name == surfaceName)
    {
      Surface surface;
      surface = foundation.surfaces[s];

      totalArea += surface.area;
    }
  }

  return totalArea;
}

double Ground::getSurfaceAverageTemperature(std::string surfaceName)
{
  double TA = 0;
  double totalArea = 0;

  // Find surface(s)
  for (size_t s = 0; s < foundation.surfaces.size(); s++)
  {
    if (foundation.surfaces[s].name == surfaceName)
    {
      for (std::size_t index = 0; index < foundation.surfaces[s].indices.size(); index++)
      {
        std::size_t i = boost::get<0>(foundation.surfaces[s].indices[index]);
        std::size_t j = boost::get<1>(foundation.surfaces[s].indices[index]);
        std::size_t k = boost::get<2>(foundation.surfaces[s].indices[index]);

        double A = domain.cell[i][j][k].area;

        totalArea += A;
        TA += TNew[i][j][k]*A;

      }
    }
  }

  return TA/totalArea;

}

double getArrayValue(boost::multi_array<double, 3> Mat, std::size_t i, std::size_t j, std::size_t k)
{
  return Mat[i][j][k];
}

#endif

