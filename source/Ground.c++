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

#include "Ground.h"

static const double PI = 4.0*atan(1.0);

Ground::Ground(WeatherData &weatherData, Foundation &foundation,
		       SimulationControl &simulationControl) : foundation(foundation),
		       simulationControl(simulationControl), weatherData(weatherData)
{
	// Build Domain Object
	buildDomain();

	// Initial Conditions
	initializeConditions();

	initializePlots();
}

Ground::~Ground()
{
	for (std::size_t p = 0; p < plots.size(); p++)
	{
		plots[p].gr.CloseGIF();
	}
}

void Ground::buildDomain()
{
	timestep = simulationControl.timestep.total_seconds();

	// Create mesh

	foundation.setMeshData();

	// Build matrices for PDE term coefficients
	domain.setDomain(foundation);

	nX = domain.meshX.centers.size();
	nY = domain.meshY.centers.size();
	nZ = domain.meshZ.centers.size();

	//domain.printCellTypes();
}

void Ground::initializeConditions()
{

	tNow = 0.0;
	annualAverageDryBulbTemperature = weatherData.dryBulbTemp.getAverage();

	// Initialize matices
	if (foundation.numericalScheme == Foundation::NS_ADE)
	{
		U.resize(boost::extents[nX][nY][nZ]);
		UOld.resize(boost::extents[nX][nY][nZ]);

		V.resize(boost::extents[nX][nY][nZ]);
		VOld.resize(boost::extents[nX][nY][nZ]);
	}

	if (foundation.numericalScheme == Foundation::NS_CRANK_NICOLSON ||
			 foundation.numericalScheme == Foundation::NS_ADI ||
			 foundation.numericalScheme == Foundation::NS_IMPLICIT ||
			 foundation.numericalScheme == Foundation::NS_STEADY_STATE ||
			 foundation.initializationMethod == Foundation::IM_STEADY_STATE)
	{
		Amat = boost::numeric::ublas::compressed_matrix<double,
		boost::numeric::ublas::column_major, 0,
		boost::numeric::ublas::unbounded_array<int>,
		boost::numeric::ublas::unbounded_array<double> >(nX*nY*nZ,nX*nY*nZ);  // Coefficient Matrix
		b = boost::numeric::ublas::vector<double> (nX*nY*nZ);  // constant and unknown vectors
		x = boost::numeric::ublas::vector<double> (nX*nY*nZ);  // constant and unknown vectors
	}

	TNew.resize(boost::extents[nX][nY][nZ]);
	TOld.resize(boost::extents[nX][nY][nZ]);

	if (foundation.initializationMethod == Foundation::IM_STEADY_STATE)
	{
		calculateMatrix(Foundation::NS_STEADY_STATE);
	}
	else if (foundation.initializationMethod == Foundation::IM_IMPLICIT_ACCELERATION)
	{
		Foundation initialFoundation = foundation;
		initialFoundation.initializationMethod = Foundation::IM_STEADY_STATE;
		initialFoundation.numericalScheme = Foundation::NS_IMPLICIT;
		initialFoundation.outputAnimations.clear();
		SimulationControl initialSimControl;
		initialSimControl.timestep = boost::posix_time::hours(foundation.implicitAccelTimestep);
		initialSimControl.endDate = simulationControl.startDate;
		boost::posix_time::ptime endTime(initialSimControl.endDate);
		endTime = endTime - simulationControl.timestep;
		initialSimControl.startTime = endTime - boost::posix_time::hours(foundation.implicitAccelTimestep*foundation.implicitAccelPeriods);
		initialSimControl.startDate = initialSimControl.startTime.date();

		boost::posix_time::ptime simEnd(endTime);
		boost::posix_time::ptime simStart(initialSimControl.startTime);
		boost::posix_time::time_duration simDuration =  simEnd  - simStart;

		double tstart = 0.0; // [s] Simulation start time
		double tend = simDuration.total_seconds(); // [s] Simulation end time
		double initialTimestep = initialSimControl.timestep.total_seconds();

		Ground initialGround(weatherData,initialFoundation,initialSimControl);

		for (double t = tstart; t <= tend; t = t + initialTimestep)
		{
			initialGround.calculate(t);
		}

		TOld = initialGround.TNew;

	}
	else
	{
		for (size_t k = 0; k < nZ; ++k)
		{
			for (size_t j = 0; j < nY; ++j)
			{
				for (size_t i = 0; i < nX; ++i)
				{
					TOld[i][j][k]= getInitialTemperature(tNow,
							domain.meshZ.centers[k]);
				}
			}
		}
	}
}

void Ground::initializePlots()
{
	for (std::size_t p = 0; p < foundation.outputAnimations.size(); p++)
	{
		if (!foundation.outputAnimations[p].startDateSet)
			foundation.outputAnimations[p].startDate = simulationControl.startDate;

		if (!foundation.outputAnimations[p].endDateSet)
			foundation.outputAnimations[p].endDate = simulationControl.endDate;

		if (foundation.coordinateSystem == Foundation::CS_3D)
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
		plots[p].tEnd = untilEnd.total_seconds();
	}
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
	boost::thread up(boost::bind(&Ground::calculateADEUpwardSweep, this));
	boost::thread down(boost::bind(&Ground::calculateADEDownwardSweep, this));

	up.join();
	down.join();

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
						double q = 0;

						double hc = getConvectionCoeff(TOld[i][j][k],
										Tair,0.0,1.0,false,tilt);
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
						double hc = getConvectionCoeff(TOld[i][j][k],Tair,v,1.0,true,tilt);
						double hr = getExteriorIRCoeff(domain.cell[i][j][k].surface.emissivity,TOld[i][j][k],Tair,eSky,tilt);
						double q = domain.cell[i][j][k].surface.absorptivity*weatherData.globalHorizontalSolar.getValue(getSimTime(tNow));

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

					if (foundation.coordinateSystem == Foundation::CS_3D)
						U[i][j][k] = (UOld[i][j][k]*(1.0 - CXP - CZP - CYP)
								- U[i-1][j][k]*CXM
								+ UOld[i+1][j][k]*CXP
								- U[i][j][k-1]*CZM
								+ UOld[i][j][k+1]*CZP
								- U[i][j-1][k]*CYM
								+ UOld[i][j+1][k]*CYP) /
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
								+ UOld[i][j][k+1]*CZP) /
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
										Tair,0.0,1.0,false,tilt);
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
						double hc = getConvectionCoeff(TOld[i][j][k],Tair,v,1.0,true,tilt);
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

					if (foundation.coordinateSystem == Foundation::CS_3D)
						V[i][j][k] = (VOld[i][j][k]*(1.0 + CXM + CZM + CYM)
								- VOld[i-1][j][k]*CXM
								+ V[i+1][j][k]*CXP
								- VOld[i][j][k-1]*CZM
								+ V[i][j][k+1]*CZP
								- VOld[i][j-1][k]*CYM
								+ V[i][j+1][k]*CYP) /
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
								+ V[i][j][k+1]*CZP) /
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
										Tair,0.0,1.0,false,tilt);
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
						double hc = getConvectionCoeff(TOld[i][j][k],Tair,v,1.0,true,tilt);
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

					if (foundation.coordinateSystem == Foundation::CS_3D)
						TNew[i][j][k] = TOld[i][j][k]*(1.0 + CXM + CZM + CYM - CXP - CZP - CYP)
								- TOld[i-1][j][k]*CXM
								+ TOld[i+1][j][k]*CXP
								- TOld[i][j][k-1]*CZM
								+ TOld[i][j][k+1]*CZP
								- TOld[i][j-1][k]*CYM
								+ TOld[i][j+1][k]*CYP;
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
								+ TOld[i][j][k+1]*CZP;
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
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*k) = 1.0;
							Amat(i + nX*j + nX*nY*k, (i+1) + nX*j + nX*nY*k) = -1.0;

							b(i + nX*j + nX*nY*k) = 0;
							break;
						case Surface::X_POS:
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*k) = 1.0;
							Amat(i + nX*j + nX*nY*k, (i-1) + nX*j + nX*nY*k) = -1.0;

							b(i + nX*j + nX*nY*k) = 0.0;
							break;
						case Surface::Y_NEG:
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*k) = 1.0;
							Amat(i + nX*j + nX*nY*k, i + nX*(j+1) + nX*nY*k) = -1.0;

							b(i + nX*j + nX*nY*k) = 0;
							break;
						case Surface::Y_POS:
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*k) = 1.0;
							Amat(i + nX*j + nX*nY*k, i + nX*(j-1) + nX*nY*k) = -1.0;

							b(i + nX*j + nX*nY*k) = 0.0;
							break;
						case Surface::Z_NEG:
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*k) = 1.0;
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*(k+1)) = -1.0;

							b(i + nX*j + nX*nY*k) = 0;
							break;
						case Surface::Z_POS:
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*k) = 1.0;
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*(k-1)) = -1.0;

							b(i + nX*j + nX*nY*k) = 0;
							break;
						}
						}
						break;

					case Surface::CONSTANT_TEMPERATURE:
						{
						Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*k) = 1.0;

						b(i + nX*j + nX*nY*k) = domain.cell[i][j][k].surface.temperature;
						}
						break;

					case Surface::INTERIOR_TEMPERATURE:
						{
						Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*k) = 1.0;

						b(i + nX*j + nX*nY*k) = foundation.indoorAirTemperature;
						}
						break;

					case Surface::EXTERIOR_TEMPERATURE:
						{
						Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*k) = 1.0;

						b(i + nX*j + nX*nY*k) = getOutdoorTemperature();
						}
						break;

					case Surface::INTERIOR_FLUX:
						{
						double Tair = foundation.indoorAirTemperature;
						double q = 0;

						double hc = getConvectionCoeff(TOld[i][j][k],
										Tair,0.0,1.0,false,tilt);
						double hr = getSimpleInteriorIRCoeff(domain.cell[i][j][k].surface.emissivity,
															 TOld[i][j][k],Tair);

						switch (domain.cell[i][j][k].surface.orientation)
						{
						case Surface::X_NEG:
							{
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*k) = domain.getKXP(i,j,k)/domain.getDXP(i) + (hc + hr);
							Amat(i + nX*j + nX*nY*k, (i+1) + nX*j + nX*nY*k) = -domain.getKXP(i,j,k)/domain.getDXP(i);

							b(i + nX*j + nX*nY*k) = (hc + hr)*Tair + q;
							}
							break;
						case Surface::X_POS:
							{
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*k) = domain.getKXM(i,j,k)/domain.getDXM(i) + (hc + hr);
							Amat(i + nX*j + nX*nY*k, (i-1) + nX*j + nX*nY*k) = -domain.getKXM(i,j,k)/domain.getDXM(i);

							b(i + nX*j + nX*nY*k) = (hc + hr)*Tair + q;
							}
							break;
						case Surface::Y_NEG:
							{
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*k) = domain.getKYP(i,j,k)/domain.getDYP(j) + (hc + hr);
							Amat(i + nX*j + nX*nY*k, i + nX*(j+1) + nX*nY*k) = -domain.getKYP(i,j,k)/domain.getDYP(j);

							b(i + nX*j + nX*nY*k) = (hc + hr)*Tair + q;
							}
							break;
						case Surface::Y_POS:
							{
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*k) = domain.getKYM(i,j,k)/domain.getDYM(j) + (hc + hr);
							Amat(i + nX*j + nX*nY*k, i + nX*(j-1) + nX*nY*k) = -domain.getKYM(i,j,k)/domain.getDYM(j);

							b(i + nX*j + nX*nY*k) = (hc + hr)*Tair + q;
							}
							break;
						case Surface::Z_NEG:
							{
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*k) = domain.getKZP(i,j,k)/domain.getDZP(k) + (hc + hr);
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*(k+1)) = -domain.getKZP(i,j,k)/domain.getDZP(k);

							b(i + nX*j + nX*nY*k) = (hc + hr)*Tair + q;
							}
							break;
						case Surface::Z_POS:
							{
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*k) = domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr);
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*(k-1)) = -domain.getKZM(i,j,k)/domain.getDZM(k);

							b(i + nX*j + nX*nY*k) = (hc + hr)*Tair + q;
							}
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
						double hc = getConvectionCoeff(TOld[i][j][k],Tair,v,1.0,true,tilt);
						double hr = getExteriorIRCoeff(domain.cell[i][j][k].surface.emissivity,TOld[i][j][k],Tair,eSky,tilt);
						double q = domain.cell[i][j][k].surface.absorptivity*weatherData.globalHorizontalSolar.getValue(getSimTime(tNow));

						switch (domain.cell[i][j][k].surface.orientation)
						{
						case Surface::X_NEG:
							{
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*k) = domain.getKXP(i,j,k)/domain.getDXP(i) + (hc + hr);
							Amat(i + nX*j + nX*nY*k, (i+1) + nX*j + nX*nY*k) = -domain.getKXP(i,j,k)/domain.getDXP(i);

							b(i + nX*j + nX*nY*k) = (hc + hr*pow(F,0.25))*Tair + q;
							}
							break;
						case Surface::X_POS:
							{
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*k) = domain.getKXM(i,j,k)/domain.getDXM(i) + (hc + hr);
							Amat(i + nX*j + nX*nY*k, (i-1) + nX*j + nX*nY*k) = -domain.getKXM(i,j,k)/domain.getDXM(i);

							b(i + nX*j + nX*nY*k) = (hc + hr*pow(F,0.25))*Tair + q;
							}
							break;
						case Surface::Y_NEG:
							{
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*k) = domain.getKYP(i,j,k)/domain.getDYP(j) + (hc + hr);
							Amat(i + nX*j + nX*nY*k, i + nX*(j+1) + nX*nY*k) = -domain.getKYP(i,j,k)/domain.getDYP(j);

							b(i + nX*j + nX*nY*k) = (hc + hr*pow(F,0.25))*Tair + q;
							}
							break;
						case Surface::Y_POS:
							{
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*k) = domain.getKYM(i,j,k)/domain.getDYM(j) + (hc + hr);
							Amat(i + nX*j + nX*nY*k, i + nX*(j-1) + nX*nY*k) = -domain.getKYM(i,j,k)/domain.getDYM(j);

							b(i + nX*j + nX*nY*k) = (hc + hr*pow(F,0.25))*Tair + q;
							}
							break;
						case Surface::Z_NEG:
							{
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*k) = domain.getKZP(i,j,k)/domain.getDZP(k) + (hc + hr);
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*(k+1)) = -domain.getKZP(i,j,k)/domain.getDZP(k);

							b(i + nX*j + nX*nY*k) = (hc + hr*pow(F,0.25))*Tair + q;
							}
							break;
						case Surface::Z_POS:
							{
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*k) = domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr);
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*(k-1)) = -domain.getKZM(i,j,k)/domain.getDZM(k);

							b(i + nX*j + nX*nY*k) = (hc + hr*pow(F,0.25))*Tair + q;
							}
							break;
						}
						}
						break;
					}
					}
					break;
				case Cell::INTERIOR_AIR:
					Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*k) = 1.0;

					b(i + nX*j + nX*nY*k) = foundation.indoorAirTemperature;
					break;
				case Cell::EXTERIOR_AIR:
					Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*k) = 1.0;

					b(i + nX*j + nX*nY*k) = getOutdoorTemperature();
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

						if (foundation.coordinateSystem == Foundation::CS_3D)
						{
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*k) = (CXM + CZM + CYM - CXP - CZP - CYP);
							Amat(i + nX*j + nX*nY*k, (i-1) + nX*j + nX*nY*k) = -CXM;
							Amat(i + nX*j + nX*nY*k, (i+1) + nX*j + nX*nY*k) = CXP;
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*(k-1)) = -CZM;
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*(k+1)) = CZP;
							Amat(i + nX*j + nX*nY*k, i + nX*(j-1) + nX*nY*k) = -CYM;
							Amat(i + nX*j + nX*nY*k, i + nX*(j+1) + nX*nY*k) = CYP;

							b(i + nX*j + nX*nY*k) = 0;
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
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*k) = (CXMC + CXM + CZM - CXPC - CXP - CZP);
							Amat(i + nX*j + nX*nY*k, (i-1) + nX*j + nX*nY*k) = (-CXMC - CXM);
							Amat(i + nX*j + nX*nY*k, (i+1) + nX*j + nX*nY*k) = (CXPC + CXP);
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*(k-1)) = -CZM;
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*(k+1)) = CZP;

							b(i + nX*j + nX*nY*k) = 0;
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

						if (foundation.coordinateSystem == Foundation::CS_3D)
						{
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*k) = (1.0 + f*(CXP + CZP + CYP - CXM - CZM - CYM));
							Amat(i + nX*j + nX*nY*k, (i-1) + nX*j + nX*nY*k) = f*CXM;
							Amat(i + nX*j + nX*nY*k, (i+1) + nX*j + nX*nY*k) = f*(-CXP);
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*(k-1)) = f*CZM;
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*(k+1)) = f*(-CZP);
							Amat(i + nX*j + nX*nY*k, i + nX*(j-1) + nX*nY*k) = f*CYM;
							Amat(i + nX*j + nX*nY*k, i + nX*(j+1) + nX*nY*k) = f*(-CYP);

							b(i + nX*j + nX*nY*k)
									= TOld[i][j][k]*(1.0 + (1-f)*(CXM + CZM + CYM - CXP - CZP - CYP))
									- TOld[i-1][j][k]*(1-f)*CXM
									+ TOld[i+1][j][k]*(1-f)*CXP
									- TOld[i][j][k-1]*(1-f)*CZM
									+ TOld[i][j][k+1]*(1-f)*CZP
									- TOld[i][j-1][k]*(1-f)*CYM
									+ TOld[i][j+1][k]*(1-f)*CYP;
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
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*k) = (1.0 + f*(CXPC + CXP + CZP - CXMC - CXM - CZM));
							Amat(i + nX*j + nX*nY*k, (i-1) + nX*j + nX*nY*k) = f*(CXMC + CXM);
							Amat(i + nX*j + nX*nY*k, (i+1) + nX*j + nX*nY*k) = f*(-CXPC - CXP);
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*(k-1)) = f*CZM;
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*(k+1)) = f*(-CZP);

							b(i + nX*j + nX*nY*k)
									= TOld[i][j][k]*(1.0 + (1-f)*(CXMC + CXM + CZM - CXPC - CXP - CZP))
									- TOld[i-1][j][k]*(1-f)*(CXMC + CXM)
									+ TOld[i+1][j][k]*(1-f)*(CXPC + CXP)
									- TOld[i][j][k-1]*(1-f)*CZM
									+ TOld[i][j][k+1]*(1-f)*CZP;
						}
					}
					}
					break;
				}
			}
		}
	}

	umf::umf_solve(Amat,x,b);

	for (size_t k = 0; k < nZ; ++k)
	{
		for (size_t j = 0; j < nY; ++j)
		{
			for (size_t i = 0; i < nX; ++i)
			{
				// Read solution into temperature matrix
				TNew[i][j][k] = x(i + nX*j + nX*nY*k);

				// Update old values for next timestep
				TOld[i][j][k] = TNew[i][j][k];
			}
		}
	}
	Amat.clear();
	b.clear();
	x.clear();

}

void Ground::calculateADI(int dim)
{
	for (size_t k = 0; k < nZ; k++)
	{
		for (size_t j = 0; j < nY; j++)
		{
			for (size_t i = 0; i < nX; i++)
			{

				std::size_t index;
				if (dim == 1)
					index = i + nX*j + nX*nY*k;
				else if (dim == 2)
					index = j + nY*i + nY*nX*k;
				else if (dim == 3)
					index = k + nZ*i + nZ*nX*j;

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
							Amat(index, index) = 1.0;
							if (dim == 1)
							{
								Amat(index, index+1) = -1.0;
								b(index) = 0;
							}
							else
							{
								b(index) = TOld[i+1][j][k];
							}
							break;
						case Surface::X_POS:
							Amat(index, index) = 1.0;
							if (dim == 1)
							{
								Amat(index, index-1) = -1.0;
								b(index) = 0;
							}
							else
							{
								b(index) = TOld[i-1][j][k];
							}
							break;
						case Surface::Y_NEG:
							Amat(index, index) = 1.0;
							if (dim == 2)
							{
								Amat(index, index+1) = -1.0;
								b(index) = 0;
							}
							else
							{
								b(index) = TOld[i][j+1][k];
							}
							break;
						case Surface::Y_POS:
							Amat(index, index) = 1.0;
							if (dim == 2)
							{
								Amat(index, index-1) = -1.0;
								b(index) = 0;
							}
							else
							{
								b(index) = TOld[i][j-1][k];
							}
							break;
						case Surface::Z_NEG:
							Amat(index, index) = 1.0;
							if (dim == 3)
							{
								Amat(index, index+1) = -1.0;
								b(index) = 0;
							}
							else
							{
								b(index) = TOld[i][j][k+1];
							}
							break;
						case Surface::Z_POS:
							Amat(index, index) = 1.0;
							if (dim == 3)
							{
								Amat(index, index-1) = -1.0;
								b(index) = 0;
							}
							else
							{
								b(index) = TOld[i][j][k-1];
							}
							break;
						}
						}
						break;

					case Surface::CONSTANT_TEMPERATURE:
						{
						Amat(index, index) = 1.0;

						b(index) = domain.cell[i][j][k].surface.temperature;
						}
						break;

					case Surface::INTERIOR_TEMPERATURE:
						{
						Amat(index, index) = 1.0;

						b(index) = foundation.indoorAirTemperature;
						}
						break;

					case Surface::EXTERIOR_TEMPERATURE:
						{
						Amat(index, index) = 1.0;

						b(index) = getOutdoorTemperature();
						}
						break;

					case Surface::INTERIOR_FLUX:
						{
						double Tair = foundation.indoorAirTemperature;
						double q = 0;

						double hc = getConvectionCoeff(TOld[i][j][k],
										Tair,0.0,1.0,false,tilt);
						double hr = getSimpleInteriorIRCoeff(domain.cell[i][j][k].surface.emissivity,
															 TOld[i][j][k],Tair);

						switch (domain.cell[i][j][k].surface.orientation)
						{
						case Surface::X_NEG:
							{
							Amat(index, index) = domain.getKXP(i,j,k)/domain.getDXP(i) + (hc + hr);
							if (dim == 1)
							{
								Amat(index, index+1) = -domain.getKXP(i,j,k)/domain.getDXP(i);
								b(index) = (hc + hr)*Tair + q;
							}
							else
							{
								b(index) = TOld[i+1][j][k]*domain.getKXP(i,j,k)/domain.getDXP(i) + (hc + hr)*Tair + q;
							}
							}
							break;
						case Surface::X_POS:
							{
							Amat(index, index) = domain.getKXM(i,j,k)/domain.getDXM(i) + (hc + hr);
							if (dim == 1)
							{
								Amat(index, index-1) = -domain.getKXM(i,j,k)/domain.getDXM(i);
								b(index) = (hc + hr)*Tair + q;
							}
							else
							{
								b(index) = TOld[i-1][j][k]*domain.getKXM(i,j,k)/domain.getDXM(i) + (hc + hr)*Tair + q;
							}
							}
							break;
						case Surface::Y_NEG:
							{
							Amat(index, index) = domain.getKYP(i,j,k)/domain.getDYP(j) + (hc + hr);
							if (dim == 2)
							{
								Amat(index, index+1) = -domain.getKYP(i,j,k)/domain.getDYP(j);
								b(index) = (hc + hr)*Tair + q;
							}
							else
							{
								b(index) = TOld[i][j+1][k]*domain.getKYP(i,j,k)/domain.getDYP(j) + (hc + hr)*Tair + q;
							}
							}
							break;
						case Surface::Y_POS:
							{
							Amat(index, index) = domain.getKYM(i,j,k)/domain.getDYM(j) + (hc + hr);
							if (dim == 2)
							{
								Amat(index, index-1) = -domain.getKYM(i,j,k)/domain.getDYM(j);
								b(index) = (hc + hr)*Tair + q;
							}
							else
							{
								b(index) = TOld[i][j-1][k]*domain.getKYM(i,j,k)/domain.getDYM(j) + (hc + hr)*Tair + q;
							}
							}
							break;
						case Surface::Z_NEG:
							{
							Amat(index, index) = domain.getKZP(i,j,k)/domain.getDZP(k) + (hc + hr);
							if (dim == 3)
							{
								Amat(index, index+1) = -domain.getKZP(i,j,k)/domain.getDZP(k);
								b(index) = (hc + hr)*Tair + q;
							}
							else
							{
								b(index) = TOld[i][j][k+1]*domain.getKZP(i,j,k)/domain.getDZP(k) + (hc + hr)*Tair + q;
							}
							}
							break;
						case Surface::Z_POS:
							{
							Amat(index, index) = domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr);
							if (dim == 3)
							{
								Amat(index, index-1) = -domain.getKZM(i,j,k)/domain.getDZM(k);
								b(index) = (hc + hr)*Tair + q;
							}
							else
							{
								b(index) = TOld[i][j][k-1]*domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr)*Tair + q;
							}
							}
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
						double hc = getConvectionCoeff(TOld[i][j][k],Tair,v,1.0,true,tilt);
						double hr = getExteriorIRCoeff(domain.cell[i][j][k].surface.emissivity,TOld[i][j][k],Tair,eSky,tilt);
						double q = domain.cell[i][j][k].surface.absorptivity*weatherData.globalHorizontalSolar.getValue(getSimTime(tNow));

						switch (domain.cell[i][j][k].surface.orientation)
						{
						case Surface::X_NEG:
							{
							Amat(index, index) = domain.getKXP(i,j,k)/domain.getDXP(i) + (hc + hr);
							if (dim == 1)
							{
								Amat(index, index+1) = -domain.getKXP(i,j,k)/domain.getDXP(i);
								b(index) = (hc + hr*pow(F,0.25))*Tair + q;
							}
							else
							{
								b(index) = TOld[i+1][j][k]*domain.getKXP(i,j,k)/domain.getDXP(i) + (hc + hr*pow(F,0.25))*Tair + q;
							}
							}
							break;
						case Surface::X_POS:
							{
							Amat(index, index) = domain.getKXM(i,j,k)/domain.getDXM(i) + (hc + hr);
							if (dim == 1)
							{
								Amat(index, index-1) = -domain.getKXM(i,j,k)/domain.getDXM(i);
								b(index) = (hc + hr*pow(F,0.25))*Tair + q;
							}
							else
							{
								b(index) = TOld[i-1][j][k]*domain.getKXM(i,j,k)/domain.getDXM(i) + (hc + hr*pow(F,0.25))*Tair + q;
							}
							}
							break;
						case Surface::Y_NEG:
							{
							Amat(index, index) = domain.getKYP(i,j,k)/domain.getDYP(j) + (hc + hr);
							if (dim == 2)
							{
								Amat(index, index+1) = -domain.getKYP(i,j,k)/domain.getDYP(j);
								b(index) = (hc + hr*pow(F,0.25))*Tair + q;
							}
							else
							{
								b(index) = TOld[i][j+1][k]*domain.getKYP(i,j,k)/domain.getDYP(j) + (hc + hr*pow(F,0.25))*Tair + q;
							}
							}
							break;
						case Surface::Y_POS:
							{
							Amat(index, index) = domain.getKYM(i,j,k)/domain.getDYM(j) + (hc + hr);
							if (dim == 2)
							{
								Amat(index, index-1) = -domain.getKYM(i,j,k)/domain.getDYM(j);
								b(index) = (hc + hr*pow(F,0.25))*Tair + q;
							}
							else
							{
								b(index) = TOld[i][j-1][k]*domain.getKYM(i,j,k)/domain.getDYM(j) + (hc + hr*pow(F,0.25))*Tair + q;
							}
							}
							break;
						case Surface::Z_NEG:
							{
							Amat(index, index) = domain.getKZP(i,j,k)/domain.getDZP(k) + (hc + hr);
							if (dim == 3)
							{
								Amat(index, index+1) = -domain.getKZP(i,j,k)/domain.getDZP(k);
								b(index) = (hc + hr*pow(F,0.25))*Tair + q;
							}
							else
							{
								b(index) = TOld[i][j][k+1]*domain.getKZP(i,j,k)/domain.getDZP(k) + (hc + hr*pow(F,0.25))*Tair + q;
							}
							}
							break;
						case Surface::Z_POS:
							{
							Amat(index, index) = domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr);
							if (dim == 3)
							{
								Amat(index, index-1) = -domain.getKZM(i,j,k)/domain.getDZM(k);
								b(index) = (hc + hr*pow(F,0.25))*Tair + q;
							}
							else
							{
								b(index) = TOld[i][j][k-1]*domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr*pow(F,0.25))*Tair + q;
							}
							}
							break;
						}
						}
						break;
					}
					}
					break;
				case Cell::INTERIOR_AIR:
					Amat(index, index) = 1.0;

					b(index) = foundation.indoorAirTemperature;
					break;
				case Cell::EXTERIOR_AIR:
					Amat(index, index) = 1.0;

					b(index) = getOutdoorTemperature();
					break;
				default:
					{
					double theta;
					if (foundation.coordinateSystem == Foundation::CS_3D)
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

					double f = foundation.fADI;

					if (foundation.coordinateSystem == Foundation::CS_3D)
					{
						if (dim == 1) // x
						{
							Amat(index, index) = 1.0 + (3 - 2*f)*(CXP - CXM);
							Amat(index, index-1) = (3 - 2*f)*CXM;
							Amat(index, index+1) = (3 - 2*f)*(-CXP);

							b(index)
									= TOld[i][j][k]*(1.0 + f*(CZM + CYM - CZP - CYP))
									- TOld[i][j][k-1]*f*CZM
									+ TOld[i][j][k+1]*f*CZP
									- TOld[i][j-1][k]*f*CYM
									+ TOld[i][j+1][k]*f*CYP;
						}
						else if (dim == 2) // y
						{
							Amat(index, index) = (1.0 + (3 - 2*f)*(CYP - CYM));
							Amat(index, index-1) = (3 - 2*f)*CYM;
							Amat(index, index+1) = (3 - 2*f)*(-CYP);

							b(index)
									= TOld[i][j][k]*(1.0 + f*(CXM + CZM - CXP - CZP))
									- TOld[i-1][j][k]*f*CXM
									+ TOld[i+1][j][k]*f*CXP
									- TOld[i][j][k-1]*f*CZM
									+ TOld[i][j][k+1]*f*CZP;
						}
						else if (dim == 3) // z
						{
							Amat(index, index) = (1.0 + (3 - 2*f)*(CZP - CZM));
							Amat(index, index-1) = (3 - 2*f)*CZM;
							Amat(index, index+1) = (3 - 2*f)*(-CZP);

							b(index)
									= TOld[i][j][k]*(1.0 + f*(CXM + CYM - CXP - CYP))
									- TOld[i-1][j][k]*f*CXM
									+ TOld[i+1][j][k]*f*CXP
									- TOld[i][j-1][k]*f*CYM
									+ TOld[i][j+1][k]*f*CYP;
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
							Amat(index, index) = 1.0 + (2 - f)*(CXPC + CXP - CXMC - CXM);
							Amat(index, index-1) = (2 - f)*(CXMC + CXM);
							Amat(index, index+1) = (2 - f)*(-CXPC - CXP);

							b(index)
									= TOld[i][j][k]*(1.0 + f*(CZM - CZP))
									- TOld[i][j][k-1]*f*CZM
									+ TOld[i][j][k+1]*f*CZP;
						}
						else if (dim == 3) // z
						{
							Amat(index, index) = 1.0 + (2 - f)*(CZP - CZM);
							Amat(index, index-1) = (2 - f)*CZM;
							Amat(index, index+1) = (2 - f)*(-CZP);

							b(index)
									= TOld[i][j][k]*(1.0 + f*(CXMC + CXM - CXPC - CXP))
									- TOld[i-1][j][k]*f*(CXMC + CXM)
									+ TOld[i+1][j][k]*f*(CXPC + CXP);
						}
					}
					}
					break;
				}
			}
		}
	}

	umf::umf_solve(Amat,x,b);

	for (size_t k = 0; k < nZ; ++k)
	{
		for (size_t j = 0; j < nY; ++j)
		{
			for (size_t i = 0; i < nX; ++i)
			{
				std::size_t index;
				if (dim == 1)
					index = i + nX*j + nX*nY*k;
				else if (dim == 2)
					index = j + nY*i + nY*nX*k;
				else if (dim == 3)
					index = k + nZ*i + nZ*nX*j;

				// Read solution into temperature matrix
				TNew[i][j][k] = x(index);
				// Update old values for next timestep
				TOld[i][j][k] = TNew[i][j][k];
			}
		}
	}
	Amat.clear();
	b.clear();
	x.clear();

}

void Ground::calculate(double t)
{
	tNow = t;

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
		if (foundation.coordinateSystem == Foundation::CS_3D)
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
	// Calculate Heat Fluxes
	QSlabTotal = getSurfaceAverageHeatFlux("Slab Interior");

	plot();
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

double Ground::getSurfaceAverageHeatFlux(std::string surfaceName)
{
	// Find surface
	Surface surface;
	for (size_t s = 0; s < foundation.surfaces.size(); s++)
	{
		if (foundation.surfaces[s].name == surfaceName)
		{
			surface = foundation.surfaces[s];
		}
	}

	size_t iMin, iMax, jMin, jMax, kMin, kMax;
	double tilt;
	double totalArea;

	// Find bounding indices
	if (surface.orientation == Surface::X_POS ||
		surface.orientation == Surface::X_NEG)
	{
		iMin = domain.meshX.getNearestIndex(surface.xMin);
		iMax = domain.meshX.getNearestIndex(surface.xMax);
		jMin = domain.meshY.getNextIndex(surface.yMin);
		jMax = domain.meshY.getPreviousIndex(surface.yMax);
		kMin = domain.meshZ.getNextIndex(surface.zMin);
		kMax = domain.meshZ.getPreviousIndex(surface.zMax);
	}
	else if (surface.orientation == Surface::Y_POS ||
		surface.orientation == Surface::Y_NEG)
	{
		iMin = domain.meshX.getNextIndex(surface.xMin);
		iMax = domain.meshX.getPreviousIndex(surface.xMax);
		jMin = domain.meshY.getNearestIndex(surface.yMin);
		jMax = domain.meshY.getNearestIndex(surface.yMax);
		kMin = domain.meshZ.getNextIndex(surface.zMin);
		kMax = domain.meshZ.getPreviousIndex(surface.zMax);
	}
	else // if (surface.orientation == Surface::Z_POS ||
		 // surface.orientation == Surface::Z_NEG)
	{
		iMin = domain.meshX.getNextIndex(surface.xMin);
		iMax = domain.meshX.getPreviousIndex(surface.xMax);
		jMin = domain.meshY.getNextIndex(surface.yMin);
		jMax = domain.meshY.getPreviousIndex(surface.yMax);
		kMin = domain.meshZ.getNearestIndex(surface.zMin);
		kMax = domain.meshZ.getNearestIndex(surface.zMax);
	}

	// Find total area
	if (foundation.coordinateSystem == Foundation::CS_2DAXIAL)
	{
		if (surface.orientation == Surface::X_POS ||
			surface.orientation == Surface::X_NEG)
		{
			totalArea = 2.0*PI*surface.xMax*(surface.zMax - surface.zMin);
		}
		else // if (surface.orientation == Surface::Z_POS ||
			 // surface.orientation == Surface::Z_NEG)
		{
			totalArea = PI*pow(surface.xMax,2.0) - PI*pow(surface.xMin,2.0);
		}
	}
	else if (foundation.coordinateSystem == Foundation::CS_2DLINEAR)
	{
		if (surface.orientation == Surface::X_POS ||
			surface.orientation == Surface::X_NEG)
		{
			totalArea = surface.zMax - surface.zMin;
		}
		else // if (surface.orientation == Surface::Z_POS ||
			 // surface.orientation == Surface::Z_NEG)
		{
			totalArea = surface.xMax - surface.xMin;
		}
	}
	else // if (foundation.coordinateSystem == Foundation::CS_3D)
	{
		if (surface.orientation == Surface::X_POS ||
			surface.orientation == Surface::X_NEG)
		{
			totalArea = (surface.yMax - surface.yMin)*(surface.zMax - surface.zMin);
		}
		else if (surface.orientation == Surface::Y_POS ||
			surface.orientation == Surface::Y_NEG)
		{
			totalArea = (surface.xMax - surface.xMin)*(surface.zMax - surface.zMin);
		}
		else // if (surface.orientation == Surface::Z_POS ||
			 // surface.orientation == Surface::Z_NEG)
		{
			totalArea = (surface.xMax - surface.xMin)*(surface.yMax - surface.yMin);
		}

	}

	// Find tilt
	if (surface.orientation == Surface::Z_POS)
		tilt = 0.0;
	else if (surface.orientation == Surface::Z_NEG)
		tilt = PI;
	else
		tilt = PI/2.0;

	double Tair = foundation.indoorAirTemperature;

	std::vector<double> heatFlux;

	// Loop over cells and calculate heat loss from surface
	for (size_t k = kMin; k <= kMax; ++k)
	{
		for (size_t j = jMin; j <= jMax; ++j)
		{
			for (size_t i = iMin; i <= iMax; ++i)
			{
				double h = getConvectionCoeff(TNew[i][j][k],Tair,0.0,1.0,false,tilt)
						 + getSimpleInteriorIRCoeff(domain.cell[i][j][k].surface.emissivity,
								 TNew[i][j][k],Tair);

				// Calculate Area
				double A;
				if (foundation.coordinateSystem == Foundation::CS_2DAXIAL)
				{
					if (surface.orientation == Surface::X_POS ||
						surface.orientation == Surface::X_NEG)
					{
						A = 2.0*PI*domain.meshX.centers[i]*domain.meshZ.deltas[k];
					}
					else // if (surface.orientation == Surface::Z_POS ||
						 // surface.orientation == Surface::Z_NEG)
					{
						A = 2.0*PI*domain.meshX.deltas[i]*domain.meshX.centers[i];
					}
				}
				else if (foundation.coordinateSystem == Foundation::CS_2DLINEAR)
				{
					if (surface.orientation == Surface::X_POS ||
						surface.orientation == Surface::X_NEG)
					{
						A = domain.meshZ.deltas[k];
					}
					else // if (surface.orientation == Surface::Z_POS ||
						 // surface.orientation == Surface::Z_NEG)
					{
						A = domain.meshX.deltas[i];
					}
				}
				else  // if (foundation.coordinateSystem == Foundation::CS_3D)
				{
					if (surface.orientation == Surface::X_POS ||
						surface.orientation == Surface::X_NEG)
					{
						A = domain.meshY.deltas[j]*domain.meshZ.deltas[k];
					}
					else if (surface.orientation == Surface::Y_POS ||
						surface.orientation == Surface::Y_NEG)
					{
						A = domain.meshX.deltas[i]*domain.meshZ.deltas[k];
					}
					else // if (surface.orientation == Surface::Z_POS ||
						 // surface.orientation == Surface::Z_NEG)
					{
						A = domain.meshX.deltas[i]*domain.meshY.deltas[j];
					}
				}
				heatFlux.push_back(h*A*(Tair - TNew[i][j][k]));
			}
		}
	}

	double totalFlux = std::accumulate((heatFlux).begin(),(heatFlux).end(), 0.0);

	return totalFlux/totalArea;
}

double getArrayValue(boost::multi_array<double, 3> Mat, std::size_t i, std::size_t j, std::size_t k)
{
	return Mat[i][j][k];
}

#endif

