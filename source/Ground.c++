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

	makePlot = true;
	if (foundation.outputAnimation.name == "") makePlot = false;

	if (makePlot) initializePlot();
}

Ground::~Ground()
{
	if (makePlot) gr.CloseGIF();
}

void Ground::buildDomain()
{
	timestep = simulationControl.timestep.total_seconds();

	// Create mesh

	foundation.setMeshData();

	// Build matrices for PDE term coefficients
	domain.setDomain(foundation, timestep);

	nX = domain.meshX.centers.size();
	nY = domain.meshY.centers.size();
	nZ = domain.meshZ.centers.size();

	//cout << "Number of Cells: " << nX << " x " << nZ << " = " << nX*nZ << endl;
	domain.printCellTypes();
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
	else if (foundation.numericalScheme == Foundation::NS_CRANK_NICOLSON ||
			 foundation.numericalScheme == Foundation::NS_IMPLICIT ||
			 foundation.numericalScheme == Foundation::NS_STEADY_STATE)
	{
		Amat = boost::numeric::ublas::compressed_matrix<double,
		boost::numeric::ublas::column_major, 0,
		boost::numeric::ublas::unbounded_array<int>,
		boost::numeric::ublas::unbounded_array<double> >(nX*nZ,nX*nZ);  // Coefficient Matrix
		b = boost::numeric::ublas::vector<double> (nX*nZ);  // constant and unknown vectors
		x = boost::numeric::ublas::vector<double> (nX*nZ);  // constant and unknown vectors
	}

	TNew.resize(boost::extents[nX][nY][nZ]);
	TOld.resize(boost::extents[nX][nY][nZ]);

	QSlab = std::vector<double>(domain.slabImax - domain.slabImin + 1);

	if (foundation.initializationMethod == Foundation::IM_STEADY_STATE)
		calculateSteadyState();

	else if (foundation.initializationMethod == Foundation::IM_IMPLICIT_ACCELERATION)
	{
		Foundation initialFoundation = foundation;
		initialFoundation.initializationMethod = Foundation::IM_STEADY_STATE;
		initialFoundation.numericalScheme = Foundation::NS_IMPLICIT;
		initialFoundation.outputAnimation.name = "";
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
					switch (domain.cellType[i][j][k])
					{
					case Domain::INTERIOR_AIR:
						TOld[i][j][k] = foundation.indoorAirTemperature;
						break;
					case Domain::EXTERIOR_AIR:
						TOld[i][j][k] = getOutdoorTemperature();
						break;
					case Domain::DEEP_GROUND:
						TOld[i][j][k] = getDeepGroundTemperature();
						break;
					default:
						TOld[i][j][k]= getInitialTemperature(tNow,
														domain.meshZ.centers[k]);
						break;
					}
				}
			}
		}
	}
}

void Ground::initializePlot()
{

	startString = boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time());

	size_t contourLevels = 13;

	mglData TRef(nX, nZ),
			rRef(nX),
			zRef(nZ),
			rGridRef(nX + 1),
			zGridRef(nZ + 1),
			TGridRef(nX + 1, nZ + 1),
			cRef(contourLevels);


	TDat = TRef;
	rDat = rRef;
	zDat = zRef;
	rGrid = rGridRef;
	zGrid = zGridRef;
	TGrid = TGridRef;
	cDat = cRef;

	nextPlotTime = 0.0;
	plotFreq = foundation.outputAnimation.frequency.total_seconds();

	double aspect = 1.0;
	int height = foundation.outputAnimation.size;
	int width = height*aspect;
	gr.SetSize(width,height);

	gr.StartGIF((foundation.outputAnimation.name + ".gif").c_str(),100);

	rGrid.a[0] = domain.meshX.dividers[0];

	for(size_t i = 0; i < nX; i++)
	{
		rDat.a[i] = domain.meshX.centers[i];
		rGrid.a[i + 1] = domain.meshX.dividers[i + 1];
	}

	zGrid.a[0] = domain.meshZ.dividers[0];

	for(size_t k = 0; k < nZ; k++)
	{
		zDat.a[k] = domain.meshZ.centers[k];
		zGrid.a[k + 1] = domain.meshZ.dividers[k + 1];
	}

	for(size_t k = 0; k <= nZ; k++)
	{
		for(size_t i = 0; i <= nX; i++)
		{
			TGrid.a[i+nX*k] = 200.0;
		}
	}

	for (size_t n = 0; n < contourLevels; n++)
	{
		double mid = 50.0;
		double range = 120.0;
		double step = range / (contourLevels - 1);
		cDat.a[n] = mid - range/2 + double(n)*step;
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
				double r = domain.meshX.centers[i];
				double A = domain.cxp_c[i][j][k]*domain.theta[i][j][k]/r;
				double B = domain.cxm_c[i][j][k]*domain.theta[i][j][k]/r;
				double C = domain.cxp[i][j][k]*domain.theta[i][j][k];
				double D = domain.cxm[i][j][k]*domain.theta[i][j][k];
				double E = domain.czp[i][j][k]*domain.theta[i][j][k];
				double F = domain.czm[i][j][k]*domain.theta[i][j][k];

				switch (domain.cellType[i][j][k])
				{
				case Domain::SYMMETRY:
					// Apply zero-flux BC (U[i-1][j][k] = UOld[i+1][j][k])
					// Remove A & B at axis to avoid divide by zero
					U[i][j][k] = (UOld[i][j][k]*(1.0 - C - E) - UOld[i+1][j][k]*D
							  + UOld[i+1][j][k]*C - U[i][j][k-1]*F + UOld[i][j][k+1]*E)
							  /(1.0 - D - F);
					break;
				case Domain::FAR_FIELD:
					// Apply zero-flux BC (UOld[i+1][j][k] = U[i-1][j][k])
					U[i][j][k] = (UOld[i][j][k]*(1.0 - A - C - E) - U[i-1][j][k]*(B + D)
							  + U[i-1][j][k]*(A + C) - U[i][j][k-1]*F + UOld[i][j][k+1]*E)
							  /(1.0 - B - D - F);
					break;
				case Domain::WALL_TOP:
					// Apply zero-flux BC (UOld[i][j][k+1] = U[i][j][k-1]
					U[i][j][k] = (UOld[i][j][k]*(1.0 - A - C - E) - U[i-1][j][k]*(B + D)
							  + UOld[i+1][j][k]*(A + C) - U[i][j][k-1]*F + U[i][j][k-1]*E)
							  /(1.0 - B - D - F);
					break;
				case Domain::INTERIOR_AIR:
					U[i][j][k] = foundation.indoorAirTemperature;
					break;
				case Domain::INTERIOR_WALL:
					{
					double Tair = foundation.indoorAirTemperature;
					double hc = getConvectionCoeff(TOld[i][j][k],Tair,0.0,1.0,false,PI/2.0);
					double hr = getSimpleInteriorIRCoeff(foundation.wall.interiorEmissivity,
														 TOld[i][j][k],Tair);
					double q = 0;

					U[i][j][k] = (domain.getKXP(i,j,k)*U[i+1][j][k]/domain.getDXP(i) +
							 (hc + hr)*Tair + q)/(domain.getKXP(i,j,k)/domain.getDXP(i) + (hc + hr));
					}
					break;
				case Domain::INTERIOR_SLAB:
					{
					double Tair = foundation.indoorAirTemperature;
					double hc = getConvectionCoeff(TOld[i][j][k],Tair,0.0,1.0,false,0.0);
					double hr = getSimpleInteriorIRCoeff(foundation.wall.interiorEmissivity,
														 TOld[i][j][k],Tair);
					double q = 0;

					U[i][j][k] = (domain.getKZM(i,j,k)*U[i][j][k-1]/domain.getDZM(k) +
							 (hc + hr)*Tair + q)/(domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr));
					}
					break;
				case Domain::EXTERIOR_AIR:
					U[i][j][k] = getOutdoorTemperature();
					break;
				case Domain::EXTERIOR_GRADE:
					{
					double Tair = getOutdoorTemperature();
					double v = getLocalWindSpeed();
					double eSky = weatherData.skyEmissivity.getValue(getSimTime(tNow));
					double F = getEffectiveExteriorViewFactor(eSky,0.0);
					double hc = getConvectionCoeff(TOld[i][j][k],Tair,v,1.0,true,0.0);
					double hr = getExteriorIRCoeff(foundation.soilEmissivity,TOld[i][j][k],Tair,eSky,0.0);
					double q = foundation.soilAbsorptivity*weatherData.globalHorizontalSolar.getValue(getSimTime(tNow));

					U[i][j][k] = (domain.getKZM(i,j,k)*U[i][j][k-1]/domain.getDZM(k) +
							  (hc + hr*pow(F,0.25))*Tair + q)/(domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr));
					}
					break;
				case Domain::EXTERIOR_WALL:
					{
					double Tair = getOutdoorTemperature();
					double v = getLocalWindSpeed();
					double eSky = weatherData.skyEmissivity.getValue(getSimTime(tNow));
					double F = getEffectiveExteriorViewFactor(eSky,PI/2.0);
					double hc = getConvectionCoeff(TOld[i][j][k],Tair,v,1.0,true,PI/2.0);
					double hr = getExteriorIRCoeff(foundation.wall.exteriorEmissivity,TOld[i][j][k],Tair,eSky,PI/2.0);
					double q = foundation.wall.exteriorAbsorptivity*weatherData.globalHorizontalSolar.getValue(getSimTime(tNow));

					U[i][j][k] = (domain.getKXM(i,j,k)*U[i-1][j][k]/domain.getDXM(i) +
							  (hc + hr*pow(F,0.25))*Tair + q)/(domain.getKXM(i,j,k)/domain.getDXM(i) + (hc + hr));
					}
					break;
				case Domain::DEEP_GROUND:
					U[i][j][k] = getDeepGroundTemperature();
					break;
				default:
					U[i][j][k] = (UOld[i][j][k]*(1.0 - A - C - E) - U[i-1][j][k]*(B + D)
							  + UOld[i+1][j][k]*(A + C) - U[i][j][k-1]*F + UOld[i][j][k+1]*E)
							  /(1.0 - B - D - F);
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
				double r = domain.meshX.centers[i];
				double A = domain.cxp_c[i][j][k]*domain.theta[i][j][k]/r;
				double B = domain.cxm_c[i][j][k]*domain.theta[i][j][k]/r;
				double C = domain.cxp[i][j][k]*domain.theta[i][j][k];
				double D = domain.cxm[i][j][k]*domain.theta[i][j][k];
				double E = domain.czp[i][j][k]*domain.theta[i][j][k];
				double F = domain.czm[i][j][k]*domain.theta[i][j][k];

				switch (domain.cellType[i][j][k])
				{
				case Domain::SYMMETRY:
					// Apply zero-flux BC (VOld[i-1][j][k] = V[i+1][j][k])
					// Remove A & B at axis to avoid divide by zero
					V[i][j][k] = (VOld[i][j][k]*(1.0 + D + F) - V[i+1][j][k]*D
							  + V[i+1][j][k]*C - VOld[i][j][k-1]*F + V[i][j][k+1]*E)
							  /(1.0 + C + E);
					break;
				case Domain::FAR_FIELD:
					// Apply zero-flux BC (V[i+1][j][k] = VOld[i-1][j][k])
					V[i][j][k] = (VOld[i][j][k]*(1.0 + B + D + F) - VOld[i-1][j][k]*(B + D)
							  + VOld[i-1][j][k]*(A + C) - VOld[i][j][k-1]*F + V[i][j][k+1]*E)
							  /(1.0 + A + C + E);
					break;
				case Domain::WALL_TOP:
					// Apply zero-flux BC (V[i][j][k+1] = VOld[i][j][k-1])
					V[i][j][k] = (VOld[i][j][k]*(1.0 + B + D + F) - VOld[i-1][j][k]*(B + D)
							  + V[i+1][j][k]*(A + C) - VOld[i][j][k-1]*F + VOld[i][j][k-1]*E)
							  /(1.0 + A + C + E);
					break;
				case Domain::INTERIOR_AIR:
					V[i][j][k] = foundation.indoorAirTemperature;
					break;
				case Domain::INTERIOR_WALL:
					{
					double Tair = foundation.indoorAirTemperature;
					double hc = getConvectionCoeff(TOld[i][j][k],Tair,0.0,1.0,false,PI/2.0);
					double hr = getSimpleInteriorIRCoeff(foundation.wall.interiorEmissivity,
														 TOld[i][j][k],Tair);
					double q = 0;

					V[i][j][k] = (domain.getKXP(i,j,k)*VOld[i+1][j][k]/domain.getDXP(i) +
							  (hc + hr)*Tair + q)/(domain.getKXP(i,j,k)/domain.getDXP(i) + (hc + hr));
					}
					break;
				case Domain::INTERIOR_SLAB:
					{
					double Tair = foundation.indoorAirTemperature;
					double hc = getConvectionCoeff(TOld[i][j][k],Tair,0.0,1.0,false,0.0);
					double hr = getSimpleInteriorIRCoeff(foundation.wall.interiorEmissivity,
														 TOld[i][j][k],Tair);
					double q = 0;

					V[i][j][k] = (domain.getKZM(i,j,k)*VOld[i][j][k-1]/domain.getDZM(k) +
							  (hc + hr)*Tair + q)/(domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr));
					}
					break;
				case Domain::EXTERIOR_AIR:
					V[i][j][k] = getOutdoorTemperature();
					break;
				case Domain::EXTERIOR_GRADE:
					{
					double Tair = getOutdoorTemperature();
					double v = getLocalWindSpeed();
					double eSky = weatherData.skyEmissivity.getValue(getSimTime(tNow));
					double F = getEffectiveExteriorViewFactor(eSky,0.0);
					double hc = getConvectionCoeff(TOld[i][j][k],Tair,v,1.0,true,0.0);
					double hr = getExteriorIRCoeff(foundation.soilEmissivity,TOld[i][j][k],Tair,eSky,0.0);
					double q = foundation.soilAbsorptivity*weatherData.globalHorizontalSolar.getValue(getSimTime(tNow));

					V[i][j][k] = (domain.getKZM(i,j,k)*VOld[i][j][k-1]/domain.getDZM(k) +
							(hc + hr*pow(F,0.25))*Tair + q)/(domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr));
					}
					break;
				case Domain::EXTERIOR_WALL:
					{
					double Tair = getOutdoorTemperature();
					double v = getLocalWindSpeed();
					double eSky = weatherData.skyEmissivity.getValue(getSimTime(tNow));
					double F = getEffectiveExteriorViewFactor(eSky,PI/2.0);
					double hc = getConvectionCoeff(TOld[i][j][k],Tair,v,1.0,true,PI/2.0);
					double hr = getExteriorIRCoeff(foundation.wall.exteriorEmissivity,TOld[i][j][k],Tair,eSky,PI/2.0);
					double q = foundation.wall.exteriorAbsorptivity*weatherData.globalHorizontalSolar.getValue(getSimTime(tNow));

					V[i][j][k] = (domain.getKXM(i,j,k)*VOld[i-1][j][k]/domain.getDXM(i) +
							(hc + hr*pow(F,0.25))*Tair + q)/(domain.getKXM(i,j,k)/domain.getDXM(i) + (hc + hr));
					}
					break;
				case Domain::DEEP_GROUND:
					V[i][j][k] = getDeepGroundTemperature();
					break;
				default:
					V[i][j][k] = (VOld[i][j][k]*(1.0 + B + D + F) - VOld[i-1][j][k]*(B + D)
							  + V[i+1][j][k]*(A + C) - VOld[i][j][k-1]*F + V[i][j][k+1]*E)
							  /(1.0 + A + C + E);
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
				double r = domain.meshX.centers[i];
				double A = domain.cxp_c[i][j][k]*domain.theta[i][j][k]/r;
				double B = domain.cxm_c[i][j][k]*domain.theta[i][j][k]/r;
				double C = domain.cxp[i][j][k]*domain.theta[i][j][k];
				double D = domain.cxm[i][j][k]*domain.theta[i][j][k];
				double E = domain.czp[i][j][k]*domain.theta[i][j][k];
				double F = domain.czm[i][j][k]*domain.theta[i][j][k];

				switch (domain.cellType[i][j][k])
				{
				case Domain::SYMMETRY:
					// Apply zero-flux BC (TOld[i-1][j][k] = TOld[i+1][j][k])
					// Remove A & B at axis to avoid divide by zero
					TNew[i][j][k] = (TOld[i][j][k]*(1.0 + + D + F - C - E) - TOld[i+1][j][k]*(D)
							  + TOld[i+1][j][k]*(C) - TOld[i][j][k-1]*F + TOld[i][j][k+1]*E);
					break;
				case Domain::FAR_FIELD:
					// Apply zero-flux BC (TOld[i+1][j][k] = TOld[i-1][j][k])
					TNew[i][j][k] = (TOld[i][j][k]*(1.0 + B + D + F - A - C - E) - TOld[i-1][j][k]*(B + D)
							  + TOld[i-1][j][k]*(A + C) - TOld[i][j][k-1]*F + TOld[i][j][k+1]*E);
					break;
				case Domain::WALL_TOP:
					// Apply zero-flux BC (TOld[i][j][k+1] = TOld[i][j][k-1]
					TNew[i][j][k] = (TOld[i][j][k]*(1.0 + B + D + F - A - C - E) - TOld[i-1][j][k]*(B + D)
							  + TOld[i+1][j][k]*(A + C) - TOld[i][j][k-1]*F + TOld[i][j][k-1]*E);
					break;
				case Domain::INTERIOR_AIR:
					TNew[i][j][k] = foundation.indoorAirTemperature;
					break;
				case Domain::INTERIOR_SLAB:
					{
					double Tair = foundation.indoorAirTemperature;
					double hc = getConvectionCoeff(TOld[i][j][k],Tair,0.0,1.0,false,0.0);
					double hr = getSimpleInteriorIRCoeff(foundation.wall.interiorEmissivity,
														 TOld[i][j][k],Tair);
					double q = 0;

					TNew[i][j][k] = (domain.getKZM(i,j,k)*TOld[i][j][k-1]/domain.getDZM(k) +
							  (hc + hr)*Tair + q)/(domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr));
					}
					break;
				case Domain::INTERIOR_WALL:
					{
					double Tair = foundation.indoorAirTemperature;
					double hc = getConvectionCoeff(TOld[i][j][k],Tair,0.0,1.0,false,PI/2.0);
					double hr = getSimpleInteriorIRCoeff(foundation.wall.interiorEmissivity,
														 TOld[i][j][k],Tair);
					double q = 0;

					TNew[i][j][k] = (domain.getKXP(i,j,k)*TOld[i+1][j][k]/domain.getDXP(i) +
							  (hc + hr)*Tair + q)/(domain.getKXP(i,j,k)/domain.getDXP(i) + (hc + hr));
					}
					break;
				case Domain::EXTERIOR_AIR:
					TNew[i][j][k] = getOutdoorTemperature();
					break;
				case Domain::EXTERIOR_GRADE:
					{
					double Tair = getOutdoorTemperature();
					double v = getLocalWindSpeed();
					double eSky = weatherData.skyEmissivity.getValue(getSimTime(tNow));
					double F = getEffectiveExteriorViewFactor(eSky,0.0);
					double hc = getConvectionCoeff(TOld[i][j][k],Tair,v,1.0,true,0.0);
					double hr = getExteriorIRCoeff(foundation.soilEmissivity,TOld[i][j][k],Tair,eSky,0.0);
					double q = foundation.soilAbsorptivity*weatherData.globalHorizontalSolar.getValue(getSimTime(tNow));

					TNew[i][j][k] = (domain.getKZM(i,j,k)*TOld[i][j][k-1]/domain.getDZM(k) +
							(hc + hr*pow(F,0.25))*Tair + q)/(domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr));
					}
					break;
				case Domain::EXTERIOR_WALL:
					{
					double Tair = getOutdoorTemperature();
					double v = getLocalWindSpeed();
					double eSky = weatherData.skyEmissivity.getValue(getSimTime(tNow));
					double F = getEffectiveExteriorViewFactor(eSky,PI/2.0);
					double hc = getConvectionCoeff(TOld[i][j][k],Tair,v,1.0,true,PI/2.0);
					double hr = getExteriorIRCoeff(foundation.wall.exteriorEmissivity,TOld[i][j][k],Tair,eSky,PI/2.0);
					double q = foundation.wall.exteriorAbsorptivity*weatherData.globalHorizontalSolar.getValue(getSimTime(tNow));

					TNew[i][j][k] = (domain.getKXM(i,j,k)*TOld[i][j][k]/domain.getDXM(i) +
							(hc + hr*pow(F,0.25))*Tair + q)/(domain.getKXM(i,j,k)/domain.getDXM(i) + (hc + hr));
					}
					break;
				case Domain::DEEP_GROUND:
					TNew[i][j][k] = getDeepGroundTemperature();
					break;
				default:
					TNew[i][j][k] = (TOld[i][j][k]*(1.0 + B + D + F - A - C - E) - TOld[i-1][j][k]*(B + D)
							  + TOld[i+1][j][k]*(A + C) - TOld[i][j][k-1]*F + TOld[i][j][k+1]*E);
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

void Ground::calculateImplicit()
{
	for (size_t k = 0; k < nZ; k++)
	{
		for (size_t j = 0; j < nY; ++j)
		{
			for (size_t i = 0; i < nX; i++)
			{
				double r = domain.meshX.centers[i];
				double A = domain.cxp_c[i][j][k]*domain.theta[i][j][k]/r;
				double B = domain.cxm_c[i][j][k]*domain.theta[i][j][k]/r;
				double C = domain.cxp[i][j][k]*domain.theta[i][j][k];
				double D = domain.cxm[i][j][k]*domain.theta[i][j][k];
				double E = domain.czp[i][j][k]*domain.theta[i][j][k];
				double F = domain.czm[i][j][k]*domain.theta[i][j][k];

				switch (domain.cellType[i][j][k])
				{
				case Domain::SYMMETRY:
					// Apply zero-flux BC (T[i-1][j][k] = T[i+1][j][k])
					// Remove A & B at axis to avoid divide by zero
					Amat(i + nX*k, i + nX*k) = (1.0 + C + E - D - F);
					Amat(i + nX*k, (i+1) + nX*k) = D - C;
					Amat(i + nX*k, i + nX*(k-1)) = F;
					Amat(i + nX*k, i + nX*(k+1)) = -E;

					b(i + nX*k) = TOld[i][j][k];
					break;
				case Domain::FAR_FIELD:
					// Apply zero-flux BC (T[i+1][j][k] = T[i-1][j][k])
					Amat(i + nX*k, i + nX*k) = (1.0 + A + C + E - B - D - F);
					Amat(i + nX*k, (i-1) + nX*k) = (B + D - A - C);
					Amat(i + nX*k, i + nX*(k-1)) = F;
					Amat(i + nX*k, i + nX*(k+1)) = -E;

					b(i + nX*k) = TOld[i][j][k];
					break;
				case Domain::WALL_TOP:
					// Apply zero-flux BC (T[i][j][k+1] = T[i][j][k-1]
					Amat(i + nX*k, i + nX*k) = (1.0 + A + C + E - B - D - F);
					Amat(i + nX*k, (i-1) + nX*k) = (B + D);
					Amat(i + nX*k, (i+1) + nX*k) = (-A - C);
					Amat(i + nX*k, i + nX*(k-1)) = F - E;

					b(i + nX*k) = TOld[i][j][k];
					break;
				case Domain::INTERIOR_AIR:
					Amat(i + nX*k, i + nX*k) = 1.0;

					b(i + nX*k) = foundation.indoorAirTemperature;
					break;
				case Domain::INTERIOR_SLAB:
					{
					double Tair = foundation.indoorAirTemperature;
					double hc = getConvectionCoeff(TOld[i][j][k],Tair,0.0,1.0,false,0.0);
					double hr = getSimpleInteriorIRCoeff(foundation.wall.interiorEmissivity,
														 TOld[i][j][k],Tair);
					double q = 0;

					Amat(i + nX*k, i + nX*k) = domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr);
					Amat(i + nX*k, i + nX*(k-1)) = -domain.getKZM(i,j,k)/domain.getDZM(k);

					b(i + nX*k) = (hc + hr)*Tair + q;
					}
					break;
				case Domain::INTERIOR_WALL:
					{
					double Tair = foundation.indoorAirTemperature;
					double hc = getConvectionCoeff(TOld[i][j][k],Tair,0.0,1.0,false,PI/2.0);
					double hr = getSimpleInteriorIRCoeff(foundation.wall.interiorEmissivity,
														 TOld[i][j][k],Tair);
					double q = 0;

					Amat(i + nX*k, i + nX*k) = domain.getKXP(i,j,k)/domain.getDXP(i) + (hc + hr);
					Amat(i + nX*k, (i+1) + nX*k) = -domain.getKXP(i,j,k)/domain.getDXP(i);

					b(i + nX*k) = (hc + hr)*Tair + q;
					}
					break;
				case Domain::EXTERIOR_AIR:
					Amat(i + nX*k, i + nX*k) = 1.0;

					b(i + nX*k) = getOutdoorTemperature();
					break;
				case Domain::EXTERIOR_GRADE:
					{
					// Apply boundary condition
					double Tair = getOutdoorTemperature();
					double v = getLocalWindSpeed();
					double eSky = weatherData.skyEmissivity.getValue(getSimTime(tNow));
					double F = getEffectiveExteriorViewFactor(eSky,0.0);
					double hc = getConvectionCoeff(TOld[i][j][k],Tair,v,1.0,true,0.0);
					double hr = getExteriorIRCoeff(foundation.soilEmissivity,TOld[i][j][k],Tair,eSky,0.0);
					double q = foundation.soilAbsorptivity*weatherData.globalHorizontalSolar.getValue(getSimTime(tNow));

					Amat(i + nX*k, i + nX*k) = domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr);
					Amat(i + nX*k, i + nX*(k-1)) = -domain.getKZM(i,j,k)/domain.getDZM(k);

					b(i + nX*k) = (hc + hr*pow(F,0.25))*Tair + q;
					}
					break;
				case Domain::EXTERIOR_WALL:
					{
					// Apply boundary condition
					double Tair = getOutdoorTemperature();
					double v = getLocalWindSpeed();
					double eSky = weatherData.skyEmissivity.getValue(getSimTime(tNow));
					double F = getEffectiveExteriorViewFactor(eSky,PI/2.0);
					double hc = getConvectionCoeff(TOld[i][j][k],Tair,v,1.0,true,PI/2.0);
					double hr = getExteriorIRCoeff(foundation.wall.exteriorEmissivity,TOld[i][j][k],Tair,eSky,PI/2.0);
					double q = foundation.wall.exteriorAbsorptivity*weatherData.globalHorizontalSolar.getValue(getSimTime(tNow));

					Amat(i + nX*k, i + nX*k) = domain.getKXM(i,j,k)/domain.getDXM(i) + (hc + hr);
					Amat(i + nX*k, ((i-1) + nX*k)) = -domain.getKXM(i,j,k)/domain.getDXM(i);

					b(i + nX*k) = (hc + hr*pow(F,0.25))*Tair + q;
					}
					break;
				case Domain::DEEP_GROUND:
					Amat(i + nX*k, i + nX*k) = 1.0;

					b(i + nX*k) = getDeepGroundTemperature();
					break;
				default:
					Amat(i + nX*k, i + nX*k) = (1.0 + A + C + E - B - D - F);
					Amat(i + nX*k, (i-1) + nX*k) = (B + D);
					Amat(i + nX*k, (i+1) + nX*k) = (-A - C);
					Amat(i + nX*k, i + nX*(k-1)) = F;
					Amat(i + nX*k, i + nX*(k+1)) = -E;

					b(i + nX*k) = TOld[i][j][k];
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
				TNew[i][j][k] = x(i + nX*k);

				// Update old values for next timestep
				TOld[i][j][k] = TNew[i][j][k];
			}
		}
	}

}

void Ground::calculateCrankNicolson()
{
	for (size_t k = 0; k < nZ; k++)
	{
		for (size_t j = 0; j < nY; ++j)
		{
			for (size_t i = 0; i < nX; i++)
			{
				double r = domain.meshX.centers[i];
				double A = 0.5*domain.cxp_c[i][j][k]*domain.theta[i][j][k]/r;
				double B = 0.5*domain.cxm_c[i][j][k]*domain.theta[i][j][k]/r;
				double C = 0.5*domain.cxp[i][j][k]*domain.theta[i][j][k];
				double D = 0.5*domain.cxm[i][j][k]*domain.theta[i][j][k];
				double E = 0.5*domain.czp[i][j][k]*domain.theta[i][j][k];
				double F = 0.5*domain.czm[i][j][k]*domain.theta[i][j][k];

				switch (domain.cellType[i][j][k])
				{
				case Domain::SYMMETRY:
					// Apply zero-flux BC (T[i-1][j][k] = T[i+1][j][k])
					// Remove A & B at axis to avoid divide by zero
					Amat(i + nX*k, i + nX*k) = (1.0 + C + E - D - F);
					Amat(i + nX*k, (i+1) + nX*k) = D - C;
					Amat(i + nX*k, i + nX*(k-1)) = F;
					Amat(i + nX*k, i + nX*(k+1)) = -E;

					b(i + nX*k) = (TOld[i][j][k]*(1.0 + D + F - C - E) - TOld[i+1][j][k]*(D)
							  + TOld[i+1][j][k]*(C) - TOld[i][j][k-1]*F + TOld[i][j][k+1]*E);
					break;
				case Domain::FAR_FIELD:
					// Apply zero-flux BC (T[i+1][j][k] = T[i-1][j][k])
					Amat(i + nX*k, i + nX*k) = (1.0 + A + C + E - B - D - F);
					Amat(i + nX*k, (i-1) + nX*k) = (B + D - A - C);
					Amat(i + nX*k, i + nX*(k-1)) = F;
					Amat(i + nX*k, i + nX*(k+1)) = -E;

					b(i + nX*k) = (TOld[i][j][k]*(1.0 + B + D + F - A - C - E) - TOld[i-1][j][k]*(B + D)
							  + TOld[i-1][j][k]*(A + C) - TOld[i][j][k-1]*F + TOld[i][j][k+1]*E);
					break;
				case Domain::WALL_TOP:
					// Apply zero-flux BC (T[i][j][k+1] = T[i][j][k-1]
					Amat(i + nX*k, i + nX*k) = (1.0 + A + C + E - B - D - F);
					Amat(i + nX*k, (i-1) + nX*k) = (B + D);
					Amat(i + nX*k, (i+1) + nX*k) = (-A - C);
					Amat(i + nX*k, i + nX*(k-1)) = F - E;

					b(i + nX*k) = (TOld[i][j][k]*(1.0 + B + D + F - A - C - E) - TOld[i-1][j][k]*(B + D)
							  + TOld[i+1][j][k]*(A + C) - TOld[i][j][k-1]*F + TOld[i][j][k-1]*E);
					break;
				case Domain::INTERIOR_AIR:
					Amat(i + nX*k, i + nX*k) = 1.0;

					b(i + nX*k) = foundation.indoorAirTemperature;
					break;
				case Domain::INTERIOR_WALL:
					{
					double Tair = foundation.indoorAirTemperature;
					double hc = getConvectionCoeff(TOld[i][j][k],Tair,0.0,1.0,false,PI/2.0);
					double hr = getSimpleInteriorIRCoeff(foundation.wall.interiorEmissivity,
														 TOld[i][j][k],Tair);
					double q = 0;

					Amat(i + nX*k, i + nX*k) = domain.getKXP(i,j,k)/domain.getDXP(i) + (hc + hr);
					Amat(i + nX*k, (i+1) + nX*k) = -domain.getKXP(i,j,k)/domain.getDXP(i);

					b(i + nX*k) = (hc + hr)*Tair + q;
					}
					break;
				case Domain::INTERIOR_SLAB:
					{
					double Tair = foundation.indoorAirTemperature;
					double hc = getConvectionCoeff(TOld[i][j][k],Tair,0.0,1.0,false,0.0);
					double hr = getSimpleInteriorIRCoeff(foundation.wall.interiorEmissivity,
														 TOld[i][j][k],Tair);
					double q = 0;

					Amat(i + nX*k, i + nX*k) = domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr);
					Amat(i + nX*k, i + nX*(k-1)) = -domain.getKZM(i,j,k)/domain.getDZM(k);

					b(i + nX*k) = (hc + hr)*Tair + q;
					}
					break;
				case Domain::EXTERIOR_AIR:
					Amat(i + nX*k, i + nX*k) = 1.0;

					b(i + nX*k) = getOutdoorTemperature();
					break;
				case Domain::EXTERIOR_GRADE:
					{
					// Apply boundary condition
					double Tair = getOutdoorTemperature();
					double v = getLocalWindSpeed();
					double eSky = weatherData.skyEmissivity.getValue(getSimTime(tNow));
					double F = getEffectiveExteriorViewFactor(eSky,0.0);
					double hc = getConvectionCoeff(TOld[i][j][k],Tair,v,1.0,true,0.0);
					double hr = getExteriorIRCoeff(foundation.soilEmissivity,TOld[i][j][k],Tair,eSky,0.0);
					double q = foundation.soilAbsorptivity*weatherData.globalHorizontalSolar.getValue(getSimTime(tNow));

					Amat(i + nX*k, i + nX*k) = domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr);
					Amat(i + nX*k, i + nX*(k-1)) = -domain.getKZM(i,j,k)/domain.getDZM(k);

					b(i + nX*k) = (hc + hr*pow(F,0.25))*Tair + q;
					}
					break;
				case Domain::DEEP_GROUND:
					Amat(i + nX*k, i + nX*k) = 1.0;

					b(i + nX*k) = getDeepGroundTemperature();
					break;
				case Domain::EXTERIOR_WALL:
					{
					// Apply boundary condition
					double Tair = getOutdoorTemperature();
					double v = getLocalWindSpeed();
					double eSky = weatherData.skyEmissivity.getValue(getSimTime(tNow));
					double F = getEffectiveExteriorViewFactor(eSky,PI/2.0);
					double hc = getConvectionCoeff(TOld[i][j][k],Tair,v,1.0,true,PI/2.0);
					double hr = getExteriorIRCoeff(foundation.wall.exteriorEmissivity,TOld[i][j][k],Tair,eSky,PI/2.0);
					double q = foundation.wall.exteriorAbsorptivity*weatherData.globalHorizontalSolar.getValue(getSimTime(tNow));

					Amat(i + nX*k, i + nX*k) = domain.getKXM(i,j,k)/domain.getDXM(i) + (hc + hr);
					Amat(i + nX*k, ((i-1) + nX*k)) = -domain.getKXM(i,j,k)/domain.getDXM(i);

					b(i + nX*k) = (hc + hr*pow(F,0.25))*Tair + q;
					}
					break;
				default:
					Amat(i + nX*k, i + nX*k) = (1.0 + A + C + E - B - D - F);
					Amat(i + nX*k, (i-1) + nX*k) = (B + D);
					Amat(i + nX*k, (i+1) + nX*k) = (-A - C);
					Amat(i + nX*k, i + nX*(k-1)) = F;
					Amat(i + nX*k, i + nX*(k+1)) = -E;

					b(i + nX*k) = (TOld[i][j][k]*(1.0 + B + D + F - A - C - E) - TOld[i-1][j][k]*(B + D)
							  + TOld[i+1][j][k]*(A + C) - TOld[i][j][k-1]*F + TOld[i][j][k+1]*E);
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
				TNew[i][j][k] = x(i + nX*k);

				// Update old values for next timestep
				TOld[i][j][k] = TNew[i][j][k];
			}
		}
	}

}

void Ground::calculateSteadyState()
{
	for (size_t k = 0; k < nZ; k++)
	{
		for (size_t j = 0; j < nY; ++j)
		{
			for (size_t i = 0; i < nX; i++)
			{
				double r = domain.meshX.centers[i];
				double A = domain.cxp_c[i][j][k]/r;
				double B = domain.cxm_c[i][j][k]/r;
				double C = domain.cxp[i][j][k];
				double D = domain.cxm[i][j][k];
				double E = domain.czp[i][j][k];
				double F = domain.czm[i][j][k];

				switch (domain.cellType[i][j][k])
				{
				case Domain::SYMMETRY:
					// Apply zero-flux BC (T[i-1][j][k] = T[i+1][j][k])
					// Remove A & B at axis to avoid divide by zero
					Amat(i + nX*k, i + nX*k) = (D + F - C - E);
					Amat(i + nX*k, (i+1) + nX*k) = C - D;
					Amat(i + nX*k, i + nX*(k-1)) = -F;
					Amat(i + nX*k, i + nX*(k+1)) = E;

					b(i + nX*k) = 0;
					break;
				case Domain::FAR_FIELD:
					// Apply zero-flux BC (T[i+1][j][k] = T[i-1][j][k])
					Amat(i + nX*k, i + nX*k) = (B + D + F - A - C - E);
					Amat(i + nX*k, (i-1) + nX*k) = (-B - D) + (A + C);
					Amat(i + nX*k, i + nX*(k-1)) = -F;
					Amat(i + nX*k, i + nX*(k+1)) = E;

					b(i + nX*k) = 0;
					break;
				case Domain::WALL_TOP:
					// Apply zero-flux BC (T[i][j][k+1] = T[i][j][k-1]
					Amat(i + nX*k, i + nX*k) = (B + D + F - A - C - E);
					Amat(i + nX*k, (i-1) + nX*k) = (-B - D);
					Amat(i + nX*k, (i+1) + nX*k) = (A + C);
					Amat(i + nX*k, i + nX*(k-1)) = E - F;

					b(i + nX*k) = 0;
					break;
				case Domain::INTERIOR_AIR:
					Amat(i + nX*k, i + nX*k) = 1.0;

					b(i + nX*k) = foundation.indoorAirTemperature;
					break;
				case Domain::INTERIOR_SLAB:
					{
					double Tair = foundation.indoorAirTemperature;
					double hc = getConvectionCoeff(TOld[i][j][k],Tair,0.0,1.0,false,0.0);
					double hr = getSimpleInteriorIRCoeff(foundation.wall.interiorEmissivity,
														 TOld[i][j][k],Tair);
					double q = 0;

					Amat(i + nX*k, i + nX*k) = domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr);
					Amat(i + nX*k, i + nX*(k-1)) = -domain.getKZM(i,j,k)/domain.getDZM(k);

					b(i + nX*k) = (hc + hr)*Tair + q;
					}
					break;
				case Domain::INTERIOR_WALL:
					{
					double Tair = foundation.indoorAirTemperature;
					double hc = getConvectionCoeff(TOld[i][j][k],Tair,0.0,1.0,false,PI/2.0);
					double hr = getSimpleInteriorIRCoeff(foundation.wall.interiorEmissivity,
														 TOld[i][j][k],Tair);
					double q = 0;

					Amat(i + nX*k, i + nX*k) = domain.getKXP(i,j,k)/domain.getDXP(i) + (hc + hr);
					Amat(i + nX*k, (i+1) + nX*k) = -domain.getKXP(i,j,k)/domain.getDXP(i);

					b(i + nX*k) = (hc + hr)*Tair + q;
					}
					break;
				case Domain::EXTERIOR_AIR:
					Amat(i + nX*k, i + nX*k) = 1.0;

					b(i + nX*k) = getOutdoorTemperature();
					break;
				case Domain::EXTERIOR_GRADE:
					{
					// Apply boundary condition
					double Tair = getOutdoorTemperature();
					double v = getLocalWindSpeed();
					double eSky = weatherData.skyEmissivity.getValue(getSimTime(tNow));
					double F = getEffectiveExteriorViewFactor(eSky,0.0);
					double hc = getConvectionCoeff(TOld[i][j][k],Tair,v,1.0,true,0.0);
					double hr = getExteriorIRCoeff(foundation.soilEmissivity,TOld[i][j][k],Tair,eSky,0.0);
					double q = foundation.soilAbsorptivity*weatherData.globalHorizontalSolar.getValue(getSimTime(tNow));

					Amat(i + nX*k, i + nX*k) = domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr);
					Amat(i + nX*k, i + nX*(k-1)) = -domain.getKZM(i,j,k)/domain.getDZM(k);

					b(i + nX*k) = (hc + hr*pow(F,0.25))*Tair + q;
					}
					break;
				case Domain::EXTERIOR_WALL:
					{
					double Tair = getOutdoorTemperature();
					double v = getLocalWindSpeed();
					double eSky = weatherData.skyEmissivity.getValue(getSimTime(tNow));
					double F = getEffectiveExteriorViewFactor(eSky,PI/2.0);
					double hc = getConvectionCoeff(TOld[i][j][k],Tair,v,1.0,true,PI/2.0);
					double hr = getExteriorIRCoeff(foundation.wall.exteriorEmissivity,TOld[i][j][k],Tair,eSky,PI/2.0);
					double q = foundation.wall.exteriorAbsorptivity*weatherData.globalHorizontalSolar.getValue(getSimTime(tNow));

					Amat(i + nX*k, i + nX*k) = domain.getKXM(i,j,k)/domain.getDXM(i) + (hc + hr);
					Amat(i + nX*k, ((i-1) + nX*k)) = -domain.getKXM(i,j,k)/domain.getDXM(i);

					b(i + nX*k) = (hc + hr*pow(F,0.25))*Tair + q;
					}
					break;
				case Domain::DEEP_GROUND:
					Amat(i + nX*k, i + nX*k) = 1.0;

					b(i + nX*k) = getDeepGroundTemperature();
					break;
				default:
					Amat(i + nX*k, i + nX*k) = (B + D + F - A - C - E);
					Amat(i + nX*k, (i-1) + nX*k) = (-B - D);
					Amat(i + nX*k, (i+1) + nX*k) = (A + C);
					Amat(i + nX*k, i + nX*(k-1)) = -F;
					Amat(i + nX*k, i + nX*(k+1)) = E;

					b(i + nX*k) = 0;
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
				TNew[i][j][k] = x(i + nX*k);

				// Update old values for next timestep
				TOld[i][j][k] = TNew[i][j][k];
			}
		}
	}
}

void Ground::calculate(double t)
{
	tNow = t;

	// Calculate Temperatures
	if (foundation.numericalScheme == Foundation::NS_ADE)
		calculateADE();
	else if (foundation.numericalScheme == Foundation::NS_EXPLICIT)
		calculateExplicit();
	else if (foundation.numericalScheme == Foundation::NS_IMPLICIT)
		calculateImplicit();
	else if (foundation.numericalScheme == Foundation::NS_CRANK_NICOLSON)
		calculateCrankNicolson();
	else if (foundation.numericalScheme == Foundation::NS_STEADY_STATE)
		calculateSteadyState();

	// Calculate Heat Fluxes

	// Slab
	size_t j = 0;
	for (size_t i = domain.slabImin; i <= domain.slabImax; i++)
	{
		double Tair = foundation.indoorAirTemperature;
		double h = getConvectionCoeff(TNew[i][j][domain.slabK],Tair,0.0,1.0,false,0.0);
		double A = 2.0*PI*domain.meshX.deltas[i]*domain.meshX.centers[i]; // for cylindrical coordinates
		//double A = domain.mesher.xDeltas[i];  // for cartesian coordinates
		QSlab[i] = h*A*(Tair - TNew[i][j][domain.slabK]);
	}

	QSlabTotal = accumulate((QSlab).begin(),(QSlab).end(), 0.0)/(PI*pow(foundation.effectiveLength,2.0));  // for cylindrical coordinates
	//QSlabTotal = accumulate((QSlab).begin(),(QSlab).end(), 0.0)/foundation.radius;  // for cartesian coordinates

	if (makePlot)
		plot();
}

void Ground::plot()
{

	if (tNow >= nextPlotTime)
	{
		for(size_t k = 0; k < nZ; k++)
		{
			for(size_t i = 0; i < nX; i++)
			{
				size_t j = 0;
			    double TinF = (TNew[i][j][k] - 273.15)*9/5 + 32.0;
				TDat.a[i+nX*k] = TinF;
			}

		}

		int nX = rDat.GetNN();
		double rmin = 0.0;
		double rmax = rGrid.a[nX];
		double rrange = rmax - rmin;
		double length = foundation.effectiveLength;
		mglData rTicks(2);
		rTicks.a[0] = length;
		rTicks.a[1] = rmax;

		std::string sLength = str(boost::format("%0.2f") % length) + " m";
		std::string sRmax = str(boost::format("%0.2f") % rmax) + " m";
		std::string rTickString = sLength + "\n" + sRmax;

		int nZ = zDat.GetNN();
		double zmin = zGrid.a[0];
		double zmax = zGrid.a[nZ];
		double zrange = zmax - zmin;
		mglData zTicks(2);
		zTicks.a[0] = zmin;
		zTicks.a[1] = 0.0;

		std::string sZmin = str(boost::format("%0.2f") % zmin) + " m";
		std::string zTickString = sZmin + "\n 0,0\n";

		int nT = cDat.GetNN();
		double Tmin = cDat.a[0];
		double Tmax = cDat.a[nT - 1];
		double Tstep = cDat.a[1] - cDat.a[0];


		gr.NewFrame();

		// Plot
		//gr.MultiPlot(3,1,0,2,1,"_");
		gr.LoadFont("none");
		gr.SetFontSize(2.0);
		gr.SetOrigin(0.0, zmax);
		gr.SetRange('x', rGrid);
		gr.SetRange('y', zGrid);
		gr.SetRange('c', Tmin, Tmax);
		gr.SetRange('z', Tmin, Tmax);
		gr.SetTicks('c', Tstep, nT, Tmin);
		gr.SetTicksVal('x', rTicks, rTickString.c_str());
		gr.SetTicksVal('y', zTicks, zTickString.c_str());
		gr.SetTickLen(-0.0001);
		gr.Aspect(rrange, zrange);
		gr.Axis("yU");
		gr.Axis("x");
		gr.Colorbar("_");
		//gr.Puts(mglPoint(12.5,-10), "Temperature \\textdegree C", ":C");
		gr.Box("k",false);
		gr.Dens(rDat, zDat, TDat);
		if (foundation.outputAnimation.contours)
			gr.Cont(cDat, rDat, zDat, TDat,"H");
		if (foundation.outputAnimation.gradients)
			gr.Grad(rDat, zDat, TDat);
		if (foundation.outputAnimation.grid)
			gr.Grid(rGrid, zGrid, TGrid, "W");

		// Draw blocks
		for (size_t b = 0; b < foundation.blocks.size(); b++)
		{
			mglPoint bl = mglPoint(foundation.blocks[b].xMin,
			         	 	 	   foundation.blocks[b].zMin,
			         	 	 	   210.0);
			mglPoint br = mglPoint(foundation.blocks[b].xMax,
			         	 	 	   foundation.blocks[b].zMin,
			         	 	 	   210.0);
			mglPoint tr = mglPoint(foundation.blocks[b].xMax,
			         	 	 	   foundation.blocks[b].zMax,
			         	 	 	   210.0);
			mglPoint tl = mglPoint(foundation.blocks[b].xMin,
			         	 	 	   foundation.blocks[b].zMax,
			         	 	 	   210.0);

			gr.Line(bl, br, "k");
			gr.Line(br, tr, "k");
			gr.Line(tr, tl, "k");
			gr.Line(tl, bl, "k");

		}

		// Timestamp
		boost::posix_time::ptime tp(simulationControl.startDate,boost::posix_time::seconds(tNow));
		gr.Puts(mglPoint(-0.10,0.40), to_simple_string(tp).c_str(), ":C");

		/*
		// Text
		double y = 1.0;
		double d = 0.05;
		double x = -0.2;
		gr.SubPlot(3,1,2,"");
		gr.LoadFont("none");
		gr.SetFontSizePT(10.0);
		gr.SetRanges(0,1,0,1);

		// Grid information
		string nXs = boost::lexical_cast<string>(nX);
		string nZs = boost::lexical_cast<string>(nZ);
		string nCells = boost::lexical_cast<string>(nX*nZ);


		string gridInfo = "Number of Cells:\n\n\t(" + nXs + " x " + nZs + ") = " +
						  nCells + " Cells";

		gr.Puts(mglPoint(x,y), gridInfo.c_str(), ":L");

		// Location and time stamp
		string location = weatherData.city + ", " + weatherData.state + ", " +
				          weatherData.country;
		ptime tp(simulationControl.startDate,seconds(tNow));

		gr.SetFontSizePT(10.0);
		gr.Puts(mglPoint(-1.25,0.05), location.c_str(), ":C");
		gr.Puts(mglPoint(-1.25,0.05-d), to_simple_string(tp).c_str(), ":C");

		// Footer
		gr.SetFontSizePT(6.0);

		gr.Puts(mglPoint(1.05, -0.05), "Copyright \\textcopyright{ }Big Ladder Software", ":R");
		gr.Puts(mglPoint(1.05, -0.08), startString.c_str(), ":R");
		*/

		//gr.WriteFrame((foundation.outputAnimation.name + ".png").c_str());
		gr.EndFrame();
		nextPlotTime += plotFreq;
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

#endif

