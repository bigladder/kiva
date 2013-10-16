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

#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include "Ground.h"

double getMatrixValue(blas::matrix<double> const &m, size_t i, size_t k)
{
	return m(i,k);
};


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
		U = blas::matrix<double>(nX, nZ);
		UOld = blas::matrix<double>(nX, nZ);

		V = blas::matrix<double>(nX, nZ);
		VOld = blas::matrix<double>(nX, nZ);
	}

	TNew = blas::matrix<double>(nX, nZ);
	TOld = blas::matrix<double>(nX, nZ);
	QSlab = blas::vector<double>(domain.slabImax - domain.slabImin + 1);

	if (foundation.initializationMethod == Foundation::IM_STEADY_STATE)
		calculateSteadyState();

	else if (foundation.initializationMethod == Foundation::IM_IMPLICIT_ACCELERATION)
	{
		Foundation initialFoundation = foundation;
		initialFoundation.initializationMethod = Foundation::IM_STEADY_STATE;
		initialFoundation.numericalScheme = Foundation::NS_IMPLICIT;
		initialFoundation.outputAnimation.name = "";
		SimulationControl initialSimControl;
		initialSimControl.timestep = hours(foundation.implicitAccelTimestep);
		initialSimControl.endDate = simulationControl.startDate;
		ptime endTime(initialSimControl.endDate);
		endTime = endTime - simulationControl.timestep;
		initialSimControl.startTime = endTime - hours(foundation.implicitAccelTimestep*foundation.implicitAccelPeriods);
		initialSimControl.startDate = initialSimControl.startTime.date();

		ptime simEnd(endTime);
		ptime simStart(initialSimControl.startTime);
		time_duration simDuration =  simEnd  - simStart;

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
		// TODO: remove to incorporate full 3D
		size_t j = 0;
		for (size_t k = 0; k < nZ; ++k)
		{
			for (size_t i = 0; i < nX; ++i)
			{
				switch (domain.cellType[i][j][k])
				{
				case Domain::INTERIOR_AIR:
					TOld(i,k) = foundation.indoorAirTemperature;
					break;
				case Domain::EXTERIOR_AIR:
					TOld(i,k) = getOutdoorTemperature();
					break;
				case Domain::DEEP_GROUND:
					TOld(i,k) = getDeepGroundTemperature();
					break;
				default:
					TOld(i,k)= getInitialTemperature(tNow,
													domain.meshZ.centers[k]);
					break;
				}
			}
		}
	}
}

void Ground::initializePlot()
{

	startString = to_simple_string(second_clock::local_time());

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
		for (size_t i = 0; i < nX; ++i)
		{
			UOld(i,k) = VOld(i,k) = TOld(i,k);
		}
	}

	// Solve for new values (Main loop)
	boost::thread up(boost::bind(&Ground::calculateADEUpwardSweep, this));
	boost::thread down(boost::bind(&Ground::calculateADEDownwardSweep, this));

	up.join();
	down.join();

	for (size_t k = 0; k < nZ; ++k)
	{
		for (size_t i = 0; i < nX; ++i)
		{
			// Calculate average of sweeps
			TNew(i,k) = 0.5*(U(i,k) + V(i,k));
			// Update old values for next timestep
			TOld(i,k) = TNew(i,k);
		}
	}
}

void Ground::calculateADEUpwardSweep()
{

	// TODO: remove to incorporate full 3D
	size_t j = 0;
	// Upward sweep (Solve U Matrix starting from 1, 1)
	for (size_t k = 0; k < nZ; k++)
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
				// Apply zero-flux BC (U(i-1,k) = UOld(i+1,k))
				// Remove A & B at axis to avoid divide by zero
				U(i, k) = (UOld(i,k)*(1.0 - C - E) - UOld(i+1,k)*D
						  + UOld(i+1,k)*C - U(i,k-1)*F + UOld(i,k+1)*E)
						  /(1.0 - D - F);
				break;
			case Domain::FAR_FIELD:
				// Apply zero-flux BC (UOld(i+1,k) = U(i-1,k))
				U(i, k) = (UOld(i,k)*(1.0 - A - C - E) - U(i-1,k)*(B + D)
						  + U(i-1,k)*(A + C) - U(i,k-1)*F + UOld(i,k+1)*E)
						  /(1.0 - B - D - F);
				break;
			case Domain::WALL_TOP:
				// Apply zero-flux BC (UOld(i,k+1) = U(i,k-1)
				U(i, k) = (UOld(i,k)*(1.0 - A - C - E) - U(i-1,k)*(B + D)
						  + UOld(i+1,k)*(A + C) - U(i,k-1)*F + U(i,k-1)*E)
						  /(1.0 - B - D - F);
				break;
			case Domain::INTERIOR_AIR:
				U(i,k) = foundation.indoorAirTemperature;
				break;
			case Domain::INTERIOR_WALL:
				{
				double Tair = foundation.indoorAirTemperature;
				double hc = getConvectionCoeff(TOld(i,k),Tair,0.0,1.0,false,PI/2.0);
				double hr = getSimpleInteriorIRCoeff(foundation.wall.interiorEmissivity,
						                             TOld(i,k),Tair);
				double q = 0;

				U(i,k) = (domain.getKXP(i,j,k)*U(i+1,k)/domain.getDXP(i) +
						 (hc + hr)*Tair + q)/(domain.getKXP(i,j,k)/domain.getDXP(i) + (hc + hr));
				}
				break;
			case Domain::INTERIOR_SLAB:
				{
				double Tair = foundation.indoorAirTemperature;
				double hc = getConvectionCoeff(TOld(i,k),Tair,0.0,1.0,false,0.0);
				double hr = getSimpleInteriorIRCoeff(foundation.wall.interiorEmissivity,
						                             TOld(i,k),Tair);
				double q = 0;

				U(i,k) = (domain.getKZM(i,j,k)*U(i,k-1)/domain.getDZM(k) +
						 (hc + hr)*Tair + q)/(domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr));
				}
				break;
			case Domain::EXTERIOR_AIR:
				U(i,k) = getOutdoorTemperature();
				break;
			case Domain::EXTERIOR_GRADE:
				{
				double Tair = getOutdoorTemperature();
				double v = getLocalWindSpeed();
				double eSky = weatherData.skyEmissivity.getValue(getSimTime(tNow));
				double F = getEffectiveExteriorViewFactor(eSky,0.0);
				double hc = getConvectionCoeff(TOld(i,k),Tair,v,1.0,true,0.0);
				double hr = getExteriorIRCoeff(foundation.soilEmissivity,TOld(i,k),Tair,eSky,0.0);
				double q = foundation.soilAbsorptivity*weatherData.globalHorizontalSolar.getValue(getSimTime(tNow));

				U(i,k) = (domain.getKZM(i,j,k)*U(i,k-1)/domain.getDZM(k) +
						  (hc + hr*pow(F,0.25))*Tair + q)/(domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr));
				}
				break;
			case Domain::EXTERIOR_WALL:
				{
				double Tair = getOutdoorTemperature();
				double v = getLocalWindSpeed();
				double eSky = weatherData.skyEmissivity.getValue(getSimTime(tNow));
				double F = getEffectiveExteriorViewFactor(eSky,PI/2.0);
				double hc = getConvectionCoeff(TOld(i,k),Tair,v,1.0,true,PI/2.0);
				double hr = getExteriorIRCoeff(foundation.wall.exteriorEmissivity,TOld(i,k),Tair,eSky,PI/2.0);
				double q = foundation.wall.exteriorAbsorptivity*weatherData.globalHorizontalSolar.getValue(getSimTime(tNow));

				U(i,k) = (domain.getKXM(i,j,k)*U(i-1,k)/domain.getDXM(i) +
						  (hc + hr*pow(F,0.25))*Tair + q)/(domain.getKXM(i,j,k)/domain.getDXM(i) + (hc + hr));
				}
				break;
			case Domain::DEEP_GROUND:
				U(i,k) = getDeepGroundTemperature();
				break;
			default:
				U(i,k) = (UOld(i,k)*(1.0 - A - C - E) - U(i-1,k)*(B + D)
						  + UOld(i+1,k)*(A + C) - U(i,k-1)*F + UOld(i,k+1)*E)
						  /(1.0 - B - D - F);
				break;
			}
		}
	}
}

void Ground::calculateADEDownwardSweep()
{

	// TODO: remove to incorporate full 3D
	size_t j = 0;

	// Downward sweep (Solve V Matrix starting from I, K)
	for (size_t k = nZ - 1; k >= 0 && k < nZ; k--)
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
				// Apply zero-flux BC (VOld(i-1,k) = V(i+1,k))
				// Remove A & B at axis to avoid divide by zero
				V(i,k) = (VOld(i,k)*(1.0 + D + F) - V(i+1,k)*D
						  + V(i+1,k)*C - VOld(i,k-1)*F + V(i,k+1)*E)
						  /(1.0 + C + E);
				break;
			case Domain::FAR_FIELD:
				// Apply zero-flux BC (V(i+1,k) = VOld(i-1,k))
				V(i,k) = (VOld(i,k)*(1.0 + B + D + F) - VOld(i-1,k)*(B + D)
						  + VOld(i-1,k)*(A + C) - VOld(i,k-1)*F + V(i,k+1)*E)
						  /(1.0 + A + C + E);
				break;
			case Domain::WALL_TOP:
				// Apply zero-flux BC (V(i,k+1) = VOld(i,k-1))
				V(i,k) = (VOld(i,k)*(1.0 + B + D + F) - VOld(i-1,k)*(B + D)
						  + V(i+1,k)*(A + C) - VOld(i,k-1)*F + VOld(i,k-1)*E)
						  /(1.0 + A + C + E);
				break;
			case Domain::INTERIOR_AIR:
				V(i,k) = foundation.indoorAirTemperature;
				break;
			case Domain::INTERIOR_WALL:
				{
				double Tair = foundation.indoorAirTemperature;
				double hc = getConvectionCoeff(TOld(i,k),Tair,0.0,1.0,false,PI/2.0);
				double hr = getSimpleInteriorIRCoeff(foundation.wall.interiorEmissivity,
						                             TOld(i,k),Tair);
				double q = 0;

				V(i,k) = (domain.getKXP(i,j,k)*VOld(i+1,k)/domain.getDXP(i) +
						  (hc + hr)*Tair + q)/(domain.getKXP(i,j,k)/domain.getDXP(i) + (hc + hr));
				}
				break;
			case Domain::INTERIOR_SLAB:
				{
				double Tair = foundation.indoorAirTemperature;
				double hc = getConvectionCoeff(TOld(i,k),Tair,0.0,1.0,false,0.0);
				double hr = getSimpleInteriorIRCoeff(foundation.wall.interiorEmissivity,
						                             TOld(i,k),Tair);
				double q = 0;

				V(i,k) = (domain.getKZM(i,j,k)*VOld(i,k-1)/domain.getDZM(k) +
						  (hc + hr)*Tair + q)/(domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr));
				}
				break;
			case Domain::EXTERIOR_AIR:
				V(i,k) = getOutdoorTemperature();
				break;
			case Domain::EXTERIOR_GRADE:
				{
				double Tair = getOutdoorTemperature();
				double v = getLocalWindSpeed();
				double eSky = weatherData.skyEmissivity.getValue(getSimTime(tNow));
				double F = getEffectiveExteriorViewFactor(eSky,0.0);
				double hc = getConvectionCoeff(TOld(i,k),Tair,v,1.0,true,0.0);
				double hr = getExteriorIRCoeff(foundation.soilEmissivity,TOld(i,k),Tair,eSky,0.0);
				double q = foundation.soilAbsorptivity*weatherData.globalHorizontalSolar.getValue(getSimTime(tNow));

				V(i,k) = (domain.getKZM(i,j,k)*VOld(i,k-1)/domain.getDZM(k) +
						(hc + hr*pow(F,0.25))*Tair + q)/(domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr));
				}
				break;
			case Domain::EXTERIOR_WALL:
				{
				double Tair = getOutdoorTemperature();
				double v = getLocalWindSpeed();
				double eSky = weatherData.skyEmissivity.getValue(getSimTime(tNow));
				double F = getEffectiveExteriorViewFactor(eSky,PI/2.0);
				double hc = getConvectionCoeff(TOld(i,k),Tair,v,1.0,true,PI/2.0);
				double hr = getExteriorIRCoeff(foundation.wall.exteriorEmissivity,TOld(i,k),Tair,eSky,PI/2.0);
				double q = foundation.wall.exteriorAbsorptivity*weatherData.globalHorizontalSolar.getValue(getSimTime(tNow));

				V(i,k) = (domain.getKXM(i,j,k)*VOld(i-1,k)/domain.getDXM(i) +
						(hc + hr*pow(F,0.25))*Tair + q)/(domain.getKXM(i,j,k)/domain.getDXM(i) + (hc + hr));
				}
				break;
			case Domain::DEEP_GROUND:
				V(i,k) = getDeepGroundTemperature();
				break;
			default:
				V(i,k) = (VOld(i,k)*(1.0 + B + D + F) - VOld(i-1,k)*(B + D)
						  + V(i+1,k)*(A + C) - VOld(i,k-1)*F + V(i,k+1)*E)
						  /(1.0 + A + C + E);
				break;
			}
		}
	}
}

void Ground::calculateExplicit()
{

	// TODO: remove to incorporate full 3D
	size_t j = 0;
	for (size_t k = 0; k < nZ; k++)
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
				// Apply zero-flux BC (TOld(i-1,k) = TOld(i+1,k))
				// Remove A & B at axis to avoid divide by zero
				TNew(i,k) = (TOld(i,k)*(1.0 + + D + F - C - E) - TOld(i+1,k)*(D)
						  + TOld(i+1,k)*(C) - TOld(i,k-1)*F + TOld(i,k+1)*E);
				break;
			case Domain::FAR_FIELD:
				// Apply zero-flux BC (TOld(i+1,k) = TOld(i-1,k))
				TNew(i,k) = (TOld(i,k)*(1.0 + B + D + F - A - C - E) - TOld(i-1,k)*(B + D)
						  + TOld(i-1,k)*(A + C) - TOld(i,k-1)*F + TOld(i,k+1)*E);
				break;
			case Domain::WALL_TOP:
				// Apply zero-flux BC (TOld(i,k+1) = TOld(i,k-1)
				TNew(i,k) = (TOld(i,k)*(1.0 + B + D + F - A - C - E) - TOld(i-1,k)*(B + D)
						  + TOld(i+1,k)*(A + C) - TOld(i,k-1)*F + TOld(i,k-1)*E);
				break;
			case Domain::INTERIOR_AIR:
				TNew(i,k) = foundation.indoorAirTemperature;
				break;
			case Domain::INTERIOR_SLAB:
				{
				double Tair = foundation.indoorAirTemperature;
				double hc = getConvectionCoeff(TOld(i,k),Tair,0.0,1.0,false,0.0);
				double hr = getSimpleInteriorIRCoeff(foundation.wall.interiorEmissivity,
						                             TOld(i,k),Tair);
				double q = 0;

				TNew(i,k) = (domain.getKZM(i,j,k)*TOld(i,k-1)/domain.getDZM(k) +
						  (hc + hr)*Tair + q)/(domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr));
				}
				break;
			case Domain::INTERIOR_WALL:
				{
				double Tair = foundation.indoorAirTemperature;
				double hc = getConvectionCoeff(TOld(i,k),Tair,0.0,1.0,false,PI/2.0);
				double hr = getSimpleInteriorIRCoeff(foundation.wall.interiorEmissivity,
						                             TOld(i,k),Tair);
				double q = 0;

				TNew(i,k) = (domain.getKXP(i,j,k)*TOld(i+1,k)/domain.getDXP(i) +
						  (hc + hr)*Tair + q)/(domain.getKXP(i,j,k)/domain.getDXP(i) + (hc + hr));
				}
				break;
			case Domain::EXTERIOR_AIR:
				TNew(i,k) = getOutdoorTemperature();
				break;
			case Domain::EXTERIOR_GRADE:
				{
				double Tair = getOutdoorTemperature();
				double v = getLocalWindSpeed();
				double eSky = weatherData.skyEmissivity.getValue(getSimTime(tNow));
				double F = getEffectiveExteriorViewFactor(eSky,0.0);
				double hc = getConvectionCoeff(TOld(i,k),Tair,v,1.0,true,0.0);
				double hr = getExteriorIRCoeff(foundation.soilEmissivity,TOld(i,k),Tair,eSky,0.0);
				double q = foundation.soilAbsorptivity*weatherData.globalHorizontalSolar.getValue(getSimTime(tNow));

				TNew(i,k) = (domain.getKZM(i,j,k)*TOld(i,k-1)/domain.getDZM(k) +
						(hc + hr*pow(F,0.25))*Tair + q)/(domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr));
				}
				break;
			case Domain::EXTERIOR_WALL:
				{
				double Tair = getOutdoorTemperature();
				double v = getLocalWindSpeed();
				double eSky = weatherData.skyEmissivity.getValue(getSimTime(tNow));
				double F = getEffectiveExteriorViewFactor(eSky,PI/2.0);
				double hc = getConvectionCoeff(TOld(i,k),Tair,v,1.0,true,PI/2.0);
				double hr = getExteriorIRCoeff(foundation.wall.exteriorEmissivity,TOld(i,k),Tair,eSky,PI/2.0);
				double q = foundation.wall.exteriorAbsorptivity*weatherData.globalHorizontalSolar.getValue(getSimTime(tNow));

				TNew(i,k) = (domain.getKXM(i,j,k)*TOld(i-1,j)/domain.getDXM(i) +
						(hc + hr*pow(F,0.25))*Tair + q)/(domain.getKXM(i,j,k)/domain.getDXM(i) + (hc + hr));
				}
				break;
			case Domain::DEEP_GROUND:
				TNew(i,k) = getDeepGroundTemperature();
				break;
			default:
				TNew(i,k) = (TOld(i,k)*(1.0 + B + D + F - A - C - E) - TOld(i-1,k)*(B + D)
						  + TOld(i+1,k)*(A + C) - TOld(i,k-1)*F + TOld(i,k+1)*E);
				break;
			}
		}
	}

	for (size_t k = 0; j < nZ; ++k)
	{
		for (size_t i = 0; i < nX; ++i)
		{
			// Update old values for next timestep
			TOld(i,k) = TNew(i,k);
		}

	}
}

void Ground::calculateImplicit()
{
	// TODO: remove to incorporate full 3D
	size_t j = 0;

	// TODO: This matrix doesn't need to be constructed each timestep!
	blas::compressed_matrix<double, blas::column_major, 0,
		blas::unbounded_array<int>, blas::unbounded_array<double> > Amat(nX*nZ,nX*nZ);  // Coefficient Matrix
	blas::vector<double> b(nX*nZ), x(nX*nZ);  // constant and unknown vectors

	for (size_t k = 0; k < nZ; k++)
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
				// Apply zero-flux BC (T(i-1,k) = T(i+1,k))
				// Remove A & B at axis to avoid divide by zero
				Amat(i + nX*k, i + nX*k) = (1.0 + C + E - D - F);
				Amat(i + nX*k, (i+1) + nX*k) = D - C;
				Amat(i + nX*k, i + nX*(k-1)) = F;
				Amat(i + nX*k, i + nX*(k+1)) = -E;

				b(i + nX*k) = TOld(i,k);
				break;
			case Domain::FAR_FIELD:
				// Apply zero-flux BC (T(i+1,k) = T(i-1,k))
				Amat(i + nX*k, i + nX*k) = (1.0 + A + C + E - B - D - F);
				Amat(i + nX*k, (i-1) + nX*k) = (B + D - A - C);
				Amat(i + nX*k, i + nX*(k-1)) = F;
				Amat(i + nX*k, i + nX*(k+1)) = -E;

				b(i + nX*k) = TOld(i,k);
				break;
			case Domain::WALL_TOP:
				// Apply zero-flux BC (T(i,k+1) = T(i,k-1)
				Amat(i + nX*k, i + nX*k) = (1.0 + A + C + E - B - D - F);
				Amat(i + nX*k, (i-1) + nX*k) = (B + D);
				Amat(i + nX*k, (i+1) + nX*k) = (-A - C);
				Amat(i + nX*k, i + nX*(k-1)) = F - E;

				b(i + nX*k) = TOld(i,k);
				break;
			case Domain::INTERIOR_AIR:
				Amat(i + nX*k, i + nX*k) = 1.0;

				b(i + nX*k) = foundation.indoorAirTemperature;
				break;
			case Domain::INTERIOR_SLAB:
				{
				double Tair = foundation.indoorAirTemperature;
				double hc = getConvectionCoeff(TOld(i,k),Tair,0.0,1.0,false,0.0);
				double hr = getSimpleInteriorIRCoeff(foundation.wall.interiorEmissivity,
						                             TOld(i,k),Tair);
				double q = 0;

				Amat(i + nX*k, i + nX*k) = domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr);
				Amat(i + nX*k, i + nX*(k-1)) = -domain.getKZM(i,j,k)/domain.getDZM(k);

				b(i + nX*k) = (hc + hr)*Tair + q;
				}
				break;
			case Domain::INTERIOR_WALL:
				{
				double Tair = foundation.indoorAirTemperature;
				double hc = getConvectionCoeff(TOld(i,k),Tair,0.0,1.0,false,PI/2.0);
				double hr = getSimpleInteriorIRCoeff(foundation.wall.interiorEmissivity,
						                             TOld(i,k),Tair);
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
				double hc = getConvectionCoeff(TOld(i,k),Tair,v,1.0,true,0.0);
				double hr = getExteriorIRCoeff(foundation.soilEmissivity,TOld(i,k),Tair,eSky,0.0);
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
				double hc = getConvectionCoeff(TOld(i,k),Tair,v,1.0,true,PI/2.0);
				double hr = getExteriorIRCoeff(foundation.wall.exteriorEmissivity,TOld(i,k),Tair,eSky,PI/2.0);
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

				b(i + nX*k) = TOld(i,k);
				break;
			}
		}
	}

	umf::umf_solve(Amat,x,b);

	for (size_t k = 0; k < nZ; ++k)
	{
		for (size_t i = 0; i < nX; ++i)
		{
			// Read solution into temperature matrix
			TNew(i,k) = x(i + nX*k);

			// Update old values for next timestep
			TOld(i,k) = TNew(i,k);
		}
	}

}

void Ground::calculateCrankNicolson()
{

	// TODO: remove to incorporate full 3D
	size_t j = 0;

	blas::compressed_matrix<double, blas::column_major, 0,
		blas::unbounded_array<int>, blas::unbounded_array<double> > Amat(nX*nZ,nX*nZ);  // Coefficient Matrix
	blas::vector<double> b(nX*nZ), x(nX*nZ);  // constant and unknown vectors

	for (size_t k = 0; k < nZ; k++)
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
				// Apply zero-flux BC (T(i-1,k) = T(i+1,k))
				// Remove A & B at axis to avoid divide by zero
				Amat(i + nX*k, i + nX*k) = (1.0 + C + E - D - F);
				Amat(i + nX*k, (i+1) + nX*k) = D - C;
				Amat(i + nX*k, i + nX*(k-1)) = F;
				Amat(i + nX*k, i + nX*(k+1)) = -E;

				b(i + nX*k) = (TOld(i,k)*(1.0 + D + F - C - E) - TOld(i+1,k)*(D)
						  + TOld(i+1,k)*(C) - TOld(i,k-1)*F + TOld(i,k+1)*E);
				break;
			case Domain::FAR_FIELD:
				// Apply zero-flux BC (T(i+1,k) = T(i-1,k))
				Amat(i + nX*k, i + nX*k) = (1.0 + A + C + E - B - D - F);
				Amat(i + nX*k, (i-1) + nX*k) = (B + D - A - C);
				Amat(i + nX*k, i + nX*(k-1)) = F;
				Amat(i + nX*k, i + nX*(k+1)) = -E;

				b(i + nX*k) = (TOld(i,k)*(1.0 + B + D + F - A - C - E) - TOld(i-1,k)*(B + D)
						  + TOld(i-1,k)*(A + C) - TOld(i,k-1)*F + TOld(i,k+1)*E);
				break;
			case Domain::WALL_TOP:
				// Apply zero-flux BC (T(i,k+1) = T(i,k-1)
				Amat(i + nX*k, i + nX*k) = (1.0 + A + C + E - B - D - F);
				Amat(i + nX*k, (i-1) + nX*k) = (B + D);
				Amat(i + nX*k, (i+1) + nX*k) = (-A - C);
				Amat(i + nX*k, i + nX*(k-1)) = F - E;

				b(i + nX*k) = (TOld(i,k)*(1.0 + B + D + F - A - C - E) - TOld(i-1,k)*(B + D)
						  + TOld(i+1,k)*(A + C) - TOld(i,k-1)*F + TOld(i,k-1)*E);
				break;
			case Domain::INTERIOR_AIR:
				Amat(i + nX*k, i + nX*k) = 1.0;

				b(i + nX*k) = foundation.indoorAirTemperature;
				break;
			case Domain::INTERIOR_WALL:
				{
				double Tair = foundation.indoorAirTemperature;
				double hc = getConvectionCoeff(TOld(i,k),Tair,0.0,1.0,false,PI/2.0);
				double hr = getSimpleInteriorIRCoeff(foundation.wall.interiorEmissivity,
						                             TOld(i,k),Tair);
				double q = 0;

				Amat(i + nX*k, i + nX*k) = domain.getKXP(i,j,k)/domain.getDXP(i) + (hc + hr);
				Amat(i + nX*k, (i+1) + nX*k) = -domain.getKXP(i,j,k)/domain.getDXP(i);

				b(i + nX*k) = (hc + hr)*Tair + q;
				}
				break;
			case Domain::INTERIOR_SLAB:
				{
				double Tair = foundation.indoorAirTemperature;
				double hc = getConvectionCoeff(TOld(i,k),Tair,0.0,1.0,false,0.0);
				double hr = getSimpleInteriorIRCoeff(foundation.wall.interiorEmissivity,
						                             TOld(i,k),Tair);
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
				double hc = getConvectionCoeff(TOld(i,k),Tair,v,1.0,true,0.0);
				double hr = getExteriorIRCoeff(foundation.soilEmissivity,TOld(i,k),Tair,eSky,0.0);
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
				double hc = getConvectionCoeff(TOld(i,k),Tair,v,1.0,true,PI/2.0);
				double hr = getExteriorIRCoeff(foundation.wall.exteriorEmissivity,TOld(i,k),Tair,eSky,PI/2.0);
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

				b(i + nX*k) = (TOld(i,k)*(1.0 + B + D + F - A - C - E) - TOld(i-1,k)*(B + D)
						  + TOld(i+1,k)*(A + C) - TOld(i,k-1)*F + TOld(i,k+1)*E);
				break;
			}
		}
	}

	umf::umf_solve(Amat,x,b);


	for (size_t k = 0; k < nZ; ++k)
	{
		for (size_t i = 0; i < nX; ++i)
		{
			// Read solution into temperature matrix
			TNew(i,k) = x(i + nX*k);

			// Update old values for next timestep
			TOld(i,k) = TNew(i,k);
		}
	}

}

void Ground::calculateSteadyState()
{
	// TODO: remove to incorporate full 3D
	size_t j = 0;

	// TODO: This matrix doesn't need to be constructed each timestep!
	blas::compressed_matrix<double, blas::column_major, 0,
		blas::unbounded_array<int>, blas::unbounded_array<double> > Amat(nX*nZ,nX*nZ);  // Coefficient Matrix
	blas::vector<double> b(nX*nZ), x(nX*nZ);  // constant and unknown vectors

	for (size_t k = 0; k < nZ; k++)
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
				// Apply zero-flux BC (T(i-1,k) = T(i+1,k))
				// Remove A & B at axis to avoid divide by zero
				Amat(i + nX*k, i + nX*k) = (D + F - C - E);
				Amat(i + nX*k, (i+1) + nX*k) = C - D;
				Amat(i + nX*k, i + nX*(k-1)) = -F;
				Amat(i + nX*k, i + nX*(k+1)) = E;

				b(i + nX*k) = 0;
				break;
			case Domain::FAR_FIELD:
				// Apply zero-flux BC (T(i+1,k) = T(i-1,k))
				Amat(i + nX*k, i + nX*k) = (B + D + F - A - C - E);
				Amat(i + nX*k, (i-1) + nX*k) = (-B - D) + (A + C);
				Amat(i + nX*k, i + nX*(k-1)) = -F;
				Amat(i + nX*k, i + nX*(k+1)) = E;

				b(i + nX*k) = 0;
				break;
			case Domain::WALL_TOP:
				// Apply zero-flux BC (T(i,k+1) = T(i,k-1)
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
				double hc = getConvectionCoeff(TOld(i,k),Tair,0.0,1.0,false,0.0);
				double hr = getSimpleInteriorIRCoeff(foundation.wall.interiorEmissivity,
						                             TOld(i,k),Tair);
				double q = 0;

				Amat(i + nX*k, i + nX*k) = domain.getKZM(i,j,k)/domain.getDZM(k) + (hc + hr);
				Amat(i + nX*k, i + nX*(k-1)) = -domain.getKZM(i,j,k)/domain.getDZM(k);

				b(i + nX*k) = (hc + hr)*Tair + q;
				}
				break;
			case Domain::INTERIOR_WALL:
				{
				double Tair = foundation.indoorAirTemperature;
				double hc = getConvectionCoeff(TOld(i,k),Tair,0.0,1.0,false,PI/2.0);
				double hr = getSimpleInteriorIRCoeff(foundation.wall.interiorEmissivity,
						                             TOld(i,k),Tair);
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
				double hc = getConvectionCoeff(TOld(i,k),Tair,v,1.0,true,0.0);
				double hr = getExteriorIRCoeff(foundation.soilEmissivity,TOld(i,k),Tair,eSky,0.0);
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
				double hc = getConvectionCoeff(TOld(i,k),Tair,v,1.0,true,PI/2.0);
				double hr = getExteriorIRCoeff(foundation.wall.exteriorEmissivity,TOld(i,k),Tair,eSky,PI/2.0);
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

	umf::umf_solve(Amat,x,b);


	for (size_t k = 0; k < nZ; ++k)
	{
		for (size_t i = 0; i < nX; ++i)
		{
			// Read solution into temperature matrix
			TNew(i,k) = x(i + nX*k);

			// Update old values for next timestep
			TOld(i,k) = TNew(i,k);
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
	for (size_t i = domain.slabImin; i <= domain.slabImax; i++)
	{
		double Tair = foundation.indoorAirTemperature;
		double h = getConvectionCoeff(TNew(i,domain.slabK),Tair,0.0,1.0,false,0.0);
		double A = 2.0*PI*domain.meshX.deltas[i]*domain.meshX.centers[i]; // for cylindrical coordinates
		//double A = domain.mesher.xDeltas[i];  // for cartesian coordinates
		QSlab[i] = h*A*(Tair - TNew(i,domain.slabK));
	}

	QSlabTotal = accumulate((QSlab).begin(),(QSlab).end(), 0.0)/(PI*pow(foundation.radius,2.0));  // for cylindrical coordinates
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
				double TinF = (TNew(i,k) - 273.15)*9/5 + 32.0;
				TDat.a[i+nX*k] = TinF;
			}

		}

		int nX = rDat.GetNN();
		double rmin = 0.0;
		double rmax = rGrid.a[nX];
		double rrange = rmax - rmin;
		double radius = foundation.radius;
		mglData rTicks(2);
		rTicks.a[0] = radius;
		rTicks.a[1] = rmax;

		string sRadius = str(boost::format("%0.2f") % radius) + " m";
		string sRmax = str(boost::format("%0.2f") % rmax) + " m";
		string rTickString = sRadius + "\n" + sRmax;

		int nZ = zDat.GetNN();
		double zmin = zGrid.a[0];
		double zmax = zGrid.a[nZ];
		double zrange = zmax - zmin;
		mglData zTicks(2);
		zTicks.a[0] = zmin;
		zTicks.a[1] = 0.0;

		string sZmin = str(boost::format("%0.2f") % zmin) + " m";
		string zTickString = sZmin + "\n 0,0\n";

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
			mglPoint bl = mglPoint(foundation.blocks[b].rMin,
			         	 	 	   foundation.blocks[b].zMin,
			         	 	 	   210.0);
			mglPoint br = mglPoint(foundation.blocks[b].rMax,
			         	 	 	   foundation.blocks[b].zMin,
			         	 	 	   210.0);
			mglPoint tr = mglPoint(foundation.blocks[b].rMax,
			         	 	 	   foundation.blocks[b].zMax,
			         	 	 	   210.0);
			mglPoint tl = mglPoint(foundation.blocks[b].rMin,
			         	 	 	   foundation.blocks[b].zMax,
			         	 	 	   210.0);

			gr.Line(bl, br, "k");
			gr.Line(br, tr, "k");
			gr.Line(tr, tl, "k");
			gr.Line(tl, bl, "k");

		}

		ptime tp(simulationControl.startDate,seconds(tNow));
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

ptime Ground::getSimTime(double t)
{
	ptime tp = simulationControl.startTime + seconds(t);

	return tp;
}

double Ground::getInitialTemperature(double t, double z)
{
	if (foundation.initializationMethod == Foundation::IM_KUSUDA)
	{
		double minDryBulb = weatherData.dryBulbTemp.getMin();
		double maxDryBulb = weatherData.dryBulbTemp.getMax();

		ptime tp = getSimTime(t);
		greg_year year = tp.date().year();
		date dayBegin(year,Jan,1);
		ptime tYearStart(dayBegin);
		time_duration tSinceYearStart = tp - tYearStart;
		double trel = tSinceYearStart.total_seconds();
		double tshift = 0;  // Assume min temperature occurs at the beginning of the year
		double seconds_in_day = 60.0*60.0*24.0;
		double Tamp = (maxDryBulb - minDryBulb) / 2.0;
		double diff = foundation.soil.k/(foundation.soil.rho*foundation.soil.cp);

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

