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
	domain.setDomain(foundation);

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

	if (foundation.numericalScheme == Foundation::NS_CRANK_NICOLSON ||
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
					TOld[i][j][k]= getInitialTemperature(tNow,
							domain.meshZ.centers[k]);
				}
			}
		}
	}
}

void Ground::initializePlot()
{
	startString = boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time());

	size_t contourLevels = 13;

	mglData TRef(nX, nY),
			xRef(nX),
			yRef(nY),
			xGridRef(nX + 1),
			yGridRef(nY + 1),
			TGridRef(nX + 1, nY + 1),
			cRef(contourLevels);


	TDat = TRef;
	xDat = xRef;
	yDat = yRef;
	xGrid = xGridRef;
	yGrid = yGridRef;
	TGrid = TGridRef;
	cDat = cRef;

	nextPlotTime = 0.0;
	plotFreq = foundation.outputAnimation.frequency.total_seconds();

	double aspect = 1.0;
	int height = foundation.outputAnimation.size;
	int width = height*aspect;
	gr.SetSize(width,height);

	gr.StartGIF((foundation.outputAnimation.name + ".gif").c_str(),100);

	xGrid.a[0] = domain.meshX.dividers[0];

	for(size_t i = 0; i < nX; i++)
	{
		xDat.a[i] = domain.meshX.centers[i];
		xGrid.a[i + 1] = domain.meshX.dividers[i + 1];
	}

	yGrid.a[0] = domain.meshY.dividers[0];

	for(size_t j = 0; j < nY; j++)
	{
		yDat.a[j] = domain.meshY.centers[j];
		yGrid.a[j + 1] = domain.meshY.dividers[j + 1];
	}

	for(size_t j = 0; j <= nY; j++)
	{
		for(size_t i = 0; i <= nX; i++)
		{
			TGrid.a[i+nX*j] = 200.0;
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

					double A = 0;
					double B = 0;
					if (i != 0)
					{
						double r = domain.meshX.centers[i];
						A = domain.cell[i][j][k].cxp_c*theta/r;
						B = domain.cell[i][j][k].cxm_c*theta/r;
					}
					double C = domain.cell[i][j][k].cxp*theta;
					double D = domain.cell[i][j][k].cxm*theta;
					double E = domain.cell[i][j][k].czp*theta;
					double F = domain.cell[i][j][k].czm*theta;
					double G = domain.cell[i][j][k].cyp*theta;
					double H = domain.cell[i][j][k].cym*theta;

					if (foundation.coordinateSystem == Foundation::CS_3D)
						U[i][j][k] = (UOld[i][j][k]*(1.0 - A - C - E - G)
								- U[i-1][j][k]*(B + D)
								+ UOld[i+1][j][k]*(A + C)
								- U[i][j][k-1]*F
								+ UOld[i][j][k+1]*E
								- U[i][j-1][k]*H
								+ UOld[i][j+1][k]*G) /
								(1.0 - B - D - F - H);
					else
						U[i][j][k] = (UOld[i][j][k]*(1.0 - A - C - E)
								- U[i-1][j][k]*(B + D)
								+ UOld[i+1][j][k]*(A + C)
								- U[i][j][k-1]*F
								+ UOld[i][j][k+1]*E) /
								(1.0 - B - D - F);

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

					double A = 0;
					double B = 0;
					if (i != 0)
					{
						double r = domain.meshX.centers[i];
						A = domain.cell[i][j][k].cxp_c*theta/r;
						B = domain.cell[i][j][k].cxm_c*theta/r;
					}
					double C = domain.cell[i][j][k].cxp*theta;
					double D = domain.cell[i][j][k].cxm*theta;
					double E = domain.cell[i][j][k].czp*theta;
					double F = domain.cell[i][j][k].czm*theta;
					double G = domain.cell[i][j][k].cyp*theta;
					double H = domain.cell[i][j][k].cym*theta;

					if (foundation.coordinateSystem == Foundation::CS_3D)
						V[i][j][k] = (VOld[i][j][k]*(1.0 + B + D + F + H)
								- VOld[i-1][j][k]*(B + D)
								+ V[i+1][j][k]*(A + C)
								- VOld[i][j][k-1]*F
								+ V[i][j][k+1]*E
								- VOld[i][j-1][k]*H
								+ V[i][j+1][k]*G) /
								(1.0 + A + C + E + G);
					else
						V[i][j][k] = (VOld[i][j][k]*(1.0 + B + D + F)
								- VOld[i-1][j][k]*(B + D)
								+ V[i+1][j][k]*(A + C)
								- VOld[i][j][k-1]*F
								+ V[i][j][k+1]*E) /
								(1.0 + A + C + E);
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

					double A = 0;
					double B = 0;
					if (i != 0)
					{
						double r = domain.meshX.centers[i];
						A = domain.cell[i][j][k].cxp_c*theta/r;
						B = domain.cell[i][j][k].cxm_c*theta/r;
					}
					double C = domain.cell[i][j][k].cxp*theta;
					double D = domain.cell[i][j][k].cxm*theta;
					double E = domain.cell[i][j][k].czp*theta;
					double F = domain.cell[i][j][k].czm*theta;
					double G = domain.cell[i][j][k].cyp*theta;
					double H = domain.cell[i][j][k].cym*theta;

					if (foundation.coordinateSystem == Foundation::CS_3D)
						TNew[i][j][k] = TOld[i][j][k]*(1.0 + B + D + F + H - A - C - E - G)
								- TOld[i-1][j][k]*(B + D)
								+ TOld[i+1][j][k]*(A + C)
								- TOld[i][j][k-1]*F
								+ TOld[i][j][k+1]*E
								- TOld[i][j-1][k]*H
								+ TOld[i][j+1][k]*G;
					else
						TNew[i][j][k] = TOld[i][j][k]*(1.0 + B + D + F - A - C - E)
								- TOld[i-1][j][k]*(B + D)
								+ TOld[i+1][j][k]*(A + C)
								- TOld[i][j][k-1]*F
								+ TOld[i][j][k+1]*E;
					}
					break;
				}
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
						double A = 0;
						double B = 0;
						if (i != 0)
						{
							double r = domain.meshX.centers[i];
							A = domain.cell[i][j][k].cxp_c/r;
							B = domain.cell[i][j][k].cxm_c/r;
						}
						double C = domain.cell[i][j][k].cxp;
						double D = domain.cell[i][j][k].cxm;
						double E = domain.cell[i][j][k].czp;
						double F = domain.cell[i][j][k].czm;
						double G = domain.cell[i][j][k].cyp;
						double H = domain.cell[i][j][k].cym;

						if (foundation.coordinateSystem == Foundation::CS_3D)
						{
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*k) = (B + D + F + H - A - C - E - G);
							Amat(i + nX*j + nX*nY*k, (i-1) + nX*j + nX*nY*k) = (-B - D);
							Amat(i + nX*j + nX*nY*k, (i+1) + nX*j + nX*nY*k) = (A + C);
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*(k-1)) = -F;
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*(k+1)) = E;
							Amat(i + nX*j + nX*nY*k, i + nX*(j-1) + nX*nY*k) = -H;
							Amat(i + nX*j + nX*nY*k, i + nX*(j+1) + nX*nY*k) = G;

							b(i + nX*j + nX*nY*k) = 0;
						}
						else
						{
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*k) = (B + D + F - A - C - E);
							Amat(i + nX*j + nX*nY*k, (i-1) + nX*j + nX*nY*k) = (-B - D);
							Amat(i + nX*j + nX*nY*k, (i+1) + nX*j + nX*nY*k) = (A + C);
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*(k-1)) = -F;
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*(k+1)) = E;

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

						double A = 0;
						double B = 0;
						if (i != 0)
						{
							double r = domain.meshX.centers[i];
							A = domain.cell[i][j][k].cxp_c*theta/r;
							B = domain.cell[i][j][k].cxm_c*theta/r;
						}
						double C = domain.cell[i][j][k].cxp*theta;
						double D = domain.cell[i][j][k].cxm*theta;
						double E = domain.cell[i][j][k].czp*theta;
						double F = domain.cell[i][j][k].czm*theta;
						double G = domain.cell[i][j][k].cyp*theta;
						double H = domain.cell[i][j][k].cym*theta;

						if (foundation.coordinateSystem == Foundation::CS_3D)
						{
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*k) = (1.0 + f*(A + C + E + G - B - D - F - H));
							Amat(i + nX*j + nX*nY*k, (i-1) + nX*j + nX*nY*k) = f*(B + D);
							Amat(i + nX*j + nX*nY*k, (i+1) + nX*j + nX*nY*k) = f*(-A - C);
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*(k-1)) = f*F;
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*(k+1)) = f*(-E);
							Amat(i + nX*j + nX*nY*k, i + nX*(j-1) + nX*nY*k) = f*H;
							Amat(i + nX*j + nX*nY*k, i + nX*(j+1) + nX*nY*k) = f*(-G);

							b(i + nX*j + nX*nY*k)
									= TOld[i][j][k]*(1.0 + (1-f)*(B + D + F + H - A - C - E - G))
									- TOld[i-1][j][k]*(1-f)*(B + D)
									+ TOld[i+1][j][k]*(1-f)*(A + C)
									- TOld[i][j][k-1]*(1-f)*F
									+ TOld[i][j][k+1]*(1-f)*E
									- TOld[i][j-1][k]*(1-f)*H
									+ TOld[i][j+1][k]*(1-f)*G;
						}
						else
						{
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*k) = (1.0 + f*(A + C + E - B - D - F));
							Amat(i + nX*j + nX*nY*k, (i-1) + nX*j + nX*nY*k) = f*(B + D);
							Amat(i + nX*j + nX*nY*k, (i+1) + nX*j + nX*nY*k) = f*(-A - C);
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*(k-1)) = f*F;
							Amat(i + nX*j + nX*nY*k, i + nX*j + nX*nY*(k+1)) = f*(-E);

							b(i + nX*j + nX*nY*k)
									= TOld[i][j][k]*(1.0 + (1-f)*(B + D + F - A - C - E))
									- TOld[i-1][j][k]*(1-f)*(B + D)
									+ TOld[i+1][j][k]*(1-f)*(A + C)
									- TOld[i][j][k-1]*(1-f)*F
									+ TOld[i][j][k+1]*(1-f)*E;
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
			}
		}
	}
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
	// Calculate Heat Fluxes
	QSlabTotal = getSurfaceAverageHeatFlux("Slab Interior");

	if (makePlot)
		plot();
}

void Ground::plot()
{

	if (tNow >= nextPlotTime)
	{
		for(size_t j = 0; j < nY; j++)
		{
			for(size_t i = 0; i < nX; i++)
			{
				size_t k = nZ/2;
			    double TinF = (TNew[i][j][k] - 273.15)*9/5 + 32.0;
				TDat.a[i+nX*j] = TinF;
			}

		}

		//int nX = xDat.GetNN();
		//double xmin = 0.0;
		//double xmax = xGrid.a[nX];
		//double xrange = xmax - xmin;
		//double length = foundation.effectiveLength;
		//mglData xTicks(2);
		//xTicks.a[0] = length;
		//xTicks.a[1] = xmax;

		//std::string sLength = str(boost::format("%0.2f") % length) + " m";
		//std::string sRmax = str(boost::format("%0.2f") % xmax) + " m";
		//std::string rTickString = sLength + "\n" + sRmax;

		//int nY = yDat.GetNN();
		//double ymin = yGrid.a[0];
		//double ymax = yGrid.a[nY];
		//double yrange = ymax - ymin;
		//mglData yTicks(2);
		//yTicks.a[0] = ymin;
		//yTicks.a[1] = 0.0;

		//std::string sYmin = str(boost::format("%0.2f") % ymin) + " m";
		//std::string yTickString = sYmin + "\n 0,0\n";

		int nT = cDat.GetNN();
		double Tmin = cDat.a[0];
		double Tmax = cDat.a[nT - 1];
		double Tstep = cDat.a[1] - cDat.a[0];

		gr.NewFrame();

		// Plot
		//gr.MultiPlot(3,1,0,2,1,"_");
		gr.LoadFont("none");
		gr.SetFontSize(2.0);
		//gr.SetOrigin(0.0, ymax);
		gr.SetRange('x', xGrid);
		gr.SetRange('y', yGrid);
		gr.SetRange('c', Tmin, Tmax);
		gr.SetRange('z', Tmin, Tmax);
		gr.SetTicks('c', Tstep, nT, Tmin);
		//gr.SetTicksVal('x', rTicks, rTickString.c_str());
		//gr.SetTicksVal('y', zTicks, zTickString.c_str());
		//gr.SetTickLen(-0.0001);
		//gr.Aspect(xrange, yrange);
		gr.Axis("yU");
		gr.Axis("x");
		gr.Colorbar("_");
		//gr.Puts(mglPoint(12.5,-10), "Temperature \\textdegree C", ":C");
		gr.Box("k",false);
		gr.Dens(xDat, yDat, TDat);
		if (foundation.outputAnimation.contours)
			gr.Cont(cDat, xDat, yDat, TDat,"H");
		if (foundation.outputAnimation.gradients)
			gr.Grad(xDat, yDat, TDat);
		if (foundation.outputAnimation.grid)
			gr.Grid(xGrid, yGrid, TGrid, "W");

		// Draw blocks
		if (false)
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
		gr.Puts(mglPoint(-0.10,3.0), to_simple_string(tp).c_str(), ":C");

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

