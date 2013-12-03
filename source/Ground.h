/* Ground.h is part of Kiva (Written by Neal Kruis)
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

#ifndef GROUND_HPP_
#define GROUND_HPP_

#include "Mesher.h"
#include "Domain.h"
#include "WeatherData.h"
#include "Input.h"
#include "Algorithms.h"
#include "GroundPlot.h"

#include <cmath>
#include <vector>

#include <mgl2/mgl.h>

#include <boost/thread.hpp>
#include <boost/bind.hpp>

#include <boost/multi_array.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_sparse.hpp>
#include <boost/numeric/bindings/umfpack/umfpack.hpp>

namespace umf = boost::numeric::bindings::umfpack;

class Ground
{
private:

	double timestep;
	Foundation &foundation;
	SimulationControl &simulationControl;
	WeatherData &weatherData;

	double tNow;
	size_t nX, nY, nZ;
	double annualAverageDryBulbTemperature;

	Domain domain;

	// Data structures
	boost::multi_array<double, 3> U; // ADE upper sweep, n+1
	boost::multi_array<double, 3> UOld; // ADE upper sweep, n
	boost::multi_array<double, 3> V; // ADE lower sweep, n+1
	boost::multi_array<double, 3> VOld; // ADE lower sweep, n

	boost::numeric::ublas::compressed_matrix<double,
			boost::numeric::ublas::column_major, 0,
			boost::numeric::ublas::unbounded_array<int>,
			boost::numeric::ublas::unbounded_array<double> > Amat;

	boost::numeric::ublas::vector<double> b, x;

	boost::multi_array<double, 3> TOld; // solution, n

	std::vector<GroundPlot> plots;

public:

	boost::multi_array<double, 3> TNew; // solution, n+1

	// Constructor
	Ground(WeatherData &weatherData, Foundation &foundation, SimulationControl &simulationControl);

	// Destructor
	virtual ~Ground();

	// Calculator
	void calculate(double t);

	double getSurfaceAverageHeatFlux(std::string surfaceName);

	double getBelowSlabTemperature(std::string surfaceName);

private:

	// Initializers (Called from constructor)
	void buildDomain();

	void initializeConditions();

	void initializePlots();

	// Calculators (Called from main calculator)
	void calculateADE();

	void calculateADEUpwardSweep();

	void calculateADEDownwardSweep();

	void calculateExplicit();

	void calculateMatrix(Foundation::NumericalScheme scheme);

	void calculateADI(int dim);

	void plot();

	// Misc. Functions
	boost::posix_time::ptime getSimTime(double t);

	double getInitialTemperature(double t, double z);

	double getDeepGroundTemperature();

	double getConvectionCoeff(double Tsurf,
							  double Tamb,
							  double Vair,
							  double roughness,
							  bool isExterior,
							  double tilt);

	double getOutdoorTemperature();

	double getLocalWindSpeed();
};



#endif /* GROUND_HPP_ */
