/*
 * ADESolver.hpp
 *
 *  Created on: Oct 18, 2012
 *      Author: nkruis
 */

#ifndef GROUND_HPP_
#define GROUND_HPP_

#include "Mesher.h"
#include "Domain.h"
#include "WeatherData.h"
#include "Input.h"
#include "Algorithms.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <cmath>
#include <mgl2/mgl.h>
#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_sparse.hpp>
#include <boost/numeric/bindings/umfpack/umfpack.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

static const double PI = 4.0*atan(1.0);

using namespace std;
using namespace boost::posix_time;
namespace blas = boost::numeric::ublas;
namespace umf = boost::numeric::bindings::umfpack;

class Ground
{
private:

	double timestep;
	Foundation &foundation;
	SimulationControl &simulationControl;
	WeatherData &weatherData;

	double tNow;
	size_t nR, nZ;
	double annualAverageDryBulbTemperature;

	// Data structures
	blas::matrix<double> U; // ADE upper sweep, n+1
	blas::matrix<double> UOld; // ADE upper sweep, n
	blas::matrix<double> V; // ADE lower sweep, n+1
	blas::matrix<double> VOld; // ADE lower sweep, n
	blas::matrix<double> TNew; // solution, n+1
	blas::matrix<double> TOld; // solution, n

	blas::vector<double> QSlab; // Heat Fluxes through the slab


	// Plotting variables
	mglData TDat, rDat, zDat, cDat, rGrid, zGrid, TGrid;
	mglGraph gr;		// class for plot drawing

	double nextPlotTime, plotFreq;
	string startString;
	bool makePlot;
	Domain domain;

public:

	double QSlabTotal;

	// Constructor
	Ground(WeatherData &weatherData, Foundation &foundation, SimulationControl &simulationControl);

	// Destructor
	virtual ~Ground();

	// Calculator
	void calculate(double t);

private:

	// Initializers (Called from constructor)
	void buildDomain();

	void initializeConditions();

	void initializePlot();

	// Calculators (Called from main calculator)
	void calculateADE();

	void calculateADEUpwardSweep();

	void calculateADEDownwardSweep();

	void calculateExplicit();

	void calculateImplicit();

	void calculateCrankNicolson();

	void calculateSteadyState();

	void plot();

	// Misc. Functions
	ptime getSimTime(double t);

	double getInitialTemperature(double r, double z);

	double getDeepGroundTemperature();

	double getConvectionCoeff(double Tsurf,
							  double Tamb,
							  double Vair,
							  double roughness,
							  bool isExterior);

	double getOutdoorTemperature();
	double getLocalWindSpeed();
};



#endif /* GROUND_HPP_ */
