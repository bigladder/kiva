#ifndef ConvectionAlgorithms_HPP
#define ConvectionAlgorithms_HPP

#include <cmath>

using namespace std;

// TODO: use defaulting inputs
double getDOE2ConvectionCoeff(double tilt,
		  	  	  	  	  	  double azimuth,
		  	  	  	  	  	  double windDirection,
		  	  	  	  	  	  double Tsurf,
		  	  	  	  	  	  double Tamb,
		  	  	  	  	  	  double Vair,
		  	  	  	  	  	  double roughness);

bool isWindward(double tilt, double azimuth, double windDirection);

double getIRCoeff(double eSurf,
		          double Tsurf,
		          double Tamb,
		          double eSky,
		          double tilt);

#endif
