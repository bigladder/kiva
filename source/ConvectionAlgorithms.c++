#ifndef ConvectionAlgorithms_CPP
#define ConvectionAlgorithms_CPP

#include "ConvectionAlgorithms.h"

static const double PI = 4.0*atan(1.0);

double getDOE2ConvectionCoeff(double tilt,
						  	  double azimuth,
						  	  double windDirection,
						  	  double Tsurf,
						  	  double Tamb,
						  	  double Vair,
						  	  double roughness)
{
	/* Based on the DOE-2 convection model as used in EnergyPlus
	 *
	 * Roughness factors:
	 * Very Rough = 2.17
	 * Rough = 1.67
	 * Medium Rough = 1.52
	 * Medium Smooth = 1.13
	 * Smooth = 1.11
	 * Very Smooth = 1.0
	*/

	double hn;

	if (cos(tilt) == 0.0)
	{
		hn = 1.31*(pow(fabs(Tsurf - Tamb),1.0/3.0));
	}
	else if ((cos(tilt) < 0.0 && Tsurf < Tamb) ||
			 (cos(tilt) > 0.0 && Tsurf > Tamb))
	{
		hn = 9.482*(pow(fabs(Tsurf - Tamb),1.0/3.0))/
				(7.238 - fabs(cos(tilt)));
	}
	else if ((cos(tilt) < 0.0 && Tsurf > Tamb) ||
			 (cos(tilt) > 0.0 && Tsurf < Tamb))
	{
		hn = 1.810*(pow(fabs(Tsurf - Tamb),1.0/3.0))/
				(1.382 - fabs(cos(tilt)));
	}

	double hcGlass;

	if (isWindward(tilt, azimuth, windDirection))
	{
		hcGlass = sqrt(pow(hn,2) + pow(pow(3.26*Vair,0.89),2));
	}
	else
	{
		hcGlass = sqrt(pow(hn,2) + pow(pow(3.55*Vair,0.617),2));
	}

	double hf = roughness*(hcGlass - hn);
	return hn + hf;
}

bool isWindward(double tilt, double azimuth, double windDirection)
{
	if (fabs(cos(tilt) < 0.98))
	{
		double diff = fabs(windDirection - azimuth);
		if ((diff - PI) > 0.001)
		{
			diff -= 2*PI;
		}
		if (fabs((diff) - 100.0*PI/180.0) > 0.001)
		{
			return false;
		}
		else return true;
	}
	else return true;
}

#endif
