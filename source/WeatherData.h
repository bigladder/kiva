/* WeatherData.h is part of Kiva (Written by Neal Kruis)
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

#ifndef WEATHERDATA_H_
#define WEATHERDATA_H_

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <numeric>

using namespace std;
using namespace boost::posix_time;
using namespace boost::gregorian;

class HourlyData: public vector<double>
{
public:

	double getValue(ptime t);
	double getAverage();
	double getMin();
	double getMax();
};

class WeatherData
{
public:

	string city;
	string state;
	string country;

	HourlyData dryBulbTemp;
	//HourlyData wetBulbTemp;
	HourlyData dewPointTemp;
	HourlyData atmosphericPressure;
	//HourlyData density;
	HourlyData relativeHumidity;
	//HourlyData humidityRatio;
	//HourlyData enthalpy;

	//HourlyData solarAngle;
	HourlyData globalHorizontalSolar;
	HourlyData directNormalSolar;
	HourlyData diffuseHorizontalSolar;

	HourlyData windDirection;
	HourlyData windSpeed;

	HourlyData totalSkyCover;
	HourlyData opaqueSkyCover;
	HourlyData skyEmissivity;
	HourlyData skyTemp;
	//HourlyData cloudType;
	//HourlyData snowFlag;
	//HourlyData rainFlag;

public:
	WeatherData(string weatherFile);

private:
	void importEPW(string epwFile);

};

#endif /* WEATHERDATA_H_ */

