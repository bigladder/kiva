/*
 * WeatherData.cpp
 *
 *  Created on: Oct 25, 2012
 *      Author: nkruis
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
	//HourlyData cloudType;
	//HourlyData snowFlag;
	//HourlyData rainFlag;

public:
	WeatherData(string weatherFile);

private:
	void importEPW(string epwFile);

};

#endif /* WEATHERDATA_H_ */

