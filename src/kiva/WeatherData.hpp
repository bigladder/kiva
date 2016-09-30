/* Copyright (c) 2012-2016 Big Ladder Software. All rights reserved.
* See the LICENSE file for additional terms and conditions. */

#ifndef WEATHERDATA_H_
#define WEATHERDATA_H_

#include <fstream>
#include <cmath>
#include <numeric>
#include <iostream>

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>

std::istream& safeGetline(std::istream& is, std::string& t);


class HourlyData: public std::vector<double>
{
public:

  enum DataType
  {
    DT_SOLAR,
    DT_METEOROLOGICAL
  };

  DataType dataType;

  double getValue(boost::posix_time::ptime t);
  double getAverage();
  double getMin();
  double getMax();

  HourlyData();
  HourlyData(DataType dT);

};

class WeatherData
{
public:

  std::string city;
  std::string state;
  std::string country;
  double timezone;
  double latitude;
  double longitude;
  double elevation;

  double hourOfMinimumTemperature;
  std::vector<double> dailyAverageTemperatures;
  std::vector<double> monthlyAverageTemperatures;
  double minimumAverageMontlyTemperature;
  double maximumAverageMontlyTemperature;

  HourlyData dryBulbTemp;
  //HourlyData wetBulbTemp;
  HourlyData dewPointTemp;
  HourlyData atmosphericPressure;
  //HourlyData density;
  HourlyData relativeHumidity;
  //HourlyData humidityRatio;
  //HourlyData enthalpy;

  //HourlyData solarAngle;
  //HourlyData globalHorizontalSolar;
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
  HourlyData altitude;
  HourlyData azimuth;

public:
  WeatherData(std::string weatherFile);

private:
  void importEPW(std::string epwFile);

};

#endif /* WEATHERDATA_H_ */
