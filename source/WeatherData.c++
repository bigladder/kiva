/* WeatherData.c++ is part of Kiva (Written by Neal Kruis)
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

#ifndef WeatherData_CPP
#define WeatherData_CPP

#include "WeatherData.h"

static const double PI = 4.0*atan(1.0);

HourlyData::HourlyData()
{
  dataType = DT_METEOROLOGICAL;
}

HourlyData::HourlyData(DataType dT)
{
  dataType = dT;
}

double HourlyData::getValue(boost::posix_time::ptime t)
{
  long year = t.date().year();
  boost::gregorian::date newYearDate(year,boost::gregorian::Jan,1);
  boost::posix_time::ptime newYear(newYearDate);
  boost::posix_time::time_duration duration = t - newYear;
  long second = duration.total_seconds();
  double hour = second/60.0/60.0;
  long hour_before = hour;
  long hour_after = hour_before + 1;

  if (boost::gregorian::gregorian_calendar::is_leap_year(year) && hour >= 1416.0)
  {
      hour_before -= 24;
      hour_after -= 24;
      hour -= 24.0;
  }

  if (hour_after == 8760) hour_after = 0;
  if (hour_before == -1) hour_before = 8759;

  double value;
  if (dataType == DT_METEOROLOGICAL)
  {
    // Meteorological measurements are made at the hour indicated (Hour 1 = 1:00 AM)

    double frac = hour - hour_before;
    // Decrement hours since the zeroeth element corresponds to the value at the
    // first hour
    hour_before -= 1;
    hour_after -= 1;

    if (hour_before == -1) hour_before = 8759;
    if (hour_after == -1) hour_after = 8759;

    value = (*this)[hour_before]*(1.0 - frac) + (*this)[hour_after]*frac;
  }
  else // if (dataType == DT_SOLAR)
  {
    // Solar radiation values represent the energy received during the 60 minutes
    // preceding the hour indicated (Hour 1 = 12:00 AM - 1:00AM)
    value = (*this)[hour_before];
  }
  return value;

}

double HourlyData::getAverage()
{
  return accumulate((*this).begin(),(*this).end(), 0.0) / (*this).size();
}

double HourlyData::getMin()
{
  return *min_element((*this).begin(),(*this).end());
}

double HourlyData::getMax()
{
  return *max_element((*this).begin(),(*this).end());
}

WeatherData::WeatherData(std::string weatherFile)
{
  globalHorizontalSolar.dataType = HourlyData::DT_SOLAR;
  directNormalSolar.dataType = HourlyData::DT_SOLAR;
  diffuseHorizontalSolar.dataType = HourlyData::DT_SOLAR;
  importEPW(weatherFile);
}

void WeatherData::importEPW(std::string epwFile)
{

  std::ifstream inf;
  inf.open(epwFile.c_str());

    if (!inf)
    {
        // Print an error and exit
        std::cerr << "Unable to read EPW file" << std::endl;
        exit(1);
    }

    // While there's still stuff left to read
    std::string line;
    std::vector <std::string> columns;

    typedef boost::tokenizer< boost::escaped_list_separator<char> > Tokenizer;

    int row = 0;
  long hour = 0;

    while (getline(inf,line))
    {
      row += 1;
      Tokenizer tok(line, boost::escaped_list_separator<char>("\\",",","\""));

      for (boost::tokenizer< boost::escaped_list_separator<char> >::iterator
          i(tok.begin()); i!=tok.end(); ++i)
      {
        columns.push_back(*i);
      }

      if (row == 1)
      {
        // Site information
        city = columns[1];
        state = columns[2];
        country = columns[3];
        latitude = double(boost::lexical_cast<double>(columns[6]));
        longitude = double(boost::lexical_cast<double>(columns[7]));
        timezone = double(boost::lexical_cast<double>(columns[8]));
        //elevation = double(boost::lexical_cast<double>(columns[9]));
      }

      boost::gregorian::date newYear(2013,boost::gregorian::Jan,1);
      boost::posix_time::ptime newYearDate(newYear);

      if (row > 8)
      {
          boost::posix_time::ptime dateTime = newYearDate + boost::posix_time::hours(hour);
          double hourOfDay = dateTime.time_of_day().total_seconds()/3600.0;

        double Tdb = double(boost::lexical_cast<double>(columns[6])) +
                  273.15;

        dryBulbTemp.push_back(Tdb);  // [K]

        double Tdp = double(boost::lexical_cast<double>(columns[7])) +
                  273.15;

        dewPointTemp.push_back(Tdp);  // [K]

        relativeHumidity.push_back(
            double(boost::lexical_cast<double>(columns[8]))/100.0);  // [frac]

        atmosphericPressure.push_back(
            double(boost::lexical_cast<double>(columns[9])));  // [Pa]

        // Note: global horizontal solar can be calculated using solar position,
        // direct normal solar and diffuse solar.
        // double qGH = double(boost::lexical_cast<double>(columns[13]));

        double qDN = double(boost::lexical_cast<double>(columns[14]));

        directNormalSolar.push_back(qDN);  // [W/m2]

        double qDH = double(boost::lexical_cast<double>(columns[15]));

        diffuseHorizontalSolar.push_back(qDH);  // [W/m2]

        windDirection.push_back(double(boost::lexical_cast<double>(columns[20]))*
              PI/180.0);  // [rad]

        windSpeed.push_back(
            double(boost::lexical_cast<double>(columns[21])));  // [m/s]

        double fc = double(boost::lexical_cast<double>(columns[22]));

        totalSkyCover.push_back(fc);  // [tenths]

        opaqueSkyCover.push_back(
            double(boost::lexical_cast<double>(columns[23])));  // [tenths]

        double eSky = 0.787 + 0.764*log(Tdp/Tdb)*
            (1 + 0.0224*fc + 0.0035*pow(fc,2) + 0.00028*pow(fc,3));

        skyEmissivity.push_back(eSky);  // [frac]

        double Tsky = pow(eSky,0.25)*Tdb;

        skyTemp.push_back(Tsky);

        // Solar calculations (from Duffie & Beckman)
        double day = dateTime.date().day_of_year();
        double B = (day - 1)*360/365*PI/180;
        double EOT = 229.2*(0.000075 + 0.001868*cos(B) - 0.032077*sin(B)
                 - 0.014615*cos(2*B) - 0.04089*sin(2*B));
        double declination = 0.006918 - 0.399912*cos(B) + 0.070257*sin(B)
                             - 0.006758*cos(2*B) + 0.000907*sin(2*B)
                             - 0.002697*cos(3*B) + 0.00148*sin(3*B);
        double solarHour = (hourOfDay - 0.5) + (4*(timezone*15.0 - longitude) + EOT)/60.0;
        double hourAngle = (solarHour - 12.00)*15.0;

        double sinLat = sin(latitude*PI/180);
        double cosLat = cos(latitude*PI/180);
        double sinDec = sin(declination*PI/180);
        double cosDec = cos(declination*PI/180);
        double cosHA = cos(hourAngle*PI/180);

        double alt = asin(sinLat*sinDec + cosLat*cosDec*cosHA);
        altitude.push_back(alt);

        double qGH = cos(PI/2 - alt)*qDN + qDH;
        globalHorizontalSolar.push_back(qGH);  // [W/m2]

        double azi = PI + acos((sin(alt)*sinLat - sinDec)/(cos(alt)*cosLat));

        if (hourAngle < 0.0)
        {
          azi = 2*PI - azi;
        }

        azimuth.push_back(azi);

        hour++;
      }
      columns.clear();

    }

    inf.close();
}

#endif
