/* Copyright (c) 2012-2017 Big Ladder Software LLC. All rights reserved.
* See the LICENSE file for additional terms and conditions. */

#ifndef WeatherData_CPP
#define WeatherData_CPP

#include "WeatherData.hpp"

static const double PI = 4.0*atan(1.0);

std::istream& safeGetline(std::istream& is, std::string& t)
{
    t.clear();

    // The characters in the stream are read one-by-one using a std::streambuf.
    // That is faster than reading them one-by-one using the std::istream.
    // Code that uses streambuf this way must be guarded by a sentry object.
    // The sentry object performs various tasks,
    // such as thread synchronization and updating the stream state.

    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();

    for(;;) {
        int c = sb->sbumpc();
        switch (c) {
        case '\n':
            return is;
        case '\r':
            if(sb->sgetc() == '\n')
                sb->sbumpc();
            return is;
        case EOF:
            // Also handle the case when the last line has no line ending
            if(t.empty())
                is.setstate(std::ios::eofbit);
            return is;
        default:
            t += (char)c;
        }
    }
}

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

  // Handle leap year days
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
  //globalHorizontalSolar.dataType = HourlyData::DT_SOLAR;
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
      exit(EXIT_FAILURE);
  }

  // While there's still stuff left to read
  std::string line;
  std::vector <std::string> columns;

  typedef boost::tokenizer< boost::escaped_list_separator<char> > Tokenizer;

  int row = 0;
  long hour = 0;

  hourOfMinimumTemperature = hour;
  double Tmin = 9999;

  std::vector<int> months;
  std::vector<int> days;

  int dayOfYear = 1;
  int previousDayOfMonth = 1;

  while (!safeGetline(inf,line).eof())
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

      double Tdb = double(boost::lexical_cast<double>(columns[6])) + 273.15;

      months.push_back(int(boost::lexical_cast<int>(columns[1])));

      int dayOfMonth = int(boost::lexical_cast<int>(columns[2]));

      if (dayOfMonth != previousDayOfMonth)
        dayOfYear += 1;

      days.push_back(dayOfYear);

      previousDayOfMonth = dayOfMonth;

      if (Tdb < Tmin)
      {
        Tmin = Tdb;
        hourOfMinimumTemperature = hour;
      }

      dryBulbTemp.push_back(Tdb);  // [K]

      double Tdp = double(boost::lexical_cast<double>(columns[7])) +
                273.15;

      dewPointTemp.push_back(Tdp);  // [K]

      relativeHumidity.push_back(
          double(boost::lexical_cast<double>(columns[8]))/100.0);  // [frac]

      atmosphericPressure.push_back(
          double(boost::lexical_cast<double>(columns[9])));  // [Pa]

      double qDN = double(boost::lexical_cast<double>(columns[14]));

      directNormalSolar.push_back(qDN);  // [W/m2]

      double qDH = double(boost::lexical_cast<double>(columns[15]));

      diffuseHorizontalSolar.push_back(qDH);  // [W/m2]

      windDirection.push_back(double(boost::lexical_cast<double>(columns[20]))*PI/180.0);  // [rad]

      windSpeed.push_back(
          double(boost::lexical_cast<double>(columns[21])));  // [m/s]

      double fc = double(boost::lexical_cast<double>(columns[22]));

      totalSkyCover.push_back(fc);  // [tenths]

      opaqueSkyCover.push_back(
          double(boost::lexical_cast<double>(columns[23])));  // [tenths]

      double eSky = 0.787 + 0.764*log(Tdp/Tdb)*
          (1 + 0.0224*fc - 0.0035*pow(fc,2) + 0.00028*pow(fc,3));

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

      // clockwise from north
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

  // Calculate average daily and monthly temperatures
  dailyAverageTemperatures.reserve(365);
  monthlyAverageTemperatures.reserve(12);

  double dailyTemperatureSum = 0;
  double monthlyTemperatureSum = 0;

  int dayCount = 0;
  int monthCount = 0;

  for (std::size_t h = 0; h < days.size(); ++h)
  {
    dailyTemperatureSum += dryBulbTemp[h];
    monthlyTemperatureSum += dryBulbTemp[h];
    dayCount += 1;
    monthCount += 1;

    if (h != days.size() - 1)
    {
      if (days[h] != days[h+1])
      {
        dailyAverageTemperatures.push_back(dailyTemperatureSum/dayCount);
        dailyTemperatureSum = 0;
        dayCount = 0;
      }
      if (months[h] != months[h+1])
      {
        monthlyAverageTemperatures.push_back(monthlyTemperatureSum/monthCount);
        monthlyTemperatureSum = 0;
        monthCount = 0;
      }
    }
    else
    {
      dailyAverageTemperatures.push_back(dailyTemperatureSum/dayCount);
      dailyTemperatureSum = 0;
      dayCount = 0;
      monthlyAverageTemperatures.push_back(monthlyTemperatureSum/monthCount);
      monthlyTemperatureSum = 0;
      monthCount = 0;
    }
  }

  minimumAverageMontlyTemperature = *min_element(monthlyAverageTemperatures.begin(),monthlyAverageTemperatures.end());
  maximumAverageMontlyTemperature = *max_element(monthlyAverageTemperatures.begin(),monthlyAverageTemperatures.end());

}

#endif
