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
	double frac = hour - hour_before;

	if (boost::gregorian::gregorian_calendar::is_leap_year(year))
	{

		if (hour < 1416.0)
		{
			return (*this)[hour_before]*(1.0 - frac) + (*this)[hour_after]*frac;

		}
		else
		{
			return (*this)[hour_before - 24]*(1.0 - frac)
					+ (*this)[hour_after - 24]*frac;
		}
	}
	else
	{
		return (*this)[hour_before]*(1.0 - frac) + (*this)[hour_after]*frac;

	}
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
    	}

    	if (row > 8)
    	{
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

    		globalHorizontalSolar.push_back(
    				double(boost::lexical_cast<double>(columns[13])));  // [W/m2]

    		directNormalSolar.push_back(
    				double(boost::lexical_cast<double>(columns[14])));  // [W/m2]

    		diffuseHorizontalSolar.push_back(
    				double(boost::lexical_cast<double>(columns[15])));  // [W/m2]

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

    	}
    	columns.clear();

    }

    inf.close();
}

#endif
