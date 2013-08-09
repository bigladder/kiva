/*
 * WeatherData.cpp
 *
 *  Created on: Nov 7, 2012
 *      Author: nkruis
 */

#ifndef WeatherData_CPP
#define WeatherData_CPP

#include <fstream>
#include "WeatherData.h"
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <cmath>

using namespace std;
using namespace boost;

static const double PI = 4.0*atan(1.0);

double HourlyData::getValue(ptime t)
{
	long year = t.date().year();
	date newYearDate(year,Jan,1);
	ptime newYear(newYearDate);
	time_duration duration = t - newYear;
	long second = duration.total_seconds();
	double hour = second/60.0/60.0;
	long hour_before = hour;
	long hour_after = hour_before + 1;
	double frac = hour - hour_before;

	if (gregorian_calendar::is_leap_year(year))
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

WeatherData::WeatherData(string weatherFile)
{
	importEPW(weatherFile);
}

void WeatherData::importEPW(string epwFile)
{

	ifstream inf;
	inf.open(epwFile.c_str());

    if (!inf)
    {
        // Print an error and exit
        cerr << "Unable to read EPW file" << endl;
        exit(1);
    }

    // While there's still stuff left to read
    string line;
    vector <string> columns;

    typedef tokenizer< escaped_list_separator<char> > Tokenizer;

    int row = 0;

    while (getline(inf,line))
    {
    	row += 1;
    	Tokenizer tok(line, escaped_list_separator<char>("\\",",","\""));

    	for (tokenizer< escaped_list_separator<char> >::iterator
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
    		double Tdb = double(lexical_cast<double>(columns[6])) +
		              273.15;

    		dryBulbTemp.push_back(Tdb);  // [K]

    		double Tdp = double(lexical_cast<double>(columns[7])) +
		              273.15;

    		dewPointTemp.push_back(Tdp);  // [K]

    		relativeHumidity.push_back(
    				double(lexical_cast<double>(columns[8]))/100.0);  // [frac]

    		atmosphericPressure.push_back(
    				double(lexical_cast<double>(columns[9])));  // [Pa]

    		globalHorizontalSolar.push_back(
    				double(lexical_cast<double>(columns[13])));  // [W/m2]

    		directNormalSolar.push_back(
    				double(lexical_cast<double>(columns[14])));  // [W/m2]

    		diffuseHorizontalSolar.push_back(
    				double(lexical_cast<double>(columns[15])));  // [W/m2]

    		windDirection.push_back(double(lexical_cast<double>(columns[20]))*
    					PI/180.0);  // [rad]

    		windSpeed.push_back(
    				double(lexical_cast<double>(columns[21])));  // [m/s]

    		double fc = double(lexical_cast<double>(columns[22]));

    		totalSkyCover.push_back(fc);  // [tenths]

    		opaqueSkyCover.push_back(
    				double(lexical_cast<double>(columns[23])));  // [tenths]

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
