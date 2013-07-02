/*
 * Main.cpp
 *
 *  Created on: Nov 6, 2012
 *      Author: nkruis
 */

#include <iostream>
#include <fstream>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include "InputParser.h"
#include "WeatherData.h"
#include "Input.h"
#include "Ground.h"
#include <cmath>

using namespace std;
using namespace boost::posix_time;
using namespace boost::gregorian;

int main()
{

	ptime beginCalc = second_clock::local_time();
	cout << "Starting Program: " << beginCalc << endl;

	// parse input
	Input input = inputParser();

	// parse weather
	WeatherData weather(input.simulationControl.weatherFile);

	// simulation
	ptime simStart(input.simulationControl.startDate);
	ptime simEnd(input.simulationControl.endDate + days(1));
	time_duration simDuration =  simEnd - simStart;


	double tstart = 0.0; // [s] Simulation start time
	double tend = simDuration.total_seconds(); // [s] Simulation end time
	double timestep = input.simulationControl.timestep.total_seconds();

	// set up output file
	ofstream output;
	output.open("Output.csv");
	output << "Time Stamp, Heat Flux [W/m2]" << endl;

	// initialize
	Ground ground(weather,input.foundations[0],input.simulationControl);

	// loop

	ptime prevTime = second_clock::local_time();

	for (double t = tstart; t < tend; t = t + timestep)
	{
		ptime currentTime = second_clock::local_time();
		ptime simTime(input.simulationControl.startDate,seconds(t));
		date today = simTime.date();

		ground.calculate(t);

		if (currentTime - prevTime > milliseconds(500))
		{

			double percentComplete = round(t/tend*1000)/10.0;

			cout << percentComplete << "% (" << today << ")\n";

			prevTime = currentTime;
		}

		output << to_simple_string(simTime) << ", " <<
				  ground.QSlabTotal << endl;

	}

	output.close();

	ptime simTime(input.simulationControl.startDate,seconds(tend));
	date today = simTime.date();
	cout << "100% (" << today << ")\n";

	// process output
	//cout << "Total Heat Flux: "<< ground.QSlabTotal << endl;


	ptime finishCalc = second_clock::local_time();
	cout << "Finished Program: " << finishCalc << endl;

	time_duration totalCalc = finishCalc - beginCalc;
	cout << "Elapsed Time: " << totalCalc << endl;

	return 0;
}



