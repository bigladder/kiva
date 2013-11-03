/* GroundPlot.h is part of Kiva (Written by Neal Kruis)
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

#ifndef GROUNDPLOT_H_
#define GROUNDPLOT_H_

#include <mgl2/mgl.h>

#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

#include "Input.h"
#include "Domain.h"
#include "Functions.h"
#include "Geometry.h"

struct Axis
{
	std::size_t nMin;
	std::size_t nMax;
	std::size_t nN;
	Mesher mesh;
};

enum SliceType
{
	XZ_2D,
	XY,
	XZ,
	YZ
};

class GroundPlot
{
private:
	mglData TDat, hDat, vDat, cDat, hGrid, vGrid, TGrid;

	std::size_t iMin, iMax, jMin, jMax, kMin, kMax;

	double slice;
	double xMin, xMax, yMin, yMax;

	Axis hAxis;
	Axis vAxis;


	double nextPlotTime;

public:

	mglGraph gr;
	OutputAnimation outputAnimation;
	std::vector<Block> blocks;
	SliceType sliceType;

	double tStart, tEnd;
	GroundPlot(OutputAnimation &outputAnimation, Domain &domain, std::vector<Block> &blocks);
	void createFrame(boost::multi_array<double, 3> &T, std::string timeStamp);
	bool makeNewFrame(double tNow);
};


#endif /* GROUNDPLOT_H_ */
