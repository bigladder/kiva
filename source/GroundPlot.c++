/* GroundPlot.c++ is part of Kiva (Written by Neal Kruis)
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

#ifndef GroundPlot_CPP
#define GroundPlot_CPP

#include "GroundPlot.h"

GroundPlot::GroundPlot(OutputAnimation &outputAnimation, Domain &domain, std::vector<Block> &blocks) :
	outputAnimation(outputAnimation), blocks(blocks)
{

	boost::filesystem::remove_all(outputAnimation.name + "_frames");
	boost::filesystem::create_directory(outputAnimation.name + "_frames");

	nextPlotTime = 0.0;
	frameNumber = 0;

	std::size_t contourLevels = 13;

	if (isEqual(outputAnimation.xRange.first,outputAnimation.xRange.second))
	{
		iMin = domain.meshX.getNearestIndex(outputAnimation.xRange.first);
		iMax = domain.meshX.getNearestIndex(outputAnimation.xRange.second);
	}
	else
	{
		iMin = domain.meshX.getPreviousIndex(outputAnimation.xRange.first);
		iMax = domain.meshX.getNextIndex(outputAnimation.xRange.second);

		// Check for exact match
		if (isEqual(domain.meshX.centers[iMin + 1],outputAnimation.xRange.first))
			iMin += 1;
		if (isEqual(domain.meshX.centers[iMax - 1],outputAnimation.xRange.second))
			iMax -= 1;
	}

	if (isEqual(outputAnimation.yRange.first,outputAnimation.yRange.second))
	{
		jMin = domain.meshY.getNearestIndex(outputAnimation.yRange.first);
		jMax = domain.meshY.getNearestIndex(outputAnimation.yRange.second);
	}
	else
	{
		jMin = domain.meshY.getPreviousIndex(outputAnimation.yRange.first);
		jMax = domain.meshY.getNextIndex(outputAnimation.yRange.second);

		// Check for exact match
		if (isEqual(domain.meshY.centers[jMin + 1],outputAnimation.yRange.first))
			jMin += 1;
		if (isEqual(domain.meshY.centers[jMax - 1],outputAnimation.yRange.second))
			jMax -= 1;
	}

	if (isEqual(outputAnimation.zRange.first,outputAnimation.zRange.second))
	{
		kMin = domain.meshZ.getNearestIndex(outputAnimation.zRange.first);
		kMax = domain.meshZ.getNearestIndex(outputAnimation.zRange.second);
	}
	else
	{
		kMin = domain.meshZ.getPreviousIndex(outputAnimation.zRange.first);
		kMax = domain.meshZ.getNextIndex(outputAnimation.zRange.second);

		// Check for exact match
		if (isEqual(domain.meshZ.centers[kMin + 1],outputAnimation.zRange.first))
			kMin += 1;
		if (isEqual(domain.meshZ.centers[kMax - 1],outputAnimation.zRange.second))
			kMax -= 1;
	}

	std::size_t nI = iMax - iMin + 1;
	std::size_t nJ = jMax - jMin + 1;
	std::size_t nK = kMax - kMin + 1;

	xMin = domain.meshX.dividers[0];
	xMax = domain.meshX.dividers[domain.nX];

	yMin = domain.meshY.dividers[0];
	yMax = domain.meshY.dividers[domain.nY];


	if (nI == 1)
	{
		sliceType = YZ;

		hAxis.nN = nJ;
		hAxis.nMin = jMin;
		hAxis.nMax = jMax;
		hAxis.mesh = domain.meshY;

		vAxis.nN = nK;
		vAxis.nMin = kMin;
		vAxis.nMax = kMax;
		vAxis.mesh = domain.meshZ;

		slice = domain.meshX.centers[iMin];
	}
	else if (nJ == 1)
	{
		if (domain.meshY.centers.size() == 1)
			sliceType = XZ_2D;
		else
			sliceType = XZ;

		hAxis.nN = nI;
		hAxis.nMin = iMin;
		hAxis.nMax = iMax;
		hAxis.mesh = domain.meshX;

		vAxis.nN = nK;
		vAxis.nMin = kMin;
		vAxis.nMax = kMax;
		vAxis.mesh = domain.meshZ;

		slice = domain.meshY.centers[jMin];
	}
	else if (nK == 1)
	{
		sliceType = XY;

		hAxis.nN = nI;
		hAxis.nMin = iMin;
		hAxis.nMax = iMax;
		hAxis.mesh = domain.meshX;

		vAxis.nN = nJ;
		vAxis.nMin = jMin;
		vAxis.nMax = jMax;
		vAxis.mesh = domain.meshY;

		slice = domain.meshZ.centers[kMin];
	}

	mglData TDatRef(hAxis.nN,vAxis.nN),
			hDatRef(hAxis.nN),
			vDatRef(vAxis.nN),
			hGridRef(hAxis.nN + 1),
			vGridRef(vAxis.nN + 1),
			TGridRef(hAxis.nN + 1, vAxis.nN + 1),
			cDatRef(contourLevels);

	TDat = TDatRef;
	hDat = hDatRef;
	vDat = vDatRef;
	hGrid = hGridRef;
	vGrid = vGridRef;
	TGrid = TGridRef;
	cDat = cDatRef;

	hGrid.a[0] = hAxis.mesh.centers[hAxis.nMin];

	for(size_t i = 0; i < hAxis.nN - 1; i++)
	{
		hDat.a[i] = hAxis.mesh.centers[i + hAxis.nMin];
		hGrid.a[i + 1] = hAxis.mesh.dividers[i + hAxis.nMin + 1];
	}

	hDat.a[hAxis.nN - 1] = hAxis.mesh.centers[hAxis.nMax];
	hGrid.a[hAxis.nN] = hAxis.mesh.centers[hAxis.nMax];

	vGrid.a[0] = vAxis.mesh.centers[vAxis.nMin];

	for(size_t j = 0; j < vAxis.nN - 1; j++)
	{
		vDat.a[j] = vAxis.mesh.centers[j + vAxis.nMin];
		vGrid.a[j + 1] = vAxis.mesh.dividers[j + vAxis.nMin + 1];
	}

	vDat.a[vAxis.nN - 1] = vAxis.mesh.centers[vAxis.nMax];
	vGrid.a[vAxis.nN] = vAxis.mesh.centers[vAxis.nMax];

	for(size_t j = 0; j <= vAxis.nN; j++)
	{
		for(size_t i = 0; i <= hAxis.nN; i++)
		{
			TGrid.a[i+hAxis.nN*j] = 200.0;
		}
	}

	for (size_t n = 0; n < contourLevels; n++)
	{
		double mid = 50.0;
		double range = 120.0;
		double step = range / (contourLevels - 1);
		cDat.a[n] = mid - range/2 + double(n)*step;
	}

}

void GroundPlot::createFrame(boost::multi_array<double, 3> &T, std::string timeStamp)
{

	std::size_t nI = iMax - iMin + 1;
	std::size_t nJ = jMax - jMin + 1;

	for(size_t k = kMin; k <= kMax; k++)
	{
		for(size_t j = jMin; j <= jMax; j++)
		{
			for(size_t i = iMin; i <= iMax; i++)
			{
				std::size_t index = (i-iMin)+nI*(j-jMin)+nI*nJ*(k-kMin);
				double TinF = (T[i][j][k] - 273.15)*9/5 + 32.0;
				TDat.a[index] = TinF;
			}
		}
	}

	double hMin = hGrid.a[0];
	double hMax = hGrid.a[hAxis.nN];
	double hRange = hMax - hMin;

	double vMin = vGrid.a[0];
	double vMax = vGrid.a[vAxis.nN];
	double vRange = vMax - vMin;

	int nT = cDat.GetNN();
	double Tmin = cDat.a[0];
	double Tmax = cDat.a[nT - 1];
	double Tstep = cDat.a[1] - cDat.a[0];

	// Text properties
	double hText = 0.05;
	double vText = 0.95;

	double vTextSpacing = 0.05;

	mglGraph gr;

	// Plot
	gr.Clf(1,1,1);
	double aspect = 1.0;
	int height = outputAnimation.size;
	int width = height*aspect;

	gr.SetSize(width,height);

	gr.LoadFont("none");
	gr.SetFontSize(2.0);
	gr.SetRange('x', hGrid);
	gr.SetRange('y', vGrid);
	gr.SetRange('c', Tmin, Tmax);
	gr.SetRange('z', Tmin, Tmax);
	gr.SetTicks('c', Tstep, nT, Tmin);
	gr.Aspect(hRange, vRange);
	gr.Axis("yU");
	gr.Axis("x");
	gr.Colorbar("_");
	gr.Box("k",false);
	gr.Dens(hDat, vDat, TDat);
	if (outputAnimation.contours)
		gr.Cont(cDat, hDat, vDat, TDat,"H");
	if (outputAnimation.gradients)
		gr.Grad(hDat, vDat, TDat);
	if (outputAnimation.grid)
		gr.Grid(hGrid, vGrid, TGrid, "W");

	// Draw blocks
	for (size_t b = 0; b < blocks.size(); b++)
	{
		switch (sliceType)
		{
		case XZ_2D:
			{
			mglPoint bl = mglPoint(std::min(std::max(blocks[b].xMin, hMin),hMax),
								   std::min(std::max(blocks[b].zMin, vMin),vMax),
								   210.0);
			mglPoint br = mglPoint(std::max(std::min(blocks[b].xMax, hMax),hMin),
								   std::min(std::max(blocks[b].zMin, vMin),vMax),
								   210.0);
			mglPoint tr = mglPoint(std::max(std::min(blocks[b].xMax, hMax),hMin),
								   std::max(std::min(blocks[b].zMax, vMax),vMin),
								   210.0);
			mglPoint tl = mglPoint(std::min(std::max(blocks[b].xMin, hMin),hMax),
								   std::max(std::min(blocks[b].zMax, vMax),vMin),
								   210.0);

			gr.Line(bl, br, "k");
			gr.Line(br, tr, "k");
			gr.Line(tr, tl, "k");
			gr.Line(tl, bl, "k");
			}
			break;
		case XY:
			{
			// Find intersection with viewing window
			Polygon view;
			view.outer().push_back(Point(hMin,vMin));
			view.outer().push_back(Point(hMin,vMax));
			view.outer().push_back(Point(hMax,vMax));
			view.outer().push_back(Point(hMax,vMin));

			MultiPolygon intersection;
			boost::geometry::intersection(view, blocks[b].polygon, intersection);

			// loop through each polygon in resulting multi_polygon
			for (std::size_t p = 0; p < intersection.size(); p++)
			{
				// loop through each vertex in each polygon and create a line
				std::size_t nV = intersection[p].outer().size();
				for (std::size_t v = 0; v < nV - 1; v++)
				{
					gr.Line(mglPoint(intersection[p].outer()[v].get<0>(),intersection[p].outer()[v].get<1>(),210.0),
							mglPoint(intersection[p].outer()[v+1].get<0>(),intersection[p].outer()[v+1].get<1>(),210.0),
							"k");
				}
				gr.Line(mglPoint(intersection[p].outer()[nV - 1].get<0>(),intersection[p].outer()[nV - 1].get<1>(),210.0),
						mglPoint(intersection[p].outer()[0].get<0>(),intersection[p].outer()[0].get<1>(),210.0),
						"k");
			}

			std::string sliceString = "Z = " + str(boost::format("%0.2f") % slice);
			gr.Puts(hText, vText - vTextSpacing, sliceString.c_str(), ":AL");

			}
			break;
		case XZ:
			{
			// Find intersecting point pairs with slicing plane
			Line slicingPlane;
			slicingPlane.push_back(Point(xMin,slice));
			slicingPlane.push_back(Point(xMax,slice));

			MultiPoint intersection;
			boost::geometry::intersection(slicingPlane, blocks[b].polygon, intersection);

			// sort points in ascending order
			sort(intersection.begin(), intersection.end(), comparePointsX);

			for (std::size_t p = 0; p < intersection.size(); p++)
			{
				// Use point pairs and zmin/max to draw rectangles similar to 2D
				// case.
				double p1 = intersection[p].get<0>();
				double p2 = intersection[p+1].get<0>();


				mglPoint bl = mglPoint(std::min(std::max(p1, hMin),hMax),
									   std::min(std::max(blocks[b].zMin, vMin),vMax),
									   210.0);
				mglPoint br = mglPoint(std::max(std::min(p2, hMax),hMin),
									   std::min(std::max(blocks[b].zMin, vMin),vMax),
									   210.0);
				mglPoint tr = mglPoint(std::max(std::min(p2, hMax),hMin),
									   std::max(std::min(blocks[b].zMax, vMax),vMin),
									   210.0);
				mglPoint tl = mglPoint(std::min(std::max(p1, hMin),hMax),
									   std::max(std::min(blocks[b].zMax, vMax),vMin),
									   210.0);

				gr.Line(bl, br, "k");
				gr.Line(br, tr, "k");
				gr.Line(tr, tl, "k");
				gr.Line(tl, bl, "k");

				p += 1; // skip one point, on to the next pair
			}

			std::string sliceString = "Y = " + str(boost::format("%0.2f") % slice);
			gr.Puts(hText, vText - vTextSpacing, sliceString.c_str(), ":AL");

			}
			break;
		case YZ:
			{
			// Find intersecting point pairs with slicing plane
			Line slicingPlane;
			slicingPlane.push_back(Point(slice,yMin));
			slicingPlane.push_back(Point(slice,yMax));

			MultiPoint intersection;
			boost::geometry::intersection(slicingPlane, blocks[b].polygon, intersection);

			// sort points in ascending order
			sort(intersection.begin(), intersection.end(), comparePointsY);

			for (std::size_t p = 0; p < intersection.size(); p++)
			{
				// Use point pairs and zmin/max to draw rectangles similar to 2D
				// case.
				double p1 = intersection[p].get<1>();
				double p2 = intersection[p+1].get<1>();


				mglPoint bl = mglPoint(std::min(std::max(p1, hMin),hMax),
									   std::min(std::max(blocks[b].zMin, vMin),vMax),
									   210.0);
				mglPoint br = mglPoint(std::max(std::min(p2, hMax),hMin),
									   std::min(std::max(blocks[b].zMin, vMin),vMax),
									   210.0);
				mglPoint tr = mglPoint(std::max(std::min(p2, hMax),hMin),
									   std::max(std::min(blocks[b].zMax, vMax),vMin),
									   210.0);
				mglPoint tl = mglPoint(std::min(std::max(p1, hMin),hMax),
									   std::max(std::min(blocks[b].zMax, vMax),vMin),
									   210.0);

				gr.Line(bl, br, "k");
				gr.Line(br, tr, "k");
				gr.Line(tr, tl, "k");
				gr.Line(tl, bl, "k");

				p += 1; // skip one point, on to the next pair
			}

			std::string sliceString = "X = " + str(boost::format("%0.2f") % slice);
			gr.Puts(hText, vText - vTextSpacing, sliceString.c_str(), ":AL");

			}
			break;
		}
	}

	// Timestamp
	gr.Puts(hText, vText, timeStamp.c_str(), ":AL");


	gr.WritePNG((outputAnimation.name + "_frames/" + outputAnimation.name + str(boost::format("%04d") % frameNumber) + ".png").c_str(),"",false);

	frameNumber += 1;
	nextPlotTime += outputAnimation.frequency.total_seconds();
}

bool GroundPlot::makeNewFrame(double tNow)
{
	if (isGreaterOrEqual(tNow, nextPlotTime) &&
		isGreaterOrEqual(tNow, tStart) &&
		isLessOrEqual(tNow, tEnd))
		return true;
	else
		return false;
}
#endif


