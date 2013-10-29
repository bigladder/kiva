/*
 * Geometry.h
 *
 *  Created on: Oct 28, 2013
 *      Author: nkruis
 */

#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>

#include "Functions.h"


typedef boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian> Point;
typedef boost::geometry::model::polygon<Point, true, false> Polygon;
typedef boost::geometry::model::ring<Point, true, false> Ring;

namespace geom
{
	enum Direction
	{
		X_NEG,
		X_POS,
		Y_NEG,
		Y_POS
	};

	enum Turn
	{
		LEFT,
		RIGHT
	};
};

bool isRectilinear(Polygon poly);
Polygon offset(Polygon poly, double dist);
geom::Direction getDirectionIn(Polygon poly, std::size_t vertex);
geom::Direction getDirectionOut(Polygon poly, std::size_t vertex);
geom::Turn getTurn(Polygon poly, std::size_t vertex);
double getXmin(Polygon poly, std::size_t vertex);
double getYmin(Polygon poly, std::size_t vertex);
double getXmax(Polygon poly, std::size_t vertex);
double getYmax(Polygon poly, std::size_t vertex);



#endif /* GEOMETRY_H_ */
