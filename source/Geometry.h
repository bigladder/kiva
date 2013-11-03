/* Geometry.c++ is part of Kiva (Written by Neal Kruis)
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

#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/multi/geometries/multi_polygon.hpp>
#include <boost/geometry/multi/geometries/multi_point.hpp>

#include "Functions.h"


typedef boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian> Point;
typedef boost::geometry::model::polygon<Point, true, false> Polygon;
typedef boost::geometry::model::ring<Point, true, false> Ring;
typedef boost::geometry::model::multi_polygon<Polygon> MultiPolygon;
typedef boost::geometry::model::multi_point<Point> MultiPoint;
typedef boost::geometry::model::linestring<Point> Line;

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

bool comparePointsX(Point first, Point second);
bool comparePointsY(Point first, Point second);




#endif /* GEOMETRY_H_ */
