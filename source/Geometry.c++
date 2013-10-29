/*
 * Geometry.c++
 *
 *  Created on: Oct 28, 2013
 *      Author: nkruis
 */

#ifndef GEOMETRY_CPP_
#define GEOMETRY_CPP_

#include "Geometry.h"

bool isRectilinear(Polygon poly)
{
	for (std::size_t v = 0; v < poly.outer().size(); v++)
	{
		double x = poly.outer()[v].get<0>();
		double y = poly.outer()[v].get<1>();
		double xNext, yNext;
		if (v == poly.outer().size() - 1)
		{
			xNext = poly.outer()[0].get<0>();
			yNext = poly.outer()[0].get<1>();
		}
		else
		{
			xNext = poly.outer()[v+1].get<0>();
			yNext = poly.outer()[v+1].get<1>();
		}

		if (isEqual(x,xNext) || isEqual(y,yNext))
		{
			// Do nothing
		}
		else
		{
			return false;
		}
	}
	return true;
}

Polygon offset(Polygon poly, double dist)
{
	if (!isRectilinear(poly))
	{
		// Throw exception
	}
	// March around polygon set offsets from each vertex
	int nV = poly.outer().size();


	Polygon offset;
	// Main loop
	for (int v = 0; v < nV; v++)
	{
		// current coordinates
		double x = poly.outer()[v].get<0>();
		double y = poly.outer()[v].get<1>();

		double xNew, yNew;

		switch (getDirectionOut(poly,v))
		{
		case geom::Y_POS:
			if (getTurn(poly,v) == geom::LEFT)
			{
				xNew = x - dist;
				yNew = y + dist;
			}
			else
			{
				xNew = x - dist;
				yNew = y - dist;
			}
			break;
		case geom::X_POS:
			if (getTurn(poly,v) == geom::LEFT)
			{
				xNew = x + dist;
				yNew = y + dist;
			}
			else
			{
				xNew = x - dist;
				yNew = y + dist;
			}
			break;
		case geom::Y_NEG:
			if (getTurn(poly,v) == geom::LEFT)
			{
				xNew = x + dist;
				yNew = y - dist;
			}
			else
			{
				xNew = x + dist;
				yNew = y + dist;
			}
			break;
		case geom::X_NEG:
			if (getTurn(poly,v) == geom::LEFT)
			{
				xNew = x - dist;
				yNew = y - dist;
			}
			else
			{
				xNew = x + dist;
				yNew = y - dist;
			}
			break;
		}

		offset.outer().push_back(Point(xNew, yNew));

	}
	return offset;
}

geom::Direction getDirectionIn(Polygon poly, std::size_t vertex)
{
	if (!isRectilinear(poly))
	{
		// Throw exception
	}

	double xPrev;
	double yPrev;
	double x;
	double y;

	std::size_t nV = poly.outer().size();

	if (vertex == 0)
	{
		xPrev = poly.outer()[nV - 1].get<0>();
		yPrev = poly.outer()[nV - 1].get<1>();
	}
	else
	{
		xPrev = poly.outer()[vertex - 1].get<0>();
		yPrev = poly.outer()[vertex - 1].get<1>();
	}

	x = poly.outer()[vertex].get<0>();
	y = poly.outer()[vertex].get<1>();

	if (isLessThan(x,xPrev))
	{
		return geom::X_NEG;
	}
	else if (isGreaterThan(x,xPrev))
	{
		return geom::X_POS;
	}
	else if (isLessThan(y,yPrev))
	{
		return geom::Y_NEG;
	}
	else // if (isGreaterThan(y,yPrev))
	{
		return geom::Y_POS;
	}

}

geom::Direction getDirectionOut(Polygon poly, std::size_t vertex)
{
	if (!isRectilinear(poly))
	{
		// Throw exception
	}

	double xNext;
	double yNext;
	double x;
	double y;

	std::size_t nV = poly.outer().size();

	if (vertex == nV - 1)
	{
		xNext = poly.outer()[0].get<0>();
		yNext = poly.outer()[0].get<1>();
	}
	else
	{
		xNext = poly.outer()[vertex + 1].get<0>();
		yNext = poly.outer()[vertex + 1].get<1>();
	}

	x = poly.outer()[vertex].get<0>();
	y = poly.outer()[vertex].get<1>();

	if (isLessThan(xNext,x))
	{
		return geom::X_NEG;
	}
	else if (isGreaterThan(xNext,x))
	{
		return geom::X_POS;
	}
	else if (isLessThan(yNext,y))
	{
		return geom::Y_NEG;
	}
	else // if (isGreaterThan(yNext,y))
	{
		return geom::Y_POS;
	}

}

geom::Turn getTurn(Polygon poly, std::size_t vertex)
{
	switch (getDirectionIn(poly,vertex))
	{
	case geom::X_NEG:
		{
		if (getDirectionOut(poly,vertex) == geom::Y_POS)
			return geom::RIGHT;
		else  // if (getDirectionOut(poly,vertex) == Y_NEG)
			return geom::LEFT;
		}
		break;
	case geom::X_POS:
		{
		if (getDirectionOut(poly,vertex) == geom::Y_POS)
			return geom::LEFT;
		else  // if (getDirectionOut(poly,vertex) == Y_NEG)
			return geom::RIGHT;
		}
		break;
	case geom::Y_NEG:
		{
		if (getDirectionOut(poly,vertex) == geom::X_POS)
			return geom::LEFT;
		else  // if (getDirectionOut(poly,vertex) == X_NEG)
			return geom::RIGHT;
		}
		break;
	default:  //case geom::Y_POS:
		{
		if (getDirectionOut(poly,vertex) == geom::X_POS)
			return geom::RIGHT;
		else  // if (getDirectionOut(poly,vertex) == X_NEG)
			return geom::LEFT;
		}
		break;
	}
}

double getXmax(Polygon poly, std::size_t vertex)
{
	double x = poly.outer()[vertex].get<0>();
	double xNext;
	std::size_t nV = poly.outer().size();

	if (vertex == nV - 1)
	{
		xNext = poly.outer()[0].get<0>();
	}
	else
	{
		xNext = poly.outer()[vertex + 1].get<0>();
	}

	return std::max(x,xNext);
}

double getYmax(Polygon poly, std::size_t vertex)
{
	double y = poly.outer()[vertex].get<1>();
	double yNext;
	std::size_t nV = poly.outer().size();

	if (vertex == nV - 1)
	{
		yNext = poly.outer()[0].get<1>();
	}
	else
	{
		yNext = poly.outer()[vertex + 1].get<1>();
	}

	return std::max(y,yNext);
}

double getXmin(Polygon poly, std::size_t vertex)
{
	double x = poly.outer()[vertex].get<0>();
	double xNext;
	std::size_t nV = poly.outer().size();

	if (vertex == nV - 1)
	{
		xNext = poly.outer()[0].get<0>();
	}
	else
	{
		xNext = poly.outer()[vertex + 1].get<0>();
	}

	return std::min(x,xNext);
}

double getYmin(Polygon poly, std::size_t vertex)
{
	double y = poly.outer()[vertex].get<1>();
	double yNext;
	std::size_t nV = poly.outer().size();

	if (vertex == nV - 1)
	{
		yNext = poly.outer()[0].get<1>();
	}
	else
	{
		yNext = poly.outer()[vertex + 1].get<1>();
	}

	return std::min(y,yNext);
}

#endif /* GEOMETRY_CPP_ */



