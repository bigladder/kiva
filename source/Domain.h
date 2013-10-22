/* Domain.h is part of Kiva (Written by Neal Kruis)
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

#ifndef Domain_HPP
#define Domain_HPP

#include "Input.h"
#include "Mesher.h"
#include "Functions.h"

#include <fstream>

#include <boost/multi_array.hpp>

class Cell
{
public:

	// inherent properties
	double density;
	double specificHeat;
	double conductivity;

	// derived properties
	double cxp_c;
	double cxm_c;
	double cxp;
	double cxm;
	double cyp;
	double cym;
	double czp;
	double czm;

	// organizational properties
	enum CellType
	{
		EXTERIOR_AIR,  // 0
		INTERIOR_AIR,  // 1
		NORMAL,  // 2
		BOUNDARY,  // 3
		ZERO_THICKNESS  // 4
	};
	CellType cellType;

	size_t blockNumber;
	Block block;

	size_t surfaceNumber;
	Surface surface;
};

class Domain
{
public:

		// mesh
		Mesher meshX;
		Mesher meshY;
		Mesher meshZ;
		std::size_t nX;
		std::size_t nY;
		std::size_t nZ;

		boost::multi_array<Cell, 3> cell;
		
public:

		Domain();
		Domain(Foundation &foundation);
		void setDomain(Foundation &foundation);
		double getDXP(size_t i);
		double getDXM(size_t i);
		double getDYP(size_t j);
		double getDYM(size_t j);
		double getDZP(size_t k);
		double getDZM(size_t k);
		double getKXP(size_t i,size_t j,size_t k);
		double getKXM(size_t i,size_t j,size_t k);
		double getKYP(size_t i,size_t j,size_t k);
		double getKYM(size_t i,size_t j,size_t k);
		double getKZP(size_t i,size_t j,size_t k);
		double getKZM(size_t i,size_t j,size_t k);
		void printCellTypes();

};


#endif
