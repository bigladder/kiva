#ifndef Domain_HPP
#define Domain_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include "Input.h"
#include "Mesher.h"

namespace blas = boost::numeric::ublas;

class Domain
{
public:

		// mesher
		Mesher mesher;
		size_t nR;
		size_t nZ;

		// inherent properties
		blas::matrix<double> rho;
		blas::matrix<double> cp;
		blas::matrix<double> k;

		// derived properties
		blas::matrix<double> theta;
		blas::matrix<double> a;
		blas::matrix<double> b;
		blas::matrix<double> c;
		blas::matrix<double> d;
		blas::matrix<double> e;
		blas::matrix<double> f;

		// organizational properties
		enum CellType
		{
			EXTERIOR_AIR,
			EXTERIOR_GRADE,
			EXTERIOR_WALL,
			INTERIOR_AIR,
			INTERIOR_SLAB,
			INTERIOR_WALL,
			INTERIOR_INSULATION_EDGE, // The bottom edge of partial wall insulation
			WALL_TOP,
			AXIS,  // TODO: Change to SYMMETRY
			FAR_FIELD,
			DEEP_GROUND,
			NORMAL,
			ZERO_WIDTH_R,
			ZERO_WIDTH_Z,
			ZERO_WIDTH_RZ
		};
		blas::matrix<CellType> cellType;
		size_t slabJ;
		size_t slabImin;
		size_t slabImax;

public:

		Domain();
		Domain(Mesher &mesher, Foundation &foundation, double &timestep);
		double getDRP(size_t i);
		double getDRM(size_t i);
		double getDZP(size_t j);
		double getDZM(size_t j);
		double getKRP(size_t i,size_t j);
		double getKRM(size_t i,size_t j);
		double getKZP(size_t i,size_t j);
		double getKZM(size_t i,size_t j);
		void printCellTypes();

};


#endif
