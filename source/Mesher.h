#ifndef Mesher_HPP
#define Mesher_HPP

#include <boost/numeric/ublas/vector.hpp>
namespace blas = boost::numeric::ublas;

class Interval
{
public:

	double maxGrowthCoeff;
	double minCellDim;
	enum Growth
	{
		FORWARD,
		BACKWARD,
		UNIFORM
	};
	Growth growthDir;

};

class MeshData
{
public:
	std::vector<double> points;
	std::vector<Interval> intervals;

};

class Mesher
{
private:
		MeshData xData;
		MeshData yData;

public:
		blas::vector<double> xDividers, yDividers;
		blas::vector<double> xDeltas, yDeltas;
		blas::vector<double> xCenters, yCenters;

public:

		Mesher();
		Mesher(MeshData &xData, MeshData &yData);
		void makeMesh(MeshData data, blas::vector<double> &dividers,
				      blas::vector<double> &deltas,
				      blas::vector<double> &centers);

};


#endif
