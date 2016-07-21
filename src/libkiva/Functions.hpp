/* Copyright (c) 2012-2016 Big Ladder Software. All rights reserved.
* See the LICENSE file for additional terms and conditions. */

#ifndef Functions_HPP
#define Functions_HPP

#include <fstream>
#include <math.h>
#include <algorithm>
#include <vector>

bool isLessThan(double first, double second);
bool isLessOrEqual(double first, double second);
bool isEqual(double first, double second);
bool isEqual(double first, double second, double epsilon);
bool isGreaterThan(double first, double second);
bool isGreaterOrEqual(double first, double second);
bool isEven(int N);
bool isOdd(int N);
void solveTDM(const std::vector<double>& a1, const std::vector<double>& a2,
          std::vector<double>& a3, std::vector<double>& b,
              std::vector<double>& x);
std::istream& safeGetline(std::istream& is, std::string& t);

#endif
