/* Functions.c++ is part of Kiva (Written by Neal Kruis)
 * Copyright (C) 2012-2015 Big Ladder Software <info@bigladdersoftware.com>
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

#ifndef Functions_CPP
#define Functions_CPP

#include "Functions.h"

static const double EPSILON = 1E-7;

bool isLessThan(double first, double second)
{
  if(first - second < -EPSILON)
    return true;
  else
    return false;
}

bool isLessOrEqual(double first, double second)
{
  if(first - second < EPSILON)
    return true;
  else
    return false;
}

bool isEqual(double first, double second)
{
  if(fabs(first - second) < EPSILON)
    return true;
  else
    return false;
}

bool isEqual(double first, double second, double epsilon)
{
  if(fabs(first - second) < epsilon)
    return true;
  else
    return false;
}

bool isGreaterThan(double first, double second)
{
  if(first - second > EPSILON)
    return true;
  else
    return false;
}

bool isGreaterOrEqual(double first, double second)
{
  if(first - second > -EPSILON)
    return true;
  else
    return false;
}

bool isOdd(int N)
{
  return (N % 2 != 0);
}

bool isEven(int N)
{
  return (N % 2 == 0);
}

void solveTDM(const std::vector<double>& a1, const std::vector<double>& a2,
              std::vector<double>& a3, std::vector<double>& b,
              std::vector<double>& x)
{
  std::size_t N = b.size();

  a3[0] /= a2[0];
  b[0] /= a2[0];

  for (int i = 1; i < N; i++)
  {
    a3[i] /= a2[i] - a1[i]*a3[i-1];
    b[i] = (b[i] - a1[i]*b[i-1]) / (a2[i] - a1[i]*a3[i-1]);
  }

  x[N-1] = b[N-1];
  for (int i = N-2; i >= 0 && i < N; i--)
  {
      x[i] = b[i] - a3[i]*x[i+1];
  }
}

std::istream& safeGetline(std::istream& is, std::string& t)
{
    t.clear();

    // The characters in the stream are read one-by-one using a std::streambuf.
    // That is faster than reading them one-by-one using the std::istream.
    // Code that uses streambuf this way must be guarded by a sentry object.
    // The sentry object performs various tasks,
    // such as thread synchronization and updating the stream state.

    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();

    for(;;) {
        int c = sb->sbumpc();
        switch (c) {
        case '\n':
            return is;
        case '\r':
            if(sb->sgetc() == '\n')
                sb->sbumpc();
            return is;
        case EOF:
            // Also handle the case when the last line has no line ending
            if(t.empty())
                is.setstate(std::ios::eofbit);
            return is;
        default:
            t += (char)c;
        }
    }
}

#endif



