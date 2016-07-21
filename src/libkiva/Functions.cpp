/* Copyright (c) 2012-2016 Big Ladder Software. All rights reserved.
* See the LICENSE file for additional terms and conditions. */

#ifndef Functions_CPP
#define Functions_CPP

#include "Functions.hpp"

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
