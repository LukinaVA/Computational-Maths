#include <math.h>
#include <stdio.h>
#include <iostream>

#include "cmath.h"

const int n = 11;

double f(double x) //F(x)
{
  return (1 - exp(-x));
}

double s(double xx) //SPLINE
{
  double x[n], y[n];
  double b[n], c[n], d[n];
  int    spFlag;

  for (int i = 0; i < n; ++i)
  {
    x[i] = 0.3 * i;
    y[i] = f(x[i]);
  }

  spline (n, 0, 0, 0, 0, x, y, b, c, d, &spFlag);

  int last = 0;
  double ff = seval (n, xx, x, y, b, c, d, &last);

  return ff;
}

double l(double xx) //LAGRANGE
{
  double x[n], y[n];
  for (int i = 0; i < n; ++i)
  {
    x[i] = 0.3 * i;
    y[i] = f(x[i]);
  }

  double ll = 0;
  double numerator, denominator;
  for (int i = 0; i < n; ++i)
  {
    numerator = 1;
    denominator = 1;
    for (int j = 0; j < n; ++j)
    {
      if (i != j)
      {
        numerator *= (xx - x[j]);
        denominator *= (x[i] - x[j]);
      }
    }
    ll += ((numerator / denominator) * y[i]);
  }

  return ll;
}

int main()
{
  double left, right, Ea, Er, result, err, posn;
  int    nfe, flag;

  left = 0.0;
  right = 3.0;
  Er = 1.0e-10;
  Ea = 0.0;

  std::cout.precision(8);
  /////// F(x)
  quanc8 (f, left, right, Ea, Er, &result, &err,
          &nfe, &posn, &flag);
  std::cout << "First integral: F(x) = " << result << std::endl;

  /////// SPLINE
  quanc8 (s, left, right, Ea, Er, &result, &err,
          &nfe, &posn, &flag);
  std::cout << "Second integral: S(x) = " << result << std::endl;

  /////// LAGRANGE
  quanc8 (l, left, right, Ea, Er, &result, &err,
          &nfe, &posn, &flag);
  std::cout << "Third integral: L(x) = " << result << std::endl;

  return 0;
}
