#include <math.h>
#include <stdio.h>
#include <iostream>

#include "cmath.h"

int main()
{
  ////////// Исходные данные
  int n = 8; //порядок матрицы
  int mdim = 8;
  int ipvt[n], flag;
  double cond;
  //при p -> 0 матрица становится особой
  double p[5] = {1.0, 0.1, 0.01, 0.0001, 0.000001};

  double a[n][n] = {
    {p[4] + 13, 2, 8, -7, 7, 5, -7, -7},
    {7, 2, -4, 2, 3, 3, -1, -2},
    {-7, 2, 1, 3, 6, -6, -3, -4},
    {-2, -8, -6, -1, 6, 2, 1, -4},
    {0, 4, -7, 1, 22, 0, -6, -6},
    {0, -3, -6, 6, 4, 13, 0, 6},
    {-8, -6, -4, 7, -5, -5, -2, 1},
    {5, 5, -2, -2, -3, 0, -7, 14}
  };

  double aT[n][n];
  for (int i = 0; i < n; ++i)
  {
    for(int j = 0; j < n; ++j)
    {
      aT[j][i] = a[i][j];
    }
  }

  double aTa[n][n] = { 0 };
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < n; ++j)
    {
      for (int tmp = 0; tmp < n; ++tmp)
      {
        aTa[i][j] += aT[i][tmp] * a[tmp][j];
      }
    }
  }

  double b[n] = {4 * p[4] + 6, 36, -25, -57, 32, 62, -71, 70};

  double aTb[n] = { 0 };
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < n; ++j)
    {
      aTb[i] += aT[i][j] * b[j];
    }
  }

  ////////// ax1 = b
  decomp (n, mdim, *a, &cond, ipvt, &flag);
  solve (n, mdim, *a, b, ipvt);

  double x1 = 0; // ||x1||
  std::cout << "Решение системы ax1 = b: { ";
  for (int i = 0; i < n; ++i)
  {
    std::cout << b[i] << " ";
    x1 += b[i] * b[i]; //Евклидова норма
  }
  std::cout << "}\n";
  x1 = sqrt(x1);

  std::cout << "COND = " << cond << std::endl;

  ////////// aTax2 = aTb
  decomp (n, mdim, *aTa, &cond, ipvt, &flag);
  solve (n, mdim, *aTa, aTb, ipvt);

  double x1x2 = 0; // ||x1 - x2||
  std::cout << "Решение системы aTax2 = aTb: { ";
  for (int i = 0; i < n; ++i)
  {
    std::cout << aTb[i] << " ";
    x1x2 += (b[i] - aTb[i]) * (b[i] - aTb[i]);
  }
  std::cout << "}\n";
  x1x2 = sqrt(x1x2);

  std::cout << "COND = " << cond << std::endl;
  double s = x1x2 / x1;
  std::cout << "SIGMA = " << s << std::endl;

  return 0;
}
