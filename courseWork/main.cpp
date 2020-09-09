#include <math.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>

#include "cmath.h"

double L;
double K;
int M = 1;
const double g = 9.81;

double l(double x)
{
  return (cos(x*x));
}

//делаем замену в исходных уравениях для понижения порядка и получаем:
int f(int n, double t, double z[], double dzdt[])
{
  dzdt[0] = z[1];
  dzdt[1] = -K * z[0] / M - g * (1 - cos(z[2])) + (L + z[0]) * z[3] * z[3];
  dzdt[2] = z[3];
  dzdt[3] = -g / (L + z[0]) * sin(z[2]) - 2 / (L + z[0]) * z[1] * z[3];

  return 0;
}

double ff(double k1)
{
  K = k1;
  int n = 4;
  //Начальные значения:
  double x0[] = {0, 0, 0, 4};
  //Экспериментальные значения:
  double x[] = {0, 0.303, -0.465, 0.592, -0.409, 0.164, 0.180};
  //среднеквадратичный критерий
  double rmsCriterion = 0;

  double dxdt[n];
  int fail, nfe;
  double h, t, tt;

  int flag = 1;
  int maxfe = 5000;
  double relerr = 1.0e-8;
  double abserr = 1.0e-8;

  rkfinit(n, &fail);
  if (fail == 0)
  {
    t = 0.0;
    std::cout << "   t      x\n";
    std::cout << "-----------------------\n";
    for (int i = 0; i < 7; ++i)
    {
      tt = 0.4 * i;
      rkf45 (f, n, x0, dxdt, &t, tt, &relerr, abserr,
             &h, &nfe, maxfe, &flag);
      std::cout << std::setw(4) << t << "  " << std::setw(8) << x0[0] << '\n';
      rmsCriterion += (x0[0] - x[i]) * (x0[0] - x[i]);
    }
  }
  rkfend();

  return rmsCriterion;
}

//метод половинного деления
double FMIN(double A, double B, double (*f)(double k), double EPS)
{
  double x1, x2;
  while (abs(A - B) > EPS)
  {
    x1 = (A + B) / 2 - EPS / 2.1;
    x2 = (A + B) / 2 + EPS / 2.1;

    if (f(x1) > f(x2))
    {
      A = x1;
    }
    else if (f(x1) < f(x2))
    {
      B = x2;
    }
    else
    {
      A = x1;
      B = x2;
    }
  }

  return (A + B) / 2;
}

int main()
{
  //Вычисление L
  double left, right, Ea, Er, result, err, posn;
  int    nfe, flag;

  left = 0.0;
  right = 1.0;
  Er = 1.0e-10;
  Ea = 0.0;

  quanc8(l, left, right, Ea, Er, &result, &err,
          &nfe, &posn, &flag);
  L = result / 0.90452424;
  std::cout << "Initial spring length: L = " << L << std::endl;

  //Оценка значения K
  FMIN(36, 46, ff, 1.0e-8);
  std::cout.precision(8);
  std::cout << "Spring stiffness: K = " << K << std::endl;

  //увеличиваем длину пружины на 1%
  L *= 1.01;
  FMIN(36, 46, ff, 1.0e-8);
  std::cout.precision(8);
  std::cout << "Spring stiffness: K = " << K << std::endl;

  //возвращаем исходную длину и увеличиваем массу на 1%
  L = result / 0.90452424;
  M *= 1.01;
  FMIN(36, 46, ff, 1.0e-8);
  std::cout.precision(8);
  std::cout << "Spring stiffness: K = " << K << std::endl;
}
