#include <math.h>
#include <stdio.h>
#include <iostream>

#include "cmath.h"

int f(int n, double t, double x[], double dxdt[])
{
  dxdt[0] = -130 * x[0] + 900 * x[1] + exp(-10 * t);
  dxdt[1] = 30 * x[0] - 300 * x[1] + log(1 + 100 * t * t);
  return 0;
}

void rkf3(double t,double* zn, double* z, double h)
{
  double k[2];
  double k1[2];
	double k2[2];
	double k3[2];

	f(2, t, zn, k1);
	k1[0] = h * k1[0];
	k1[1]= h * k1[1];

	k[0] = zn[0] + k1[0] / 2;
	k[1] = zn[1] + k1[1] / 2;
	f(2, t + h / 2, k, k2);
	k2[0] = h * k2[0];
	k2[1] = h * k2[1];

	k[0] = zn[0] + 3 * k2[0] / 4;
	k[1] = zn[1] + 3 * k2[1] / 4;
	f(2, t + 3 * h / 4, k, k3);
	k3[0] = h * k3[0];
	k3[1] = h * k3[1];

	z[0] = zn[0] + (2 * k1[0] + 3 * k2[0] + 4 * k3[0]) / 9;
	z[1] = zn[1] + (2 * k1[1] + 3 * k2[1] + 4 * k3[1]) / 9;
}

int main()
{
  ///FIRST
  double dx[2];
  int nfe, fail;
  double h, t, tt;

  int n = 2;
  int flag = 1;
  int maxfe = 5000;
  double x[] = {3, -1};
  double relerr = 1.0e-4;
  double abserr = 1.0e-4;

  rkfinit(n, &fail);

  if (fail == 0)
  {
    std::cout << "RKF45" << std::endl;
    std::cout << "       t      x1       x2\n";
    std::cout << "----------------------------------\n";

    for (int i = 1; i <= 20; ++i)
    {
      tt = 0.0075 * i;
      rkf45 (f, n, x, dx, &t, tt, &relerr, abserr,
             &h, &nfe, maxfe, &flag);
      std::cout << t << "  " << x[0] << "  " << x[1] << '\n';
    }

   rkfend ();
   std::cout << '\n';
  }

  ///SECOND
  x[0] = 3; x[1] = -1;

  std::cout << "RKF3: h = 0.0075" << std::endl;
  std::cout << "       t      x1       x2\n";
  std::cout << "----------------------------------\n";

  t = 0;
  for(int i = 1; i <= 20; ++i)
  {
    rkf3(t, x, dx, 0.0075);
    std::cout << t + 0.0075 << "  " << dx[0] << "  " << dx[1] << '\n';
    x[0] = dx[0];
    x[1] = dx[1];
    t = 0.0075 * i;
  }
  std::cout << '\n';

	x[0] = 3; x[1] = -1;

	std::cout << "RKF3: h = 0.00375" << std::endl;
  std::cout << "       t      x1       x2\n";
  std::cout << "----------------------------------\n";

  t = 0;
  for(int i = 1; i <= 40; ++i)
  {
    rkf3(t, x, dx, 0.00375);
    std::cout << t + 0.00375 << "  " << dx[0] << "  " << dx[1] << '\n';
    x[0] = dx[0];
    x[1] = dx[1];
    t = 0.00375 * i;
  }

return 0;
}
