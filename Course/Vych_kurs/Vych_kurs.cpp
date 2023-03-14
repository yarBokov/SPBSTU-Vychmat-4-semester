#include "Forsythe.h"
#include <iostream>
#include <cmath>
////#include "decomp.cpp"

using namespace std;
static double E = 0.0;
static double epsrel = 1.0e-6;
static double epsabs = 1.0e-18;
const double L = 1;
double A, B, F, G;

double funE(double x)
{
  return sqrt(1 + pow(x, 5));
}

void funcRKF(double x, double Z[], double dzdx[])
{
  dzdx[0] = (x / E) * (Z[1] - Z[0]);
  dzdx[1] = Z[0] + 0 * Z[1];
}

double F_forZeroin(double param)
{
  rkf base;
  base.f = funcRKF;
  const int neqn = 2;
  //написано в документации
  unsigned char work[6 * (neqn * sizeof(Float)) + sizeof(struct rkf_inside)];
  base.neqn = neqn;
  base.work = work;
  base.ae = epsabs;
  base.re = epsrel;
  base.flag = 1;
  double Z[] = { param , A };
  base.Y = Z;
  base.t = 0;
  base.tout = L;
  rkf45(&base);
  //cout << Z[0] << "      " << Z[1] << "\n";
  cout << base.flag << "\n\n";
  return Z[1] - B;
}


int main()
{
  double end1 = 0;
  double end2 = 0.5;
  double errest;
  int nofun;
  double flag;
  double integral = Quanc8(funE, end1, end2, epsabs, epsrel, &errest, &nofun, &flag);
  E = 1.994827 * integral;
  double Brr[] = {23, 32, 36, 36};
  double Arr[] = { 5, 7, 6, 5,
    7, 10, 8, 7,
    6, 8, 10, 9,
    5, 7, 9, 10 };
  int* pivot = new int[4];
  double cond;

  Decomp(4, Arr, &cond, pivot);
  Solve(4, Arr, Brr, pivot);
  
  A = Brr[0];
  B = Brr[1];
  F = Brr[2];
  G = Brr[3];
  auto T0 = Zeroin(F_forZeroin, F, G, epsabs);
  rkf base;
  base.f = funcRKF;
  const int neqn = 2;
  //написано в документации
  unsigned char work[6 * (neqn * sizeof(Float)) + sizeof(struct rkf_inside)];
  base.neqn = neqn;
  base.work = work;
  base.ae = epsabs;
  base.re = epsrel;
  base.flag = 1;
  double Z[] = { T0, 0 };
  base.Y = Z;
  base.t = 0;
  double h_print = 0.05;
  base.tout = 0;
  rkf45(&base);
  cout << base.t << "       \t" << Z[0]
    << "          \t" << Z[1]
    << "             \t" << base.flag << endl;
  for (double h = h_print; h <= L; h += h_print)
  {
    base.tout = h;
    rkf45(&base);
    cout << base.t << "       \t" << Z[0]
      << "          \t" << Z[1]
      << "             \t" << base.flag << endl;
  }
  base.tout = 1;
  rkf45(&base);
  cout << base.t << "       \t" << Z[0]
    << "          \t" << Z[1]
    << "             \t" << base.flag << endl;
}

