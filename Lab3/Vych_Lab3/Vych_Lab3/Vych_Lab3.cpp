#include <iostream>
#include <fstream>
#include "cmath.h"
#include "Forsythe.h"
#include <cmath>
#include <iomanip>
using namespace std;
void func(double x, double y[], double dydx[])
{
  dydx[0] = (-14) * y[0] + 13 * y[1] + cos(1 + x);
  dydx[1] = 20 * y[0] - 30 * y[1] + atan(1 + pow(x, 2));
}
void runge_kutta3(void(*f)(double x, double y[], double dydx[]),
  double t, double out, double Zn[], double h)
{
  double k1[2], k2[2], k3[2];
  double func[2], temp[2];
  double tn;
  for (double t = 0.0; t <= out; t += h)
  {
    tn = t;
    temp[0] = Zn[0];
    temp[1] = Zn[1];
    f(tn, temp, func);
    k1[0] = h * func[0];
    k1[1] = h * func[1];
    tn = t + h / 2;
    temp[0] = Zn[0] + k1[0] / 2;
    temp[1] = Zn[1] + k1[1] / 2;
    f(tn, temp, func);
    k2[0] = h * func[0];
    k2[1] = h * func[1];
    tn = t + ((3 * h) / 4);
    temp[0] = Zn[0] + ((3 * k2[0]) / 4);
    temp[1] = Zn[1] + ((3 * k2[1]) / 4);
    f(tn, temp, func);
    k3[0] = h * func[0];
    k3[1] = h * func[1];
    Zn[0] = Zn[0] + (2 * k1[0] + 3 * k2[0] + 4 * k3[0]) / 9;
    Zn[1] = Zn[1] + (2 * k1[1] + 3 * k2[1] + 4 * k3[1]) / 9;
  }
}
int main()
{
  ofstream fout;
  fout.open("output.txt");
  double X[2] = {2, 0.5};
  double EPS = 0.0001;
  const double h_print = 0.075;
  const double t_start = 0;
  const double t_end = 1.5;
  const int neqn = 2;
  unsigned char work[6 * (neqn * sizeof(Float)) + sizeof(struct rkf_inside)];
  rkf argument;
  argument.f = func;
  argument.neqn = neqn;
  argument.re = EPS;
  argument.ae = EPS;
  argument.work = work;
  argument.flag = 1;
  argument.Y = X;
  argument.t = t_start;
  fout << "t         \tx[0]               \tx[1]\n\n";
  cout << "-----RKF45-----\n";
  cout << "t         \tx[0]               \tx[1]                   \tflag\n\n";
  for (double h = h_print; h <= t_end; h += 0.075)
  {
    argument.tout = h;
    rkf45(&argument);
    fout << argument.t << "       \t" << X[0]
      << "          \t" << X[1] << "\n";

    cout << argument.t << "       \t" << X[0] 
      << "          \t" << X[1] 
      << "             \t" << argument.flag << endl;
  }
  cout << "\n-----RUNGE-KUTTA-3-----\n";
  double h_1 = h_print;
  fout << "\nt         \tx[0]         \tx[1]\n\n";
  cout << "\n1) h= " << h_1 << "\n\n";
  cout << "t         \tx[0]         \tx[1]\n\n";
  for (double tMid = h_1; tMid <= t_end; tMid += h_1)
  {
    X[0] = 2;
    X[1] = 0.5;
    runge_kutta3(func, t_start, tMid, X, h_1);
    fout << tMid << "     \t" << X[0] << "       \t" << X[1] << "\n";
    cout << tMid << "     \t" << X[0] << "       \t" << X[1] << endl;
  }
  double h_2 = h_print / 10.0;
  int k = 0;
  fout << "\nt             \tx[0]             \tx[1]\n\n";
  cout << "\n2) h= " << h_2 << "\n\n";
  cout << "t             \tx[0]             \tx[1]\n\n";
  for (double tMid = h_2; tMid < t_end + h_2; tMid += h_2)
  {
    k++;
    X[0] = 2;
    X[1] = 0.5;
    runge_kutta3(func, t_start, tMid, X, h_2);
    if (k == 10)
    {
      fout << tMid << "     \t" << X[0] << "       \t" << X[1] << "\n";
      cout << tMid << "         \t" << X[0] << "         \t" << X[1] << endl;
      k=0;
    }
   }
  fout.close();
}
