//#include "decomp.cpp"
#include "cmath.h"
//#include "decomp.cpp"
#include "Forsythe.h"
#include <iomanip>
#include <iostream>

using namespace std;

void PrintMatrixes(int n, double* a1, double* a2)
{
  cout << "Matrix     " << "            \tVector D\n";
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      cout << a1[i * n + j] << " ";
    }
    cout << "     " << a2[i] << "                 \t\n";
  }
}

long double SigmaFunc(int n, double* d, double* d1)
{
  long double si = 0.0;
  long double sig = 0.0;
  for (int i = 0; i < n; i++)
  {
    si += d[i] * d[i];
    auto div = d[i] - d1[i];
    sig += div * div;
  }
  sig = sqrt(sig);
  si = sqrt(si);
  return sig / si;
}

void Work(int n)
{
  cout << "N= " << n << "\n";
  cout << "-----FIRST PART-----\n";
  double* C = new double[n * n];
  double* CT = new double[n * n];
  double* d = new double[n];
  double* C1 = new double[n * n];
  double* d1 = new double[n];
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      auto num = 1.0 / (double)(i + 1 + j);
      C[i * n + j] = num;
    }
  }
  for (int i = 0; i < n; i++)
  {
    d[i] = 0;
    for (int j = 0; j < n; j++)
    {
      auto num = 1.0 / (double)(i + 1 + j);
      d[i] += num;
    }
  }
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      CT[i * n + j] = C[j * n + i];
    }
  }
  for (int i = 0; i < n; i++)
  {
    d1[i] = 0;
    for (int j = 0; j < n; j++)
    {
      d1[i] += CT[i * n + j] * d[i];
    }
  }
  PrintMatrixes(n, C, d);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      C1[i * n + j] = 0;
      for (int k = 0; k < n; k++)
      {
        C1[i * n + j] += CT[i * n + k] * C[k * n + j];
      }
    }
  }
  int* pivot1 = new int[n];
  int* pivot2 = new int[n];
  double cond1 = 0.0;
  int flag;
  Decomp(n, C, &cond1, pivot1);
  Solve(n, C, d, pivot1);
  cout << "cond: " << cond1 << "\n";
  cout << "\n-----SECOND PART-----\n";
  flag = 0;
  double cond2 = 0.0;
  PrintMatrixes(n, C1, d1);
  Decomp(n, C1, &cond2, pivot2);
  Solve(n, C1, d1, pivot2);
  cout << "cond: " << cond2 << "\n\n";
  auto result = SigmaFunc(n, d, d1);
  cout << "SIGMA= " << result << "\n____________________________________________________________________________________\n";
  delete[] pivot1;
  delete[] pivot2;
  delete[] C;
  delete[] CT;
  delete[] C1;
  delete[] d;  
  delete[] d1;
}

int main()
{
  cout << fixed << setprecision(2);
  int n = 4;
  while (n != 14)
  {
    Work(n);
    n += 2;
  }
}

