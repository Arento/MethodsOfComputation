#include "InterpolationTools.h"
#include "Poly.h"
#include <iostream>

double rectangle_method(double f(double), double a, double b, double m) {
  auto h = (b - a) / m;
  double result = 0;
  for (int i = 0; i < m; ++i) {
    result += f(a + h / 2 + h * i);
  }

  result *= h;
  return result;
}

double trapezium_method(double f(double), double a, double b, double m) {
  auto h = (b - a) / m;
  double result = f(a) + f(b);
  for (int i = 1; i < m ; ++i) {
    result += 2 * f(a + i * h);
  }

  result *= h / 2;
  return result;
}

double simpson_method(double f(double), double a, double b, double m) {
  auto h = (b - a) / (2 * m);
  double result = f(a) + f(b);
  for (int i = 1; i < 2 * m ; ++i) {
    result += (2 + 2 * (i % 2)) * f(a + i * h);
  }

  result *= h / 3;
  return result;
}

double func(double x) {

  return x*x*x;
}

int main() {
  double a, b;
  int m;
  while (true) {
    std::cout << "Enter m, a, b:\n";
    std::cin >> m >> a >> b;
    std::cout << "rect " << rectangle_method(func, a, b, m) << "\n"
              << "trap " << trapezium_method(func, a, b, m) << "\n"
              << "simp " << simpson_method(func, a, b, m) << "\n";
  }
}