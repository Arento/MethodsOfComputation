#include "Poly.h"
#include "InterpolationTools.h"
#include <iostream>
#include <iomanip>
#include <functional>


enum OperationType {
  MONOTONIC_REVERSE_INTERPOLATION,
  NON_MONOTONIC_REVERSE_INTERPOLATION,
  NUMERIC_DIFFERENTIATION
};

std::vector<double> bisection(const Poly &f, double eps, double start, double end) {
  std::vector<std::pair<double, double> > separatedRoots;
  const int MAX_INTERVALS = 1000000;
  double step = end - start / double(MAX_INTERVALS);
  for (int i = 0; i < MAX_INTERVALS; i++) {
    if (f(start) * f(start + step) < 0) {
      separatedRoots.emplace_back(std::pair<double, double>{start, start + step});
      start += step;
    }
  }

  auto solveEq = [](double start, double end, double eps, const Poly &f) {
    while (end - start > 2 * eps) {
      auto middle = (start + end) / 2;
      if (f(start) * f(middle) < 0) {
        end = middle;
      } else {
        start = middle;
      }
    }

    return (end + start) / 2;
  };

  std::vector<double> solutions;
  for (auto pair : separatedRoots) {
    solutions.push_back(solveEq(pair.first, pair.second, eps, f));
  }

  return solutions;
}

void
solve(OperationType type, const std::vector<double> &xs, const std::vector<double> &ys, double y,
      const std::function<double(double)> &f,
      const std::function<double(double)> &f_diff,
      const std::function<double(double)> &f_diff2) {
  switch (type) {
    case MONOTONIC_REVERSE_INTERPOLATION: {
      auto x = newton(ys, xs, ys.size() - 1)(y);
      std::cout << "Absolute error for x = " << x << ": " << fabs(f(x) - y) << std::endl;
      break;
    }

    case NON_MONOTONIC_REVERSE_INTERPOLATION: {
      auto interpolated = lagrange(xs, ys, ys.size() - 1);
      auto solutions = bisection(interpolated - y, 1e-8, xs.front(), xs.back());

      for (auto x : solutions) {
        std::cout << "Absolute error for x = " << x << ": " << fabs(f(x) - y) << std::endl;
      }
      break;
    }

    case NUMERIC_DIFFERENTIATION: {
      std::cout << "x f'(x) f'(x)_error f''(x) f''(x)_error\n";
      double f_dx = (-3 * ys[0] + 4 * ys[1] - ys[2]) / (xs[2] - xs[0]);
      std::cout << xs[0] << " " << f_dx << " " << fabs(f_dx - f_diff(xs[0])) << " - -\n";
      std::cout.flush();
      for (int i = 1; i + 1 < xs.size(); ++i) {
        f_dx = (ys[i + 1] - ys[i - 1]) / (xs[i + 1] - xs[i - 1]);
        double f_dx2 = (ys[i + 1] - 2 * ys[i] + ys[i - 1]) / ((xs[i] - xs[i - 1]) * (xs[i] - xs[i - 1]));
        std::cout << xs[i] << " "
                  << f_dx << " "
                  << fabs(f_dx - f_diff(xs[i])) << " "
                  << f_dx2 << " "
                  << fabs(f_dx2 - f_diff2(xs[i])) << "\n";
        std::cout.flush();
      }

      f_dx = (3 * ys[ys.size() - 1] - 4 * ys[ys.size() - 2] + ys[ys.size() - 3]) / (xs[2] - xs[0]);
      std::cout << xs.back() << " " << f_dx << " " << fabs(f_dx - f_diff(xs.back())) << " - -\n";
      std::cout.flush();
      break;
    }
  }

  std::cout << std::endl;
}

int main() {
  std::cout << "Reverse interpolation f = sin(x) + x^2\n"
            << "Enter m, a, b:\n";

  int M;
  double A;
  double B;
  double h;
  std::cin >> M >> A >> B;
  h = (B - A) / M;
  //std::function<double(double)> f = [](double x) { return std::sin(x) + x * x; };
  //std::function<double(double)> f_diff = [](double x) { return std::cos(x) + 2 * x; };
  //std::function<double(double)> f_diff2 = [](double x) { return -std::sin(x) + 2; };
  std::function<double(double)> f = [](double x) { return x * x * x; };
  std::function<double(double)> f_diff = [](double x) { return 3 * x * x; };
  std::function<double(double)> f_diff2 = [](double x) { return 6 * x; };


  std::vector<double> xs;
  std::vector<double> ys;
  std::cout << std::fixed;
  std::cout << std::setprecision(10);
  std::cout << "x f(x)\n";
  for (int i = 0; i <= M; i++) {
    xs.push_back(A + h * i);
    ys.push_back(f(xs.back()));
    std::cout << xs.back() << " " << ys.back() << "\n";
  }

  std::cout << "\nOperation types:\n"
            << "1 - monotonic function\n"
            << "2 - non monotonic function\n"
            << "3 - numeric differentiation\n"
            << "4 - exit\n\n";
  while (true) {
    double y;
    int op;
    std::cout << "Enter operation type and value: \n";
    std::cin >> op >> y;
    if (op == 4) {
      break;
    }

    solve(OperationType(op - 1), xs, ys, y, f, f_diff, f_diff2);
  }

  return 0;
}