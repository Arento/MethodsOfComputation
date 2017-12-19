#include "Poly.h"
#include <vector>
#include <cstdint>
#include <cmath>
#include <iostream>

using Differences = std::vector<std::vector<double>>;

const double PI = std::acos(-1);

Differences getDifferences(const std::vector<double> &ys) {
  Differences dsum(ys.size(), std::vector<double>(ys.size(), 0));
  for (size_t i = 0; i < ys.size(); i++) {
    dsum[0][i] = ys[i];
  }

  for (size_t i = 1; i < ys.size(); i++) {
    for (size_t k = 0; k < ys.size() - i; k++) {
      dsum[i][k] = (dsum[i - 1][k + 1] - dsum[i - 1][k]);
    }
  }

  return dsum;
}

std::vector<double> chebRoots(double a, double b, size_t n) {
  auto c = (b - a) / 2;
  std::vector<double> res(n);
  for (size_t i = 0; i < n; ++i) {
    res[i] = a + c + c * std::cos(PI * (2 * i + 1) / (2 * n));
  }

  return res;
}

Poly newtonForward_equidistant(const std::vector<double> &xs, const std::vector<double> &ys, double h, uint32_t degree) {
  auto dsum = getDifferences(ys);
  auto q = (Poly(0, 1) - xs[0]) / h;
  auto buf = Poly(1);
  auto result = Poly(ys[0]);
  for (int i = 1; i <= degree; ++i) {
    buf = buf * (q - i + 1) / i;
    result = result + dsum[i][0] * buf;
  }

  //std::cout << "\n" << result << "\n";
  return result;
}

Poly newtonBackward_equidistant(const std::vector<double> &xs, const std::vector<double> &ys, double h, uint32_t degree) {
  auto dsum = getDifferences(ys);
  const auto M = ys.size() - 1;
  auto buf = Poly(1);
  auto q = (Poly(0, 1) - xs[M]) / h;
  auto result = Poly(ys[M]);
  for (size_t i = M - 1; i >= M - degree; --i) {
    buf = buf * (q - i + M - 1) / (M - i);
    result = result + dsum[M - i][i] * buf;
  }

  //std::cout << "\n" << result << "\n";
  return result;
}

Poly gauss_equidistant(const std::vector<double> &xs, const std::vector<double> &ys, double h, uint32_t degree) {
  auto dsum = getDifferences(ys);
  const auto M = ys.size() - 1;
  auto index = M / 2 - 1 + M % 2;
  auto buf = Poly(1);
  auto result = Poly(ys[index]);
  auto q = Poly(-xs[index], 1) / h;
  for (int i = 1; i <= degree / 2; ++i) {
    buf = buf * (q + i - 1) / (2 * i - 1);
    result = result + dsum[2 * i - 1][index] * buf;
    buf = buf * (q - i) / (2 * i);
    result = result + dsum[2 * i][index - 1] * buf;
    --index;
  }

  //std::cout << "\n" << result << "\n";
  return result;
}

std::vector<double> separatedDifferences(const std::vector<double>& xs, const std::vector<double>& ys) {
  std::vector<std::vector<double>> dsum(xs.size(), std::vector<double>(xs.size(), 0));
  for (size_t i = 0; i < xs.size(); i++) {
    dsum[0][i] = ys[i];
  }

  for (size_t i = 1; i < xs.size(); i++) {
    for (size_t k = 0; k < xs.size() - i; k++) {
      dsum[i][k] = (dsum[i - 1][k + 1] - dsum[i - 1][k]) / (xs[k + i] - xs[k]);
    }
  }

  std::vector<double> res;
  for (size_t i = 0; i < xs.size(); i++) {
    res.push_back(dsum[i][0]);
  }

  return res;
}

Poly lagrange(const std::vector<double> &xs, const std::vector<double> &ys, size_t degree) {
  degree += 1;
  Poly res(0);
  for (size_t first = 0; first < degree; first++) {
    Poly tmp(1);
    for (size_t second = 0; second < degree; second++) {
      if (second != first) {
        tmp = tmp * Poly(-xs[second], 1) / (xs[first] - xs[second]);
      }
    }

    res = res + tmp * ys[first];
  }

  std::cout << std::endl << res << std::endl;
  return res;
}

Poly newton(const std::vector<double>& xs, const std::vector<double> ys, size_t degree) {
  auto sepdiffs = separatedDifferences(xs, ys);
  Poly res(0);
  for (size_t i = 0; i < degree + 1; i++) {
    Poly buf(sepdiffs[i]);
    for (size_t j = 0; j < i; j++) {
      buf = buf * Poly(-xs[j], 1);
    }

    res = res + buf;
  }
  return res;
}
