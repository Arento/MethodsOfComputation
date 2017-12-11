#include "Poly.h"
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <map>

double f(double x) {
    return 1 - std::exp(-2 * x);
}

const double PI = std::acos(-1);
const size_t M = 15;

std::vector<double> chebRoots(double a, double b, size_t n) {
    auto c = (b - a) / 2;
    std::vector<double> res(n);
    for (size_t i = 0; i < n; ++i) {
        res[i] = a + c + c * std::cos(PI * (2 * i + 1) / (2 * n));
    }

    return res;
}

Poly lagrangeInterpolation(const std::vector<double>& xs, const std::vector<double> ys, size_t degree) {
    degree += 1;
    Poly res(0);
    for (size_t first = 0; first < degree; first++) {
        Poly tmp(1);
        for (size_t second = 0; second < degree; second++) {
            if (second != first) {
                tmp = tmp * Poly(-xs[second], 1) *Poly(1 / (xs[first] - xs[second]));
            }
        }

        res = res + tmp * Poly(ys[first]);
    }

    return res;
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

Poly newtonInterpolation(const std::vector<double>& xs, const std::vector<double> ys, size_t degree) {
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

Poly reverseInterpolation(const std::vector<double>& xs, const std::vector<double> ys, size_t degree) {
    return lagrangeInterpolation(ys,xs,degree);
}

int main() {
    double a, b, x;
    size_t degree;
    std::cout << "Enter interval and degree:\n";
    std::cin >> a >> b >> degree;
    auto xs = chebRoots(a, b, std::max(M, degree + 1));
    std::vector<double> ys;
    for (auto x : xs) {
        ys.push_back(f(x));
    }

    std::cout << "x f(x)" << std::endl;
    for (size_t i = 0; i < xs.size(); i++) {
        std::cout << xs[i] << " " << ys[i] << std::endl;
    }

    std::cout << "Enter point:\n";
    std::cin >> x;
    std::sort(xs.begin(), xs.end(), [x](double a, double b) { return abs(a - x) < abs(b - x); });
    for (size_t i = 0; i < xs.size(); i++) {
        ys[i] = f(xs[i]);
    }

    std::cout << "Sorted interpolation knots:" << std::endl;
    for (int i = 0; i < degree; i++) {
        std::cout << xs[i] << std::endl;
    }

    auto lagrange = lagrangeInterpolation(xs, ys, degree);
    auto newton = newtonInterpolation(xs, ys, degree);
    auto reversed = reverseInterpolation(xs, ys, degree);
    std::cout << "F(x): " << f(x) << std::endl;
    std::cout << "Newton error:   " << abs(newton(x) - f(x)) << std::endl;
    std::cout << "Lagrange error: " << abs(lagrange(x) - f(x)) << std::endl;
    std::cout << "Reverse  error: " << abs(reversed(f(x)) - x) << std::endl;

    return 0;
}