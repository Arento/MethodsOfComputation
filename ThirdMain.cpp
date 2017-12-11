#include "Poly.h"
#include "InterpolationTools.h"
#include <iostream>
#include <functional>

void solve(const std::vector<double>& xs, const std::vector<double>& ys, double h, uint32_t n, double x, const std::function<double(double)>& f) {
    auto nf = newtonForward(xs, ys, h, n)(x);
    auto nb = newtonBackward(xs, ys, h, n)(x);
    auto ga = gauss(xs, ys, h, n)(x);
    auto real = f(x);
    std::cout << "\n"
              << "real  " << real << "\n"
              << "forw  " << nf << "\n"
              << "back  " << nb << "\n"
              << "gauss " << ga << "\n";
    std::cout << "Forward Newton  absolute error: " << std::fabs(real - nf) << "\n";
    std::cout << "Backward Newton absolute error: " << std::fabs(real - nb) << "\n";
    std::cout << "Gauss           absolute error: " << std::fabs(real - ga) << "\n";
}

int main() {
    std::cout << "Interpolation f = sin(x) + x^2\n"
              << "Enter m, a, b:\n";

    int M;
    double A;
    double B;
    double h;
    std::cin >> M >> A >> B;
    h = (B - A) / M;
    std::function<double(double)> f = [](double x) {return std::sin(x) + x*x; };
    std::vector<double> xs;
    std::vector<double> ys;
    std::cout.precision(10);
    std::cout << "x f(x)\n";
    for (int i = 0; i <= M; i++) {
        xs.push_back(A + h*i);
        ys.push_back(f(xs.back()));
        std::cout << xs.back() << " " << ys.back() << "\n";
    }

    while (true) {
        int n;
        double x;
        std::cout << "Enter degree and point: \n";
        std::cin >> n >> x;
        if (n > M || n < 0) {
            std::cout << "Degree must me in [0,m]\n";
            continue;
        }

        solve(xs, ys, h, n, x, f);
    }
}