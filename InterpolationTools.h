#include "Poly.h"
#include <vector>
#include <cstdint>
#include <cmath>

using Differences = std::vector<std::vector<double>>;

const double PI = std::acos(-1);

Differences getDifferences(const std::vector<double>& ys) {
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

Poly newtonForward(const std::vector<double>& xs, const std::vector<double>& ys, double h, uint32_t degree) {
    auto dsum = getDifferences(ys);
    auto q = (Poly(0, 1) - xs[0]) / h;
    auto buf = Poly(1);
    auto result = Poly(ys[0]);
    for (int i = 1; i <= degree; ++i) {
        buf = buf*(q - i + 1) / i;
        result = result + dsum[i][0] * buf;
    }

    return result;
}

Poly newtonBackward(const std::vector<double>& xs, const std::vector<double>& ys, double h, uint32_t degree) {
    auto dsum = getDifferences(ys);
    const auto M = ys.size() - 1;
    auto buf = Poly(1);
    auto q = (Poly(0, 1) - xs[M]) / h;
    auto result = Poly(ys[M]);
    for (size_t i = M - 1; i >= M - degree; --i) {
        buf = buf*(q - i + M - 1) / (M - i);
        result = result + dsum[M - i][i] * buf;
    }

    return result;
}

Poly gauss(const std::vector<double>& xs, const std::vector<double>& ys, double h, uint32_t degree) {
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

    return result;
}
