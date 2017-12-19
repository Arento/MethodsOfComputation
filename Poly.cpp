#include "Poly.h"
#include <ostream>
#include <algorithm>

Poly::Poly(double zeroCoef) {
    coefs.resize(1, zeroCoef);
}

Poly::Poly(double zeroCoef, double firstCoef) {
    coefs.resize(2);
    coefs[0] = zeroCoef;
    coefs[1] = firstCoef;
    reduceZeros();
}

double Poly::operator()(double val) const {
    double x = 1;
    double sum = 0;
    for (auto coef : coefs) {
        sum += x*coef;
        x *= val;
    }

    return sum;
}

Poly Poly::operator*(const Poly& op) const {
    Poly res;
    res.coefs.resize(coefs.size() + op.coefs.size(), 0);
    for (auto firstIndex = op.coefs.size(); firstIndex > 0;) {
        --firstIndex;
        for (auto secondIndex = coefs.size(); secondIndex > 0;) {
            --secondIndex;
            res.coefs[firstIndex + secondIndex] += coefs[secondIndex] * op.coefs[firstIndex];
        }
    }

    res.reduceZeros();
    return res;
}

Poly Poly::operator+(const Poly& op) const {
    Poly res;
    res.coefs = coefs;
    res.coefs.resize(std::max(coefs.size(), op.coefs.size()));
    for (auto index = op.coefs.size(); index > 0;) {
        --index;
        res.coefs[index] += op.coefs[index];
    }

    res.reduceZeros();
    return res;
}

Poly Poly::operator-(const Poly& op) const {
    return *this + (-1)*op;
}

Poly Poly::operator*(double val) const {
    Poly res = *this;
    for (int i = 0; i < coefs.size(); i++) {
        res.coefs[i] *= val;
    }

    return res;
}

Poly Poly::operator/(double val) const {
    return *this*(1 / val);
}

Poly Poly::operator+(double val) const {
    Poly res = *this;
    res.coefs[0] += val;
    return res;
}

Poly Poly::operator-(double val) const {
    Poly res = *this;
    res.coefs[0] -= val;
    return res;
}

Poly& Poly::operator=(const Poly& op) {
    coefs = op.coefs;
    return *this;
}

std::ostream& operator<<(std::ostream& os, const Poly& poly) {
    for (auto index = poly.coefs.size(); index > 1;) {
        --index;
        os << '(' << poly.coefs[index] << ") * x^" << index << " + ";
    }

    os << poly.coefs[0];
    return os;
}

void Poly::reduceZeros() {
    auto index = coefs.size();
    while (coefs.size() > 1 && coefs.back() == 0) {
        coefs.pop_back();
    }
}

Poly operator*(double val, const Poly& poly) {
    return poly*val;
}

Poly operator/(double val, const Poly& poly) {
    return poly / val;
}

Poly operator+(double val, const Poly& poly) {
    return poly + val;
}

Poly operator-(double val, const Poly& poly) {
    return val + (-1)*poly;
}