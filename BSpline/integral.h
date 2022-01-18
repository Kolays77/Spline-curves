#pragma once
#include "Poly.h"

// TODO add cond (from - to) < EPS => 0

long double EPS = 1e-12;

template<typename T>
using Rational = std::pair<Point<Poly<T>>, Poly<T>>;


template<typename T>
T Abs(T value) {
    return value < 0 ? -value : value;
}

template<typename T>
std::complex<T> Integral_first_type(std::complex<T> NUM,
                                    std::complex<T> a, int n,  T from,  T to) {
    //if (std::abs(den.eval(to)) < EPS || std::abs( den.eval(from)) < EPS)
    //    throw std::logic_error("Integral does not converge");
    std::complex<T> A(from);
    std::complex<T> B(to);

    if(n == 1)
        return NUM * ( std::log(B - a) - std::log(A - a)) ;
    else {

        return  NUM / std::complex<T>{n-1.} * (std::complex<T>{1.0}/std::pow(A-a, n-1.) - std::complex<T>{1.0}/std::pow(B-a, n-1.));
    }


}


// A / (x + a)
    template<typename T>
T Integral(const Poly<std::complex<T>>& num,
           const Poly<std::complex<T>>& den, int n,  T from,  T to) {
    // integral elementary fractions :
    // (1.1) A/x+a
    // (1.2) A/(x+a)^n

    if (from == to) return 0;

    if (num.deg() == 0 && den.deg() == 1)
        return Integral_first_type(num, den, n, from, to);

    return 0;
}

template<typename T>
T numerical_integral(const Poly<T>& poly, T A, T B) {
    // Legendre's quadrature. Need integrate_points.in, integrate_weights.in
    // https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_quadrature
    auto points = load_vector<T>("integrate_points.in");
    auto weights = load_vector<T>("integrate_weights.in");
    int n = points.size();
    T sum = 0.0;
    for (int i = 0; i < n; ++i) {
        sum += weights[i]*poly.At((B-A)*points[i]/2 + (A+B)/2);
    }
    return (B-A)/2.0 * sum;
}


template<typename T>
T numerical_integral(const Poly<T>& num,
                     const Poly<T>& den, T A, T B) {
    // Legendre's quadrature. Need integrate_points.in, integrate_weights.in
    // https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_quadrature
    auto points  = load_vector<T>("integrate_points.txt");
    auto weights = load_vector<T>("integrate_weights.txt");
    int n = points.size();
    T sum = 0.0;
    T x = 0.0;
    for (int i = 0; i < n; ++i) {
        x = (B-A)*points[i]/2 + (A+B)/2;
        sum += weights[i]* num.At(x) / den.At(x);
    }
    return (B-A)/2.0 * sum;
}

