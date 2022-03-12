#include "../../nurbs.h"

// Integral of the NURBS curve with a fixed number of control points for powers of 2 <= p <= 10.
// The average, minimum and maximum error for the analytical and numerical integral is calculated.
// Number of tests = N_test

#define N_test 100
#define N_cv 12

template <typename T>
void test_nurbs_integrals(int type) {
    std::cout << std::setprecision(15);

    std::vector<T> weights, knots;
    std::vector<T> errors_min(9, 1e+50);
    std::vector<T> errors_max(9, -1e+50);
    std::vector<T> errors_sum(9, 0.0);
    std::vector<Point<T>> cv;
    std::complex<T> int1, int2;
    T error;

    for(int p = 2; p <= 10; ++p) {
        weights = linspace(T(1.0), T(10.0), N_cv);
        knots = create_knots<T>(N_cv, p);
        for (int i = 1; i <= N_test; ++i) {
            cv = generate_points_sorted<T>(N_cv);
            NURBS<T> NURBS(p, knots, weights, cv);
            int1 = NURBS.analytic_integral(type);
            int2 = std::complex<T>{NURBS.integral()};
            error = std::abs(int1 - int2);
            errors_sum[p-2] += error;
            errors_max[p-2] = error > errors_max[p-2] ? error : errors_max[p-2];
            errors_min[p-2] = error < errors_min[p-2] ? error : errors_min[p-2];
        }
    }

    for (auto& v : errors_sum) {
        v /= N_test;
    }
    std::cout << "Average" << "\n";
    for (int p = 2; p <= 10; ++p) {
        std::cout << p << " : " << errors_sum[p-2] << "\n"; 
    }

    std::cout << "MIN" << "\n";
    for (int p = 2; p <= 10; ++p) {
        std::cout << p << " : " << errors_min[p-2] << "\n"; 
    }

    std::cout << "MAX" << "\n";
    for (int p = 2; p <= 10; ++p) {
        std::cout << p << " : " << errors_max[p-2] << "\n"; 
    }
    std::cout << "-------------\n\n";
}


int main() {
    test_nurbs_integrals<long double>(1);
    test_nurbs_integrals<long double>(2);
}