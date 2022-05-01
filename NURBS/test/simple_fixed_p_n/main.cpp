#include "../../nurbs.h"
#include "../../../include/plot.h"

template <typename T>
void test_nurbs_integrals(int p, int n ) {

    std::vector<Point<T>> cv = generate_points_sorted<T>(n);
    save_vector(cv, "points.txt");

    std::vector<T> weights;
    weights = linspace(T(1.0), T(2.0), cv.size());
    std::complex<T> int0, int1, int2;
    NURBS<T> NURBS(p, weights, cv);
    NURBS.save_coefs();
    std::vector<Point<T>> points = NURBS.get_points(1000);

    int0 = std::complex<T>{NURBS.numerical_integral()};
    int1 = NURBS.analytic_integral1(1);
    int2 = NURBS.analytic_integral1(2);
    
    T error; 
    std::cout << std::setprecision(15);
    std::cout << "Num : "        << int0 << "\n";
    std::cout << "Analytic 1: "  << int1 << "\n";
    std::cout << "Analytic 2 : " << int2 << "\n";
    
    plot_curve(cv, points, "NURBS.png", "NURBS curve");
    PLOT_END();
}

int main(int argc, char const *argv[]) {
    int p, n;
    p = std::atoi(argv[1]);
    n = std::atoi(argv[2]);
    test_nurbs_integrals<long double>(p, n);    
    return 0;
}
