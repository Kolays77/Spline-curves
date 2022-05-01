#include "../../nurbs.h"
#include "../../nurbs2.h"

#include "../../../include/plot.h"

template <typename T>
void test_nurbs_integrals(int p, int n ) {
    std::vector<Point<T>> cv = generate_points_sorted<T>(n);
    save_vector(cv, "points.txt");
    std::vector<T> weights = linspace(T(0.00001), T(1.0), cv.size());
    
    std::complex<T> int0_1;
    std::complex<T> int0_2;
    
    std::complex<T> int1;
    std::complex<T> int2;
    
    NURBS<T> nurbs1(p, weights, cv);
    NURBS2<T> nurbs2(p, weights, cv);
    
    std::vector<Point<T>> points_1 = nurbs1.get_points(1000);
    std::vector<Point<T>> points_2 = nurbs2.get_points(1000);
    
    int0_1 = std::complex<T>{nurbs1.numerical_integral()};
    int0_2 = std::complex<T>{nurbs2.numerical_integral()};
    
    int1 = nurbs1.analytic_integral1(1);
    int2 = nurbs2.analytic_integral1(1);
    
    std::cout << std::setprecision(15);
    std::cout << "Num: "        << int0_1 << "\n";
    std::cout << "Analytic: "    << int1 << "\n";
    std::cout << "Error: " << std::abs(int0_1 - int1) << "\n";
    std::cout << "========================\n";
    
    std::cout << "Num: "        << int0_2 << "\n";
    std::cout << "Analytic: "    << int2 << "\n";
    std::cout << "Error: " << std::abs(int0_2 - int2) << "\n";

    plot_curve_(points_1, "NURBS1. p =" + std::to_string(p));
    plot_curve(cv, points_2, "plot.png", "Compare nurbs1 nurbs2", "NURBS2. p =" + std::to_string(p));
    
    PLOT_END(); 
}


int main(int argc, char const *argv[]) {
    int p, n;
    p = std::atoi(argv[1]);
    n = std::atoi(argv[2]);
    test_nurbs_integrals<long double>(p, n);    
    return 0;
}
