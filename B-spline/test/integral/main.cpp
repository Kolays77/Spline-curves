#include "../../bspline.h"
#include "../../../include/plot.h"


template <typename T>
void test(int p, int n){
    std::cout << std::setprecision(16);
    
    std::vector<Point<T>> cv = generate_points_sorted<T>(n);
    std::vector<T> knots = create_knots<T>(cv.size(), p);
    
    Bspline<T> bspline(p, knots, cv);
    auto points = bspline.get_points(1000);    
    bspline.save_coefs();
    std::cout << "---INTEGRAL---" << "\n";    
    
    T int1 = bspline.analytic_integral();
    T int2 = bspline.integral();
    
    std::cout << "analytical integral = " << int1 << "\n";
    std::cout << "numerical integral = " << int2 << "\n";
    std::cout << std::abs(int1- int2) << "\n";
    plot_curve(cv, points);
    PLOT_END();
}

int main(int argc, char const *argv[]) {
    int p = std::atoi(argv[1]);
    int n = std::atoi(argv[2]);
    test<long double>(p, n);
}
