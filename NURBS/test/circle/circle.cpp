#include "../../nurbs.h"
#include "../../../include/plot.h"

template <typename T>
void test_circle(){
    int p = 2;
    std::vector<Point<T>> cv = load_points<T>("circle_points.in");
    std::vector<T> weights({1.0, sqrt(2)/2, 1, sqrt(2)/2, 1, sqrt(2)/2, 1,  sqrt(2)/2, 1 });
    std::vector<T> knots({0.,0.,0.,0.5,0.5,1.,1.,1.5,1.5,2.,2.,2.});
    NURBS<T> NURBS(p, knots, weights, cv);
    std::vector<Point<T>> points = NURBS.get_points(1000);
    plot_curve(cv,points, "circle.png", "NURBS circle");
    save_points("circle_points.out", points);
    NURBS.save_coefs();
    
    std::cout << std::setprecision(14);
    
    std::complex<T> int1(NURBS.integral());
    std::complex<T> int2(NURBS.analytic_integral());
    std::cout << "Numerical integral :" << int1 << "\n";
    std::cout << "Analytic integral :" << int2 << "\n";
    std::cout << "Error :" << std::abs(int2 - int1) << "\n";
    PLOT_END();
}

int main() {
    test_circle<long double>();
    return 0;
}