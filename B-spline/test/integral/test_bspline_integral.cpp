#include "../../bspline.h"
#include "../../../include/plot.h"

template <typename T>
std::vector<Point<T>> test_bspline(int p, int n){
    std::cout << std::setprecision(16);
    
    std::vector<Point<T>> cv = generate_points_sorted<T>(n);
    std::vector<T> knots = create_knots<T>(cv.size(), p);
    Bspline<T> bspline(p, knots, cv);
    
    auto points = bspline.get_points(1000);
    
    save_points("points_bspline.dat", points);
    //bspline.save_coefs();
    std::cout << "---INTEGRALS---" << "\n";    
    std::cout << "analytical integral = " << bspline.analytic_integral() << "\n";
    std::cout << "numerical integral = " << bspline.integral() << "\n";
    plot_curve(cv, points);
    PLOT_END();
    return points;
}

int main() {
    test_bspline<long double>(3, 20);
}
