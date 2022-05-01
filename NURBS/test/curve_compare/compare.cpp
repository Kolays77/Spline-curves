#include "../../nurbs.h"
#include "../../num_nurbs.h"

#include "../../../include/plot.h"

#define N_plot 10000

template<typename T>
void test(int p, int n) {
    std::vector<Point<T>> cv = generate_points_sorted<T>(n);
    std::vector<T> weights = linspace<T>(T(1.0), T(5.0), cv.size());
    std::vector<T> knots = create_knots<T>(n, p);
    std::vector<Point<T>> points1 = create_curve(p, knots, weights, cv, N_plot);
    
    NURBS<T> nurbs_curve(p, weights, cv);
    std::vector<Point<T>> points2 = nurbs_curve.get_points(N_plot);
    
    plot_curve_(points1, "Analytic NURBS");
    plot_curve(cv, points2, "num_nurbs.png", "Curve comparison", "Numerical NURBS");
    
    T sum_error = compare_points(points1, points2);
    std::vector<T> errors = compare_points_vector(points1, points2);
    plot_errors(errors, "errors.png");
    std::cout << "Error : " <<  sum_error << "\n";
    PLOT_END();
}

int main(int argc, char* argv[]) {
    test<long double>(std::atoi(argv[1]), std::atoi(argv[2]));
    return 0;
}
