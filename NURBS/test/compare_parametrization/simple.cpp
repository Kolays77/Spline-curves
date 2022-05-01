#include "../../nurbs.h"
#include "../../../include/plot.h"


template <typename T>
void plot_denominators(NURBS<T>& nurbs){
    // plot D(x) where x==t
    int n_ks = nurbs.ks.size();
    for (int i = 0; i < n_ks; ++i) {
        T from = nurbs.knots[nurbs.ks[i]];
        T to = nurbs.knots[nurbs.ks[i] + 1];
        std::vector<T> xs = linspace<T>(from, to, 100);
        std::vector<T> ys = nurbs.coefs[i].second.At(xs);
        if (i == n_ks-1) {
            plot_curve(xs, ys, "plot_den.png", "Denominators");
        } else {
            plot_curve_(xs, ys);
        }
    }
}

template <typename T>
void test_nurbs_uniform(int p, int n){
    // UNIFORM PARAMETRIZATION
    std::vector<Point<T>> cv = generate_points_sorted<T>(n);

    std::vector<T> knots = create_uniform_knot_vector<T>(n, p);
    //std::vector<T> knots = generate_vector_sorted<T>(n, p); // не работет
    //std::vector<T> knots = generate_clamped_random_knot_vector<T>(n, p);
    
    //std::vector<T> weights = linspace(T(1.0), T(2.0), cv.size());
    std::vector<T> weights = generate_vector(cv.size(), T(0.001), T(1.0));
    
    NURBS<T> NURBS(p, knots, weights, cv);
    
    std::vector<Point<T>> points = NURBS.get_points(1000);
    plot_curve<T>(cv, points, "nurbs_uniform.png"); 
    NURBS.save_coefs();
    
    plot_denominators<T>(NURBS);
    PLOT_END();

    int knots_size = NURBS.knots.size();
    std::cout << "KNOT VECTOR : [";
    for (int i = 0; i < knots_size; ++i) {
        std::cout << NURBS.knots[i] << "\n";
    }

    T int0 = NURBS.numerical_integral();
    std::complex<T> int1 = NURBS.analytic_integral1(1);
    std::complex<T> int2 = NURBS.analytic_integral1(2);
    T error1 = std::abs(int1 - std::complex<T>(int0));
    T error2 = std::abs(int2 - std::complex<T>(int0));

    std::cout << "---------------INTEGRAL WITH UNIFORM PARAMETRIZATION------------\n";
    std::cout << "numerical : " << int0 << "\n";
    std::cout << "Analytical(eigenvalues) :"     << "  " << int1 << " Error : " << error1 << "\n";
    std::cout << "Analytical(laguerre) :"  << "  " << int2 << " Error : " << error2 << "\n";
    
}

int main(int argc, char* argv[]) {
    test_nurbs_uniform<long double>(std::atoi(argv[1]), 
                                    std::atoi(argv[2]));
    return 0;
}
