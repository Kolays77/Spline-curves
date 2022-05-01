#include "../../nurbs.h"
#include "../../../include/plot.h"

#define N 100

template <typename T>
void test(int p, int n) {
    T h = 0.1;
    std::vector<Point<T>> cv = generate_points_sorted<T>(n);
    std::vector<T> weights = linspace(T(1.0), T(2.0), cv.size());    
    
    NURBS<T> NURBS(p, weights, cv);
    
    T t0 = 0.0;
    std::vector<T> hs = linspace<T>(h, 1.0, N);

    std::vector<T> vec_error1(hs.size());
    std::vector<T> vec_error2(hs.size());
    
    for (int i = 0; i < N; ++i) {
        T int0 = NURBS.numerical_integral(t0, hs[i]);
        std::complex<T> int1 = NURBS.analytic_integral1(1, t0, hs[i]);
        std::complex<T> int2 = NURBS.analytic_integral1(2, t0, hs[i]);
        
        T error1 = std::abs(int1 - std::complex<T>(int0));
        T error2 = std::abs(int2 - std::complex<T>(int0));        
        
        vec_error1[i] = error1;
        vec_error2[i] = error2;
    }    
    plot_errors(hs, vec_error1, "errors1_" + std::to_string(p) + "_" + std::to_string(n) + ".png");
    plot_errors(hs, vec_error2, "errors2_" + std::to_string(p) + "_" + std::to_string(n) + ".png");
    
    PLOT_END();
}


int main(int argc, char* argv[]) {
    test<long double>(std::atoi(argv[1]), std::atoi(argv[2]));
    return 0;    
}
