#include "../../nurbs.h"
#include "../../../include/plot.h"

#define N_points 100

template <typename T>
void integral_fixed_deg(int p) {
    // integral with increasing number of points
    // для нового числа точек происходит перегенерация
    std::vector<T> errors; 
    std::vector<int> n_points;
    for(int N = p+2; N < N_points; ++N) {
        
        std::vector<Point<T>> cv = generate_points_sorted<T>(N);
        std::vector<T> weights = linspace(T(1.0), T(10.0), cv.size());
        std::vector<T> knots = create_knots<T>(cv.size(), p);
        NURBS<T> NURBS(p, knots, weights, cv);
        
        std::complex<T> int1 = NURBS.analytic_integral(1);
        std::complex<T> int2(NURBS.integral());
        T error = std::abs(int1 - int2);
        n_points.push_back(N);
        errors.push_back(error);
    } 
    plot_errors(n_points, errors, "errors_fixed_deg_" + std::to_string(p) + ".png");
}


int main() {
	for (int p = 2; p <= 7; ++p) {
		integral_fixed_deg<long double>(p);
	}
    PLOT_END();
    return 0;    
}
