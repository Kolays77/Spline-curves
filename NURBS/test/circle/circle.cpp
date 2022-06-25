#include "../../nurbs.h"
#include "../../nurbs2.h"
#include "../../../include/plot.h"

template <typename T>
void test_circle(){
    int p = 2;
    std::vector<Point<T>> cv = load_points<T>("circle_points.in");
    std::vector<T> weights({1.0, sqrt(2)/2, 1, sqrt(2)/2, 1, sqrt(2)/2, 1,  sqrt(2)/2, 1 });
    std::vector<T> knots({0.,0.,0.,0.5,0.5,1.,1.,1.5,1.5,2.,2.,2.});
    NURBS2<T> NURBS(p, knots, weights, cv);
    std::vector<Point<T>> points = NURBS.get_points(1000);
    
	// Создание окна	    
    double w = 500;
    double h = 500;
    plt::figure_size(w, h);	
    plt::xlabel("X", {{"fontsize", "large"}});
    plt::ylabel("Y", {{"fontsize", "large"}});
	
	plot_curve(cv,points, "circle.png", "NURBS окружность");
    save_points("circle_points.out", points);
    NURBS.save_coefs();
    
    std::cout << std::setprecision(14);
    int N = 10000; // количество точек для каждого сегмента для вычисления ошибки.
    std::vector<T> errors(N * NURBS.N_segments);

    for (int i = 0; i < NURBS.N_segments; ++i) {
        T A = knots[NURBS.ks[i]];
        T B = knots[NURBS.ks[i+1]]; 
        std::vector<T> ts = linspace<T>(0.0, 1.0, N);
        for (int j = 0; j < N; ++j) {
            T temp_x = NURBS.coefs[i].first[0].At(ts[j]);
            T temp_y = NURBS.coefs[i].first[1].At(ts[j]);
            T temp_den = NURBS.coefs[i].second.At(ts[j]);
            errors[i * N + j] = std::abs(temp_x * temp_x + temp_y*temp_y - temp_den*temp_den);
        }  

    }
    plt::xlabel("Параметр t", {{"fontsize", "large"}});
    plt::ylabel("Ошибка", {{"fontsize", "large"}});
    save_vector<T>(errors, "errors.txt");
	plot_errors(errors, "errors.png");
    PLOT_END();
}

int main() {
    test_circle<long double>();
    return 0;
}
