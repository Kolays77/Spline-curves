#include "../../bspline.h"
#include "../../num_bspline.h"
#include "../../../include/plot.h"

// Построение графика B-сплайн кривой (аналитич. и числен. получение) и графика ошибки численной и аналитической кривой.
// Параметры
// p - степень кривой
// n - количество контрольных точек
template<typename T>
void test(int p, int n) {
    int N = 10000;
    auto cv = generate_points_sorted<T>(n);
    std::vector<T> knots = create_knots<T>(cv.size(), p);

    Bspline<T> bspline(p, knots, cv);
    std::vector<Point<T>> points1 = bspline.get_points(N);
    std::vector<Point<T>> points2 = create_curve<T>(p, knots, cv, N);
    
    plt::xlabel("X", {{"fontsize", "large"}});
    plt::ylabel("Y", {{"fontsize", "large"}});
    
    plot_curve_(cv, points1, "Аналит.");
    plot_curve(points2, "plot.png", "", "Числен.");
    
    std::vector<T> errors = compare_points_vector<T>(points1,points2); 
    plt::xlabel("Параметр t", {{"fontsize", "large"}});
    plt::ylabel("Погрешность", {{"fontsize", "large"}});
    plot_errors(errors, "errors.png");
    PLOT_END();
}

int main(int argc, char const *argv[]) {
    int p = std::atoi(argv[1]);
    int n = std::atoi(argv[2]);
    test<long double>(p, n);
    return 0;
}
