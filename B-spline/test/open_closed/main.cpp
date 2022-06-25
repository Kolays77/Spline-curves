#include "../../bspline.h"
#include "../../num_bspline.h"
#include "../../../include/plot.h"


// Реализация открытой и замкнутой B сплайн кривой.
// Параметры
// p - степень кривой
// n - количество контрольных точек
template<typename T>
void test(int p, int n) {
    int N = 10000;
    auto cv1 = generate_points_sorted<T>(n);
    auto cv2 = cv1;

    for (int i = 0; i < p; ++i) {
        cv2.push_back(cv1[i]);    
    }

    std::vector<T> knots1 = create_knots<T>(cv1.size(), p);
    std::vector<T> knots2 = create_uniform_knot_vector<T>(cv2.size(), p);
    print(knots2);

    // Bspline<T> bspline1(p, knots1, cv1);
    // Bspline<T> bspline2(p, knots2, cv2);

    std::vector<Point<T>> points1 = create_curve<T>(p, knots1, cv1, N);
    std::vector<Point<T>> points2 = create_curve<T>(p, knots2, cv2, N);
    
    plt::xlabel("X", {{"fontsize", "large"}});
    plt::ylabel("Y", {{"fontsize", "large"}}); 
    plot_curve(cv1, points1, "plot1.png", "", "B-сплайн кривая");
    

    plt::xlabel("X", {{"fontsize", "large"}});
    plt::ylabel("Y", {{"fontsize", "large"}}); 
    plot_curve(cv2, points2, "plot2.png", "", "B-сплайн кривая");
    PLOT_END();
}

int main(int argc, char const *argv[]) {
    int p = std::atoi(argv[1]);
    int n = std::atoi(argv[2]);
    test<long double>(p, n);
    return 0;
}
