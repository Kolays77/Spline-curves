#include "../../xspline.h"
#include "../../../include/plot.h"


template<typename T>
void test(int n) {
    std::vector<T> cv_x = generate_vector_sorted<T>(n, 0.0, 5.0);
    std::vector<T> cv_y = generate_vector<T>(n, -5.0, 5.0);
    
        // std::vector<T> s = generate_vector<T>(n, -1.0, 1.0);
    std::vector<T> s = std::vector<T>(n, 1.0);

    Xspline<T> spline(cv_x, cv_y, s);
    spline.print_segments();
    
    auto res_x_y_1 = spline.get_points();
    auto res_x_y_2 = make_curve(cv_x, cv_y, s);

    plot_curve(cv_x, cv_y, res_x_y_1.first, res_x_y_1.second, "plot.png", "X-сплайн кривая", "Аналитический тип");
    plot_curve(cv_x, cv_y, res_x_y_2.first, res_x_y_2.second, "plot_num.png", "X-сплайн кривая", "Численный тип");
    
    std::complex<T> integral_1 = spline.numerical_integral();
    std::complex<T> integral_2 = spline.analytical_integral();

    std::cout << "numerical : " <<  integral_1 << "\n";
    std::cout << "analytical : " << integral_2 << "\n";
    std::cout << "Error : " <<  std::abs(integral_1 - integral_2) << "\n";
    PLOT_END();
}

int main(int argc, char* argv[]) {
    test<long double>(std::atoi(argv[1]));
}

