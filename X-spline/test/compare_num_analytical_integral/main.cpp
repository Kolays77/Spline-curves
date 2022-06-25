#include "../../xspline.h"
#include "../../../include/plot.h"


template<typename T>
void test(int n) {
    std::vector<T> cv_x = generate_vector_sorted<T>(n, 0.0, 5.0);
    std::vector<T> cv_y = generate_vector<T>(n, -5.0, 5.0);
    std::vector<T> s = generate_vector<T>(n, -1.0, 1.0);
    
    // std::vector<long double> cv_x({0.0, 2.0, 3.0, 3.5, 4.0, 5.0});
    // std::vector<long double> cv_y({5.0, -2.0, 4.0, 6.3, 7.0, -2.0});
    // std::vector<long double> s({1.0, 0.1, -0.435, 0.351, 0.5, 1.0});    

    Xspline<T> spline(cv_x, cv_y, s);
    
    auto res_x_y_1 = spline.get_points(1000);
    auto res_x_y_2 = make_curve(cv_x, cv_y, s, 1000);
    
    // make assert size_1 != size_2;
    int size = res_x_y_1.first.size();
    
    std::vector<T> errors(size);
    
    T max = 0;
    T sum = 0;
    for (int i = 0; i < size; ++i) {    
        T x1 = res_x_y_1.first[i];
        T y1 = res_x_y_1.second[i];
        
        T x2 = res_x_y_2.first[i];
        T y2 = res_x_y_2.second[i];

        T temp_err = std::sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
        errors[i] = temp_err;    
        sum += temp_err;
        if (max < temp_err) {
            max = temp_err;
        }
    }

    std::cout << "MAX error : " << max << "\n";
    std::cout << "Sum error : " << sum << "\n";
    plot_errors(errors, "errors.png");
    PLOT_END();
}

int main(int argc, char* argv[]) {
    test<long double>(std::atoi(argv[1]));
}

