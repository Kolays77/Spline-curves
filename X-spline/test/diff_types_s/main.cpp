#include "../../xspline.h"
#include "../../../include/plot.h"

template<typename T>
void test() {
    // std::vector<T> cv_x = generate_vector_sorted<T>(n, 0.0, 5.0);
    // std::vector<T> cv_y = generate_vector<T>(n, -5.0, 5.0);
    // std::vector<T> s = generate_vector<T>(n, -1.0, 1.0);

    std::vector<T> cv_x({1.0, 2.0, 3.0, 4.0, 5.0});
    std::vector<T> cv_y({2.0, 4.0, 2.0, 4.0, 2.0});

    std::vector<T> s1({ 1.0,  1.0,  1.0,  1.0,  1.0});  
    std::vector<T> s2({ 1.0,  0.0,  1.0,  -1.0,  0.0});  
    std::vector<T> s3({-1.0, -1.0, -1.0, -1.0, -1.0});  

    Xspline<T> spline1(cv_x, cv_y, s1);
    Xspline<T> spline2(cv_x, cv_y, s2);
    Xspline<T> spline3(cv_x, cv_y, s3);
    
    double w = 900;
    double h = 300;

    plt::figure_size(w, h);
    spline1.plot_segments("1_segm.png");
    plt::figure_size(w, h);
    spline2.plot_segments("2_segm.png");
    plt::figure_size(w, h);
    spline3.plot_segments("3_segm.png");
    
    
    // auto res_x_y_1 = spline1.get_points();
    // auto res_x_y_2 = spline2.get_points();
    // auto res_x_y_3 = spline3.get_points();
    
    
    // plot_curve(cv_x, cv_y, res_x_y_1.first, res_x_y_1.second, "1.png");
    // plot_curve(cv_x, cv_y, res_x_y_2.first, res_x_y_2.second, "2.png");
    // plot_curve(cv_x, cv_y, res_x_y_3.first, res_x_y_3.second, "3.png");
    PLOT_END();
}

int main(int argc, char* argv[]) {
    test<long double>();
}
