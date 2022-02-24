#include "../include/tools.h"
#include "../include/plot.h"

// compile  g++ test.cpp -I/usr/include/python3.9  -lpython3.9 

int main() {
    std::vector<Point<double>> vec1 = generate_points_sorted<double>(20, 3.14, 6.99);
    plot_curve<double>(vec1, "test_pic_A.png");
    
    std::vector<Point<double>> vec2 = generate_points_sorted<double>(20);
    plot_curve<double>(vec2, "test_pic_B.png");
    PLOT_END();

    return 0;
}