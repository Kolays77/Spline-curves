#include "../../bspline.h"


template <typename T>
void test_len_curve(int p) {
    std::initializer_list<int> N{100, 1000, 5000, 10000, 20000, 50000, 100000, 200000, 500000};
    std::vector<Point<T>> cv = load_points<T>("points.in");
    std::vector<T> knots = create_knots<T>(cv.size(), p);
    Bspline<T> bspline(p, knots, cv);
    std::ofstream out("errors_len.out");
    T len = bspline.length_curve();
    out << std::setprecision(16);
    for(int n : N) {
        std::cout << n << "\n";
        auto points = bspline.get_points(n);
        T len_segments =  curve_length(points);
        save_points("points_bspline.out", points);
        std::cout << "curve length (sum len segments) = " <<  len_segments << "\n";
        std::cout << "curve length (with eval integral) = " << len << "\n";
        out << n << "\t" << std::abs(curve_length(points) - bspline.length_curve()) << "\n";
    }
}

int main(){
	test_len_curve<double>(3);
	return 0;
}
