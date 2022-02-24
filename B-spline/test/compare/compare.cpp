#include "../../bspline.h"

void compare_points_bspline() {
    for (int i = 1; i < 23; ++i) {
        auto points1 = test_numerical<long double>(i);
        auto points2 = test_bspline<long double>(i);
        cout << "p=" << i << " " << compare_points(points1, points2) << "\n";
    }
}


template<typename T>
void compare_bspline() {
    auto cv = load_points<T>("points.in");
    for (int p=2; p < cv.size()-2; ++p) {
        std::vector<T> knots = create_knots<T>(cv.size(), p);
        Bspline<T> bspline(p, knots, cv);
        auto points1 = bspline.get_points(1000);
        std::vector<Point<T>> points2 = create_curve<T>(p, knots, cv, 1000);
        save_points("out/points_bspline_" + to_string(p) + ".out", points1);
        save_points("out/points_num_" + to_string(p) +".out", points2);
    }
}

int main(){
	return 0;
}
