#include <iostream>

#include "bspline.h"
#include "num_bspline.h"
#include "tools.h"
#include "nurbs.h"
#include "Point.h"


using namespace  std;

void compare_points_bspline() {
    for (int i = 1; i < 23; ++i) {
        auto points1 = test_numerical<long double>(i);
        auto points2 = test_bspline<long double>(i);
        cout << "p=" << i << " " << compare_points(points1, points2) << "\n";
    }
}



int main(){
    cout << setprecision(15);
    test_nurbs<long double>(3);
    test_nurbs_integral<long double>("eigenvalues"); // or type "laguerre"
    return 0;
}
