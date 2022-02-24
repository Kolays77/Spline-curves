#include "../../bspline.h"

// g++ test_cont.cpp 


template<typename T>
void test_cont(int p, int n) { 
    std::cout << std::setprecision(16);
    std::vector<Point<T>> cv = generate_points_sorted<T>(n);
    std::vector<T> knots = create_knots<T>(cv.size(), p);
    Bspline<T> bspline(p, knots, cv);
    std::vector<Point<Poly<T>>> coefs = bspline.get_coefs();
    std::vector<int> ks = bspline.get_ks();
    std::cout << coefs.size() << "\n\n";
    // continuity C0 check 
    
    
    for (int i = 0; i < ks.size()-1; ++i) {
        std::cout << "==========" << "\n";
        std::cout << coefs[i][0].At(knots[ks[i+1]]) << "\n";
        std::cout << coefs[i+1][0].At(knots[ks[i+1]]) << "\n";
    
    }
}


int main() {
    test_cont<long double>(3,20);
}
