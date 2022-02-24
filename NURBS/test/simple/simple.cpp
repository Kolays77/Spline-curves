#include "../../nurbs.h"
#include "../../../include/plot.h"

//compile  g++ simple.cpp -I/usr/include/python3.9  -lpython3.9 -o simple

template <typename T>
void test_nurbs_uniform(int p, int n){
      
    std::vector<Point<T>> cv = generate_points_sorted<T>(n);
    
    std::vector<T> weights = linspace(T(1.0), T(5.0), cv.size());
    NURBS<T> NURBS(p, weights, cv);
    std::vector<Point<T>> points = NURBS.get_points(1000);
    plot_curve<T>(cv, points, "nurbs.png"); 
    save_points("points_nurbs.out", points);

    NURBS.save_coefs();
    
    T int0 = NURBS.integral();
    std::complex<T> int1 = NURBS.analytic_integral(1);
    std::complex<T> int2 = NURBS.analytic_integral(2);

    T error1 = std::abs(int1 - std::complex<T>(int0));
    T error2 = std::abs(int2 - std::complex<T>(int0));

    std::cout << "---------------INTEGRAL------------\n";
    std::cout << "numerical : " << int0 << "\n";
    std::cout << "Analytical(laguerre) :"     << "  " << int1 << " Error : " << error1 << "\n";
    std::cout << "Analytical(eigenvalues) :"  << "  " << int2 << " Error : " << error2 << "\n";
    PLOT_END();
}

int main(int argc, char* argv[]) {
    test_nurbs_uniform<long double>(std::atoi(argv[1]), 
                                    std::atoi(argv[2]));
    return 0;
}
