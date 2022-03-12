#include "../include/integral.h"


void test1() {
    // wolfram : integral (x^2+12x + 6) (x^2-5x+4) / x^2 / (x-1)^2 / (x-2)^2 / (x-3)^2 / (x-4)^2 / (x) from 5 to 6
    std::complex<double> true_val = 0.0016424140532364779094;

    Poly<double> poly1({1, 12, 6});
    Poly<double> poly2({1, -5, 4});

    double from = 5.0;
    double to = 6.0;

    std::vector<std::complex<double>> roots({0.0, 1.0, 2.0, 3.0, 4.0});
    int k_repeat = 0;

    std::complex<double> res = integral(poly1, poly2, roots, k_repeat, from, to);
    std::cout << res << "\n";
    std::cout << std::abs(true_val - res);
}

void test2() {
    // integral 4 / (x+2)^3 from 0 to 1 
    std::complex<double> true_val = 0.27777777777778;
    std::complex<double> C = 4.0;
    std::complex<double> root = -2.0;
    double from = 0.0;
    double to = 1.0;
    std::complex<double> res = integral_type_1(C, root, 3, from, to);
    std::cout << res << "\n";
    std::cout << std::abs(true_val - res) << "\n";
}


void test3() {
    // integral 314x^4 + x^3+3x^2+4x+7 from 0 to 1 
    std::complex<double> true_val = 73.05;
    Poly<double> p({314., 1.,3.,4.,7.0});
    double from = 0.0;
    double to = 1.0;
    std::complex<double> res = p.integral(from, to);
    std::cout << res << "\n";
    std::cout << std::abs(true_val - res) << "\n"; 
}

void test4() {
    std::complex<double> A = 1.0;
    
    std::vector<std::complex<double>> roots;
    roots.push_back({0, -1});
    roots.push_back({0, 1});
    
    std::complex<double> res = integral_type_7(A, roots[0], roots[1], 0.0, 1.0);
    std::cout << res << "\n";
}


// integral 1 / (x - i) / (x + i) / (x - i) from 0 to 1
// --> 0.25 + 0.642699i

void test5() {
    std::complex<double> A = 1.0;
    
    std::complex<double> a({0, 3});
    std::complex<double> b({0, -1});

    std::complex<double> res_12 = integral_type_12(A, a, b, 0.0, 1.0);
    std::complex<double> res_1_1 = integral_type_1(A, a, 1, 0.0, 1.0);
    std::complex<double> res_1_2 = integral_type_1(A, b, 1, 0.0, 1.0);
    
    std::cout << a << " " << res_1_1 << "\n";
    std::cout << b << " " << res_1_2 << "\n";
    std::cout << "res_12" << res_12 << "\n";
}



int main() {
    test5();
    return 0;
}