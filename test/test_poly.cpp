#include "../include/Poly.h"


void test_error_root_finding(int n) {
    // TEST root finding for poly p = (x-1)(x-2)...(x-n)
    // n > 1
    std::cout << std::setprecision(16);
    Poly<long double> poly({1, -1});

    for (int i = 2; i <= n; ++i) {
        poly *= Poly<long double>({1, -(long double)i});
    }
    auto res1 = solve_numerical(poly, 1);
    auto res2 = solve_numerical(poly, 2);
    long double sum_err1 = 0;
    long double sum_err2 = 0;

    for (int i = 1; i <= n; ++i) {
        sum_err1 += std::abs(res1[i-1].second - (long double)i);
        sum_err2 += std::abs(res2[i-1].second - (long double)i);

    }
    std::cout << "Total err1 : " << sum_err1 << "\n";
    std::cout << "Total err2 : " << sum_err2 << "\n\n";
}



int main(int argc, char* argv[]) {
    // test cubic solve
    double a = std::atof(argv[1]);
    double b = std::atof(argv[2]);
    double c = std::atof(argv[3]);
    double d = std::atof(argv[4]);
    
    std::cout << std::setprecision(15);

    Poly<double> poly({a,b,c,d});
    std::cout << poly << "\n";
    auto roots1 = solve_numerical(poly, 1);
    for (auto root : roots1) {
         std::cout << root.first << " " << root.second << "\n";
    }

    std::cout << "\n\n";
    auto roots2 = solve_numerical(poly, 2);
    for (auto root : roots2) {
         std::cout << root.first << " " << root.second << "\n";
    }

    std::cout << "\n\n";
    auto roots3 = solve_cubic(poly);
    for (auto root : roots3) {
         std::cout << root.first << " " << root.second << "\n";
    }

    
    return 0;
}

