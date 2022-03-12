#include "../include/integral.h"

template <typename T>
void test() {
    Poly<T> poly({1, -3});
    std::vector < std::complex<T> > roots;
    roots.push_back(6);
    roots.push_back(2);
    roots.push_back(1);
    roots.push_back(20);
    
    std::vector<std::complex<T>> res = decomp(poly, roots);
    print(res);
}

int main() {
    test<double>();   
    return 0;
}