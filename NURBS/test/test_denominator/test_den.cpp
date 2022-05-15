#include "../../nurbs.h"
#include "../../nurbs2.h"
#include "../../../include/plot.h"

#define N_start 10
#define N_end 1000


// TODO уточнение корня 
template < typename T>
std::complex<T> Newton_solve(Poly<T>& fx, std::complex<T> x0) {
    T eps = 1e-14;
    Poly<T> dfx = fx.der();
    std::complex<T> x1  = x0 - fx.At(x0) / dfx.At(x0); 
    while (std::abs(x1 - x0) > eps) { 
        x0 = x1;
        x1 = x0 - fx.At(x0) / dfx.At(x0); 
    }
    return x1;
}


template <typename T>
void test(int p) {
    std::vector<T> errors1(N_end - N_start + 1, 0.0);
    std::vector<T> errors2(N_end - N_start + 1, 0.0);
    std::vector<int> NS;
    for (int n = N_start; n <= N_end; ++n) {
        std::cout << n << "\n";
        NS.push_back(n);
        
        std::vector<Point<T>> cv = generate_points_sorted<T>(n);
        //std::vector<T> weights = linspace(T(1.0), T(2.0), cv.size()); // знаменатель не зависит от вектора узлов
        std::vector<T> weights = generate_vector<T>(cv.size(), T(0.5), T(10.0)); // знаменатель не зависит от вектора узлов
    
        NURBS<T>  nurbs1(p, weights, cv);    
        NURBS2<T> nurbs2(p, weights, cv);

        //nurbs1.save_denominators("src/den_nurbs1_" + std::to_string(n));
        //nurbs2.save_denominators("src/den_nurbs2_" + std::to_string(n));

        T max_error = 0.0;
        int n_c = nurbs1.coefs.size();

        for (int i = 0; i < n_c; ++i) {
            auto roots = nurbs1.coefs[i].second.solve();
            for (auto & root : roots) {
                T err = std::abs(nurbs1.coefs[i].second.At_complex(root.second)); 
                if (max_error < err) max_error = err;
            }
        }   
        errors1[n - N_start] = max_error;
        
        max_error = 0.0;
        n_c = nurbs2.coefs.size();
        for (int i = 0; i < n_c; ++i) {
            auto roots = nurbs2.coefs[i].second.solve();
            for (auto & root : roots) {
                T err = std::abs(nurbs2.coefs[i].second.At_complex(root.second)); 
                if (max_error < err) max_error = err;
            }
        }   
        errors2[n - N_start] = max_error;
        
    }
    save_vector_errors<T>(NS, errors1, "errors_nurbs1_" + std::to_string(p) + ".out");
    save_vector_errors<T>(NS, errors2, "errors_nurbs2_" + std::to_string(p) + ".out");
    
    plot_errors_(NS, errors1, "Стандартная параметризация");
    //plot_errors(NS, errors2, "compare_errors_den_" + std::to_string(p) + ".png", "Максимальная ошибка корня знаменателя кривых", "nurbs2");
    plot_errors(NS, errors2, "plot" +  std::to_string(p) + ".png", "Максимальная ошибка корня знаменателя кривых", "Параметризация [0, 1]");
    
    PLOT_END();
}

int main(int argc, char const *argv[]) {
    int p = std::atoi(argv[1]);
    test<long double>(p);    
    return 0;
}
