#include "../../nurbs.h"
#include "../../../include/plot.h"

// Сравнение двух методов интегрирования. 
// Кривая степени P для N={P+2, ... } контрольных точек.
// Вектор узлов - uniform

#define N_start 10

template<typename T>
void compare_integration(int p, int N_end) {
    std::cout << N_start << " " << N_end << "\n";
    std::vector<int> vec_N;
    std::vector<T> vec_error1(N_end - N_start + 1, 0.0);
    std::vector<T> vec_error2(N_end - N_start + 1, 0.0);
    
    for(int N = N_start; N <= N_end; ++N) {
        std::cout << N << "\n";
        vec_N.push_back(N);
        std::vector<Point<T>> cv = generate_points_sorted<T>(N);
        //std::vector<T> weights = linspace<T>(T(1), T(2), N);
        std::vector<T> weights = generate_vector<T>(N, T(0.5), T(1.0));
        
        NURBS<T> NURBS(p, weights,cv);
        std::complex<T> int1 = NURBS.analytic_integral1(1);
        std::complex<T> int2 = NURBS.analytic_integral2(1);
        
        T int0(NURBS.numerical_integral());
        T error1 = std::abs(int1 - std::complex<T>(int0));
        T error2 = std::abs(int2 - std::complex<T>(int0));
        vec_error1[N-N_start] = error1;
        vec_error2[N-N_start] = error2;
    }
    save_vector_errors<T>(vec_N, vec_error1, "errors_1.out");
    save_vector_errors<T>(vec_N, vec_error2, "errors_2.out");
    
    plot_errors_(vec_N, vec_error1, "Новый метод");
    plot_errors(vec_N, vec_error2, 
                "compare_int_methods_" + std::to_string(p) + ".png", " ", 
                "Старый метод");
    PLOT_END();
}   


int main(int argc, char const *argv[]) {
    compare_integration<long double>(std::atoi(argv[1]), std::atoi(argv[2]));
    return 0;
}
