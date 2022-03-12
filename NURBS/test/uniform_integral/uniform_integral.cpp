#include "../../nurbs.h"
#include "../../../include/plot.h"

#define N_start 10
#define N_end 100

#define P_start 2
#define P_end 8

template <typename T>
void integral_fixed_deg() {
    // integral with increasing number of points
    // для нового числа точек происходит перегенерация

    std::vector<std::vector<T>> vec_errors(P_end - P_start + 1, std::vector<T>(N_end - N_start + 1));
    std::vector<int> vec_N;
    
    for(int N = N_start; N <= N_end; ++N) {
        std::cout << N << "\n";
        vec_N.push_back(N);
        std::vector<Point<T>> cv = generate_points_sorted<T>(N);
        std::vector<T> weights = linspace(T(1.0), T(2.0), N);

        for (int p=P_start; p <= P_end; ++p) {        
            NURBS<T> NURBS(p, weights,cv);
            std::complex<T> int1 = NURBS.analytic_integral1(1);
            T int2(NURBS.numerical_integral());
            
            T error = std::abs(int1 - std::complex<T>(int2));
            vec_errors[p-P_start][N-N_start] = error;
        }
    }
    for (int p=P_start; p < P_end; ++p) { 
        save_vector_errors<T>(vec_N, vec_errors[p-P_start], "errors_" + std::to_string(p) + ".out");
        plot_errors_(vec_N, vec_errors[p-P_start], "deg = " + std::to_string(p));
    }

    save_vector_errors<T>(vec_N, vec_errors[P_end-P_start], "errors_" + std::to_string(P_end) + ".out");
    plot_errors(vec_N, vec_errors[P_end-P_start], "errors_fixed_deg.png", "Вектор весов в интервале [1.0, 2.0]", "deg = " + std::to_string(P_end));
    
    PLOT_END();
}


int main() {
    
    integral_fixed_deg<long double>();
    return 0;    
}