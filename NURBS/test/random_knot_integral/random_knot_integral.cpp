#include "../../nurbs2.h"
#include "../../../include/plot.h"

#include <thread>
#include <mutex>

const int N_threads = std::thread::hardware_concurrency();

#define N_tests 10  

#define N_start 7   
#define N_end 100

#define P_start 2
#define P_end 7


template <typename T>
void worker(std::vector<T> & vec_errors, 
            int k,
            int p, 
            int N, 
            std::vector<Point<T>>& cv, 
            std::vector<T>& weights) {

    std::vector<T> knots = generate_clamped_random_knot_vector<T>(N, p);
    NURBS2<T> NURBS(p, knots, weights,cv);
    std::complex<T> int1 = NURBS.analytic_integral1(1);
    T int2(NURBS.numerical_integral());
    T error = std::abs(int1 - std::complex<T>(int2));
    vec_errors[k] = error;
}


template <typename T>
void integral_fixed_deg() {
    std::cout << "Number of threads : " << N_threads << "\n";
    // integral with increasing number of points
    // для нового числа точек происходит перегенерация
    
    std::vector<std::vector<T>> vec_errors(P_end - P_start + 1, std::vector<T>(N_end - N_start + 1));
    std::vector<int> vec_N;
    
    std::vector<std::thread> vec_thread(N_threads); 
    std::vector<T> vec_temp_errors(N_threads, 0.0);
    
    for(int N = N_start; N <= N_end; ++N) {
        std::cout << N << "\n";
        vec_N.push_back(N);
        std::vector<Point<T>> cv = generate_points_sorted<T>(N);
        std::vector<T> weights = generate_vector<T>(N, T(0.00001), T(10.0));

        for (int p=P_start; p <= P_end; ++p) {  
            
            // Pool threading 
            for (int i = 0; i < N_tests; ++i) {    
                for (int k_thread=0; k_thread < N_threads; ++k_thread) {
                    vec_thread[k_thread] = std::thread(worker<T>, 
                                                    std::ref(vec_temp_errors), 
                                                    k_thread,
                                                    p, 
                                                    N, 
                                                    std::ref(cv), 
                                                    std::ref(weights));
                }
                for (int k_thread=0; k_thread < N_threads; ++k_thread) {
                    vec_thread[k_thread].join();
                    vec_errors[p-P_start][N-N_start] += vec_temp_errors[k_thread];
                }
            }
        }
    }

    for(int N = N_start; N <= N_end; ++N) {
        for (int p=P_start; p <= P_end; ++p) { 
            vec_errors[p-P_start][N-N_start] /= N_tests;
        }
    }

    for (int p=P_start; p < P_end; ++p) { 
        save_vector_errors<T>(vec_N, vec_errors[p-P_start], "errors_" + std::to_string(p) + ".out");
        plot_errors_(vec_N, vec_errors[p-P_start], "deg = " + std::to_string(p));
    }

    save_vector_errors<T>(vec_N, vec_errors[P_end-P_start], "errors_" + std::to_string(P_end) + ".out");
    plot_errors(vec_N, vec_errors[P_end-P_start], "errors_fixed_deg.png", "Сравнение численного и аналитического интеграла", "deg = " + std::to_string(P_end));
    
    PLOT_END();
}


int main() {
    integral_fixed_deg<long double>();
    return 0;    
}
