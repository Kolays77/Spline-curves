#include "../../nurbs2.h"
#include "../../nurbs.h"

#include "../../../include/plot.h"

#include <thread>
#include <mutex>


#define N_tests 2
#define N_threads 8

#define N_start 10  
#define N_end 100

#define P_start 2
#define P_end 6

template <typename T>
void worker(std::vector<T> & vec_errors1, 
            std::vector<T> & vec_errors2,
            int k,
            int p, 
            int N) {

    std::vector<T> knots = generate_clamped_random_knot_vector<T>(N, p);
    // std::vector<T> weights = generate_vector<T>(N, T(0.5), T(1.0));
    
    std::vector<Point<T>> cv = generate_points_sorted<T>(N);
    
    NURBS<T>  nurbs_1(p, knots, 1.0, 2.0, cv);
    NURBS2<T> nurbs_2(p, knots, 1.0, 2.0, cv);
    
    std::complex<T> int1_1 = nurbs_1.analytic_integral1(1);
    std::complex<T> int1_2 = nurbs_2.analytic_integral1(1);
    
    T int2_1 = nurbs_1.numerical_integral();
    T int2_2 = nurbs_2.numerical_integral();

    T error1 = std::abs(int1_1 - std::complex<T>(int2_1));
    T error2 = std::abs(int1_2 - std::complex<T>(int2_2));

    vec_errors1[k] = error1;
    vec_errors2[k] = error2;    
}


template <typename T>
void integral_fixed_deg() {
    std::cout << "Number of threads : " << N_threads << "\n";
    // integral with increasing number of points
    // для нового числа точек происходит перегенерация
    
    std::vector<std::vector<T>> vec_errors1(P_end - P_start + 1, std::vector<T>(N_end - N_start + 1));
    std::vector<std::vector<T>> vec_errors2(P_end - P_start + 1, std::vector<T>(N_end - N_start + 1));

    std::vector<int> vec_N;
    
    std::vector<std::thread> vec_thread(N_threads); 
    
    std::vector<T> vec_temp_errors1(N_threads, 0.0);
    std::vector<T> vec_temp_errors2(N_threads, 0.0);

    
    for(int N = N_start; N <= N_end; ++N) {
        std::cout << N << "\n";
        vec_N.push_back(N);

        for (int p=P_start; p <= P_end; ++p) {  
            
            // Pool threading  N_tests * N_threads
            for (int i = 0; i < N_tests; ++i) {                    

                for (int k_thread=0; k_thread < N_threads; ++k_thread) {
                    vec_thread[k_thread] = std::thread(worker<T>, 
                                                    std::ref(vec_temp_errors1), 
                                                    std::ref(vec_temp_errors2), 
                                                    k_thread,
                                                    p, 
                                                    N);
                }

                for (int k_thread=0; k_thread < N_threads; ++k_thread) {
                    vec_thread[k_thread].join();
                    vec_errors1[p-P_start][N-N_start] += vec_temp_errors1[k_thread];
                    vec_errors2[p-P_start][N-N_start] += vec_temp_errors2[k_thread];

                }
            }
        }
    }

    std::ofstream out1("matrix1");
    out1 << P_start << " " << P_end << " " << N_start << " " << N_end << "\n";

    std::ofstream out2("matrix2");
    out2 << P_start << " " << P_end << " " << N_start << " " << N_end << "\n";

    for (int p=P_start; p <= P_end; ++p) { 
        for(int N = N_start; N <= N_end; ++N) {
            vec_errors1[p-P_start][N-N_start] /= (N_tests *  N_threads);
            vec_errors2[p-P_start][N-N_start] /= (N_tests * N_threads);

            out1 << vec_errors1[p-P_start][N-N_start] << " ";
            out2 << vec_errors2[p-P_start][N-N_start] << " ";

        }
        out1 << "\n";
        out2 << "\n";
        
    }
    out1.close();
    out2.close();

    for (int p=P_start; p < P_end; ++p) { 
        plot_errors_(vec_N, vec_errors1[p-P_start], "deg = " + std::to_string(p));
    }
    plot_errors(vec_N, vec_errors1[P_end-P_start], "errors1.png", "Сравнение численного и аналитического интеграла", "deg = " + std::to_string(P_end));
        
    for (int p=P_start; p < P_end; ++p) { 
        plot_errors_(vec_N, vec_errors2[p-P_start], "deg = " + std::to_string(p));
    }
    plot_errors(vec_N, vec_errors2[P_end-P_start], "errors2.png", "Сравнение численного и аналитического интеграла", "deg = " + std::to_string(P_end));
    
    PLOT_END();
}


int main() {
    integral_fixed_deg<long double>();
    return 0;    
}
