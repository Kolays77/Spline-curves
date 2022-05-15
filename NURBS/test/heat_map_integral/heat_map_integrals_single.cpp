#include "../../nurbs.h"
#include "../../nurbs2.h"

#include "../../../include/plot.h"

#define N_tests 20

#define N_start 10
#define N_end 100

#define P_start 2
#define P_end 6


template <typename T>
void test() {
    std::vector<std::vector<T>> vec_errors1(P_end - P_start + 1, std::vector<T>(N_end - N_start + 1));
    std::vector<std::vector<T>> vec_errors2(P_end - P_start + 1, std::vector<T>(N_end - N_start + 1));

    std::vector<int> vec_N;
    
    for(int N = N_start; N <= N_end; ++N) {
        std::cout << N << "\n";
        vec_N.push_back(N);
        std::vector<Point<T>> cv = generate_points_sorted<T>(N);
        for (int p=P_start; p <= P_end; ++p) {  
            for (int k=0; k < N_tests; ++k) {

                // std::vector<T> knots = generate_clamped_random_knot_vector<T>(N, p);
                // std::vector<T> weights = generate_vector<T>(N, 0.5, 1.0);
                
                std::vector<T> weights = generate_vector<T>(N, 0.5, 1.0);


                NURBS<T>  nurbs_1(p, weights,cv);
                NURBS2<T> nurbs_2(p, weights,cv);
                
                std::complex<T> int1_1 = nurbs_1.analytic_integral1(1);
                std::complex<T> int1_2 = nurbs_2.analytic_integral1(1);
                
                T int2_1 = nurbs_1.numerical_integral();
                T int2_2 = nurbs_2.numerical_integral();

                T error1 = std::abs(int1_1 - std::complex<T>(int2_1));
                T error2 = std::abs(int1_2 - std::complex<T>(int2_2));

                vec_errors1[p-P_start][N-N_start] += error1;
                vec_errors2[p-P_start][N-N_start] += error2;
            }         
        }
        PLOT_END();
    }

    std::ofstream out1("matrix1");
    out1 << P_start << " " << P_end << " " << N_start << " " << N_end << "\n";

    std::ofstream out2("matrix2");
    out2 << P_start << " " << P_end << " " << N_start << " " << N_end << "\n";

    for (int p=P_start; p <= P_end; ++p) { 
        for(int N = N_start; N <= N_end; ++N) {
            vec_errors1[p-P_start][N-N_start] /= N_tests;
            vec_errors2[p-P_start][N-N_start] /= N_tests ;
            out1 << vec_errors1[p-P_start][N-N_start] << " ";
            out2 << vec_errors2[p-P_start][N-N_start] << " ";

        }
        out1 << "\n";
        out2 << "\n";
    }

    out1.close();
    out2.close();
    PLOT_END();
}


int main() {
    test<long double>();
    return 0;    
}
