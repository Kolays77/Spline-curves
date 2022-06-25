#include "../../bspline.h"
#include "../../../include/plot.h"

#define N_tests 15

#define N_start 10
#define N_end 100

#define P_start 2
#define P_end 5


template <typename T>
void test() {
    std::vector<std::vector<T>> vec_errors1(P_end - P_start + 1, std::vector<T>(N_end - N_start + 1));

    std::vector<int> vec_N;
    
    for(int N = N_start; N <= N_end; ++N) {
        std::cout << N << "\n";
        vec_N.push_back(N);
        std::vector<Point<T>> cv = generate_points_sorted<T>(N);
        for (int p=P_start; p <= P_end; ++p) {  
            for (int k=0; k < N_tests; ++k) {
                std::vector<T> knots = create_knots<T>(N, p);
                Bspline<T>  spline(p, knots, cv);
                T int1_1 = spline.analytic_integral();
                T int2_1 = spline.integral();
                T error1 = std::abs(int1_1 - int2_1);
                vec_errors1[p-P_start][N-N_start] += error1;
            }         
        }
    }

    std::ofstream out1("matrix1");
    out1 << P_start << " " << P_end << " " << N_start << " " << N_end << "\n";

    for (int p=P_start; p <= P_end; ++p) { 
        for(int N = N_start; N <= N_end; ++N) {
            vec_errors1[p-P_start][N-N_start] /= N_tests;
            out1 << vec_errors1[p-P_start][N-N_start] << " ";

        }
        out1 << "\n";
    }

    out1.close();
    PLOT_END();
}


int main() {
    test<long double>();
    return 0;    
}
