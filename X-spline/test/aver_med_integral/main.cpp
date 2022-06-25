#include "../../xspline.h"
#include "../../../include/plot.h"

#define N_repeat 10

#define N_start 4
#define N_end 500

template<typename T>
void plot_errors_av(const std::vector<int>& Ns, const std::vector<T>& errors, std::string path) {
    plt::ylabel("Средняя ошибка", {{"fontsize", "x-large"}});
    plt::xlabel("Число сегментов сплайна", {{"fontsize", "x-large"}});
    plot_errors<T>(Ns, errors , path);
}


template<typename T>
void plot_errors_med(const std::vector<int>& Ns, const std::vector<T>& errors, std::string path) {
    plt::ylabel("Медианная ошибка", {{"fontsize", "x-large"}});
    plt::xlabel("Число сегментов сплайна", {{"fontsize", "x-large"}});
    plot_errors<T>(Ns, errors , path);
}


template<typename T>
void test() {
    // Вектор s_k -- случайный вектор. s_k in [-1, 1] 
    std::vector<T> errors_av1(N_end - N_start + 1, 0.0);
    std::vector<T> errors_med1(N_end - N_start + 1, 0.0);

    // Вектор s_k = 0  для всех k
    std::vector<T> errors_av2(N_end - N_start + 1, 0.0);
    std::vector<T> errors_med2(N_end - N_start + 1, 0.0);

    // Вектор s_k = -1  для всех k
    std::vector<T> errors_av3(N_end - N_start + 1, 0.0);
    std::vector<T> errors_med3(N_end - N_start + 1, 0.0);

    // Вектор s_k = 1  для всех k
    std::vector<T> errors_av4(N_end - N_start + 1, 0.0);
    std::vector<T> errors_med4(N_end - N_start + 1, 0.0);

    std::vector<int> Ns;
    for (int n = N_start; n <= N_end; ++n) {
        std::cout << n << "\n";
        Ns.push_back(n);
        std::vector<T> errors_rep1(N_repeat);
        std::vector<T> errors_rep2(N_repeat);
        std::vector<T> errors_rep3(N_repeat);
        std::vector<T> errors_rep4(N_repeat);

        std::vector<T> s2(n, -1.0);
        std::vector<T> s3(n,  0.0);
        std::vector<T> s4(n,  1.0);
        
        for (int j = 0; j < N_repeat; ++j) {
            std::vector<T> cv_x = generate_vector_sorted<T>(n, 0.0, 1.0);
            std::vector<T> cv_y = generate_vector<T>(n, 0.0, 1.0);

            std::vector<T> s1 = generate_vector<T>(n, -1.0, 1.0);
        
            Xspline<T> spline1(cv_x, cv_y, s1);
            Xspline<T> spline2(cv_x, cv_y, s2);
            Xspline<T> spline3(cv_x, cv_y, s3);
            Xspline<T> spline4(cv_x, cv_y, s4);
        
            std::complex<T> integral_1 = spline1.numerical_integral();
            std::complex<T> integral_2 = spline1.analytical_integral();
            errors_rep1[j] = std::abs(integral_1 - integral_2);

            integral_1 = spline2.numerical_integral();
            integral_2 = spline2.analytical_integral();
            errors_rep2[j] = std::abs(integral_1 - integral_2);


            integral_1 = spline3.numerical_integral();
            integral_2 = spline3.analytical_integral();
            errors_rep3[j] = std::abs(integral_1 - integral_2);

            integral_1 = spline4.numerical_integral();
            integral_2 = spline4.analytical_integral();
            errors_rep4[j] = std::abs(integral_1 - integral_2);
        }  
                 
        std::sort(errors_rep1.begin(), errors_rep1.end());
        std::sort(errors_rep2.begin(), errors_rep2.end());
        std::sort(errors_rep3.begin(), errors_rep3.end());
        std::sort(errors_rep4.begin(), errors_rep4.end());
        
        
        T res1 = std::accumulate(errors_rep1.begin(), errors_rep1.end(), T(0.0));
        T res2 = std::accumulate(errors_rep2.begin(), errors_rep2.end(), T(0.0));
        T res3 = std::accumulate(errors_rep3.begin(), errors_rep3.end(), T(0.0));
        T res4 = std::accumulate(errors_rep4.begin(), errors_rep4.end(), T(0.0));
        
        errors_med1[n - N_start] = errors_rep1[N_repeat / 2];
        errors_av1[n - N_start] = res1 / T(N_repeat);
    
        errors_med2[n - N_start] = errors_rep2[N_repeat / 2];
        errors_av2[n - N_start] = res2 / T(N_repeat);
        
        errors_med3[n - N_start] = errors_rep3[N_repeat / 2];
        errors_av3[n - N_start] = res3 / T(N_repeat);
        
        errors_med4[n - N_start] = errors_rep4[N_repeat / 2];
        errors_av4[n - N_start] = res4 / T(N_repeat);
    }
    plot_errors_av(Ns, errors_av1, "1/av.png");
    plot_errors_av(Ns, errors_av2, "2/av.png");
    plot_errors_av(Ns, errors_av3, "3/av.png");
    plot_errors_av(Ns, errors_av4, "4/av.png");

    plot_errors_med(Ns, errors_med1, "1/med.png");
    plot_errors_med(Ns, errors_med2, "2/med.png");
    plot_errors_med(Ns, errors_med3, "3/med.png");
    plot_errors_med(Ns, errors_med4, "4/med.png");
    PLOT_END();
}

int main(int argc, char* argv[]) {
    test<long double>();
}

