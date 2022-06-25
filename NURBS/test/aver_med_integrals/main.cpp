#include "../../nurbs.h"
#include "../../nurbs2.h"
#include "../../../include/plot.h"

#define N_test 1000


template <typename T>
void test_nurbs_integrals(int p, int n ) {
    std::vector<T> weights = linspace(T(1.0), T(2.0), n);
    std::complex<T> int0_1;
    std::complex<T> int0_2;
    
    std::complex<T> int1;
    std::complex<T> int2;
    
    std::vector<T> errors1(N_test, 0.0);
    std::vector<int> Ns(N_test);

    std::vector<T> errors2(N_test, 0.0);
    for (int i = 0; i < N_test; ++i) {
        Ns[i] = i + 1;
    }

    std::vector<T> x(n+p+1, 0.0);
    for (int i = 0; i < n; ++i) {
        x[i] = T(i);
    }    


    //std::vector<T> knots = generate_clamped_random_knot_vector<T>(n, p);
    for (int i = 0; i < N_test; ++i) { 
        std::cout << i << "\n";
        
        std::vector<T> weights = generate_vector<T>(n, 1.0, 10.0);
        std::vector<T> knots = generate_clamped_random_knot_vector<T>(n, p);
        
        std::vector<Point<T>> cv = generate_points_sorted2<T>(n); 
        
        save_points("cv/" + std::to_string(i), cv);
        
        NURBS<T> nurbs1(p, knots, weights, cv);
        NURBS2<T> nurbs2(p, knots, weights, cv);

        nurbs1.save_denominators("dens/1_" + std::to_string(i));
        nurbs2.save_denominators("dens/2_" + std::to_string(i));



        int0_1 = std::complex<T>{nurbs1.numerical_integral()};
        int0_2 = std::complex<T>{nurbs2.numerical_integral()};
        
        int1 = nurbs1.analytic_integral1(1);
        int2 = nurbs2.analytic_integral1(1);
        errors1[i] = std::abs(int0_1 - int1);
        errors2[i] = std::abs(int0_2 - int2);   
    }
    
    T av1 = count_average_error(errors1);
    T av2 = count_average_error(errors2);

    T median1 = count_median(errors1);
    T median2 = count_median(errors2);

    std::cout << "Average errors1 : " <<  av1 << "\n";
    std::cout << "Average errors2 : " <<  av2 << "\n\n";

    std::cout << "Median errors1 : " <<  median1 << "\n";
    std::cout << "Median errors2 : " <<  median2 << "\n\n";
    
    std::cout << "Min errors1 : " << *std::min_element(errors1.begin(), errors1.end()) << "\n";
    std::cout << "Min errors1 : " << *std::min_element(errors2.begin(), errors2.end()) << "\n\n";
    
    std::cout << "Max errors1 : " << *std::max_element(errors1.begin(), errors1.end()) << "\n";
    std::cout << "Max errors1 : " << *std::max_element(errors2.begin(), errors2.end()) << "\n";
    save_vector<T>(errors1, "errors1.txt");
    save_vector<T>(errors2, "errors2.txt");

    std::map<std::string, std::string> keywords_av;
    keywords_av["color"] = "red";
    keywords_av["label"] = "Среднее";

    std::map<std::string, std::string> keywords_median;
    keywords_median["color"] = "orange";
    keywords_median["label"] = "Медиана";
    
    plt::axhline(av1, 0, N_test, keywords_av);
    plt::axhline(median1, 0, N_test, keywords_median);
    plot_errors(Ns, errors1, "errors1.png", " Ошибки вычисления аналитического интеграла для NURBS1", 
                    "p : " + std::to_string(p) + " n : " + std::to_string(n));

    plt::axhline(av2, 0, N_test, keywords_av);
    plt::axhline(median2, 0, N_test, keywords_median);
    plot_errors(Ns, errors2, "errors2.png", " Ошибки вычисления аналитического интеграла для NURBS2", 
                    "p : " + std::to_string(p) + " n : " + std::to_string(n));

    PLOT_END();   
}


int main(int argc, char const *argv[]) {
    int p, n;
    p = std::atoi(argv[1]);
    n = std::atoi(argv[2]);
    test_nurbs_integrals<long double>(p, n);    
    return 0;
}
