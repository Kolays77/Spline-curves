#include "../../nurbs.h"
#include "../../nurbs2.h"

#define P_start 2
#define P_end 8


#define N_start 12
#define N_end 100

template <typename T>
T max_element(Poly<T>& poly) {
    T max = 0.0;
    int deg = poly.deg();
    for (int i = 0; i <= deg; ++i) {
        if (max < std::abs(poly[i])) max = std::abs(poly[i]); 
    }
    return max;
}

template <typename T>
void test() {
    
    std::ofstream out1("matrix1");
    std::ofstream out2("matrix2");
    out1 << P_start << " " << P_end << " " << N_start << " " << N_end << "\n";
    out2 << P_start << " " << P_end << " " << N_start << " " << N_end << "\n";
    
    for (int p=P_start; p <= P_end; ++p) {  
        for (int n = N_start; n <= N_end; ++n) {
            std::cout << "(" << p << ", " << n << ")\n";
            std::vector<Point<T>> cv = generate_points_sorted<T>(n);
            
            //std::vector<T> weights = linspace(T(1.0), T(2.0), cv.size()); // знаменатель не зависит от вектора узлов
            std::vector<T> knots = generate_clamped_random_knot_vector<T>(n, p);
            std::vector<T> weights = generate_vector<T>(n, T(1.0), T(2.0));
            
            // NURBS<T>  nurbs1(p, 1.0, 2.0, cv);    
            // NURBS2<T> nurbs2(p, 1.0, 2.0, cv);

            NURBS<T>  nurbs1(p, knots, weights, cv);    
            NURBS2<T> nurbs2(p, knots, weights, cv);



            T max_value = 0.0;
            int n_c = nurbs1.coefs.size();
            for (int i = 0; i < n_c; ++i) {
                
                T max_temp_den = max_element(nurbs1.coefs[i].second);
                T max_temp_x = max_element(nurbs1.coefs[i].first[0]);
                T max_temp_y = max_element(nurbs1.coefs[i].first[1]);

                if (max_value < max_temp_den) max_value = max_temp_den;
                if (max_value < max_temp_x) max_value = max_temp_x;
                if (max_value < max_temp_y) max_value = max_temp_y;
            
            }   
            out1 << max_value << " ";
            
            max_value = 0.0;
            n_c = nurbs2.coefs.size();
            for (int i = 0; i < n_c; ++i) {
                T max_temp_den = max_element(nurbs2.coefs[i].second);
                T max_temp_x = max_element(nurbs2.coefs[i].first[0]);
                T max_temp_y = max_element(nurbs2.coefs[i].first[1]);

                if (max_value < max_temp_den) max_value = max_temp_den;
                if (max_value < max_temp_x) max_value = max_temp_x;
                if (max_value < max_temp_y) max_value = max_temp_y;
            }   
            out2 << max_value<< " ";
        }
        out1 << "\n";
        out2 << "\n";
    }    
    out1.close();
    out2.close();
}

int main(int argc, char const *argv[]) {
    test<long double>();    
    return 0;
}
