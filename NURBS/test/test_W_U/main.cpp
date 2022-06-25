#include "../../nurbs.h"
#include "../../nurbs2.h"
#include "../../../include/plot.h"

template <typename T>
void plot_den(NURBS2<T>& nurbs, 
            const std::string& path) {
    
    const int N = 100;
    std::vector<T> ts = linspace(T(0.0), T(1.0), N);        
    std::vector<T> ts_t;

    std::vector<T> points(N);  
 	plt::xlabel("X", {{"fontsize", "large"}});
    plt::ylabel("Y", {{"fontsize", "large"}});
     
    for(int j = 0; j < nurbs.N_segments - 1; ++j) {
        ts_t = linspace(T(j) / nurbs.N_segments, T(j+1) / nurbs.N_segments, N);
        for (int i = 0; i < N; ++i) {
            points[i] = nurbs.coefs[j].second.At(ts[i]);
        }
        plot_curve_<T>(ts_t, points);
    }
    
    for (int i = 0; i < N; ++i) {
        points[i] = nurbs.coefs[nurbs.N_segments - 1].second.At(ts[i]);
    }   
    ts_t = linspace(T(nurbs.N_segments - 1) / nurbs.N_segments, T(1.0), N);      
 	plt::xlabel("t", {{"fontsize", "large"}});
    plt::ylabel("Q(t)", {{"fontsize", "large"}});
        
	plot_curve<T>(ts_t, points, path, "Знаменатель Q(t) NURBS кривой");
}

template <typename T>
void plot_den_curvature(NURBS2<T>& nurbs, 
            const std::string& path) {
    
    const int N = 100;
    std::vector<T> ts = linspace(T(0.0), T(1.0), N);        
    std::vector<T> ts_t;

    std::vector<T> k_t(N);  

    for(int j = 0; j < nurbs.N_segments- 1; ++j) {
        ts_t = linspace(T(j) / nurbs.N_segments, T(j+1) / nurbs.N_segments, N);
        
        for (int i = 0; i < N; ++i) {
            Poly<T> q_der = nurbs.coefs[j].second.der();
            T t = ts[i];
            k_t[i] = std::abs(q_der.der().At(t)) 
                        / std::pow( 1 + std::pow(q_der.At(t), 2), 1.5);
        }
        plt::xlabel("t", {{"fontsize", "large"}});
        plt::ylabel("Кривизна Q(t)", {{"fontsize", "large"}});
        plot_curve_<T>(ts_t, k_t);
    }
    
    for (int i = 0; i < N; ++i) {
        Poly<T> q_der = nurbs.coefs[nurbs.N_segments - 1].second.der();
        T t = ts[i];
        k_t[i] = std::abs(q_der.der().At(t)) 
                    / std::pow( 1 + std::pow(q_der.At(t), 2), 1.5);
    }   
    ts_t = linspace(T(nurbs.N_segments - 1) / nurbs.N_segments, T(1.0), N);      
    plt::xlabel("t", {{"fontsize", "large"}});
    plt::ylabel("Кривизна Q(t)", {{"fontsize", "large"}});
    plot_curve<T>(ts_t, k_t, path, "График кривизны Q(t) NURBS кривой");
}



template <typename T>
void plot_den(NURBS<T>& nurbs, 
            const std::string& path) {
    int N = 100;
    std::vector<T> points(N);
    std::vector<T> ts = linspace(nurbs.knots[nurbs.domain.first], nurbs.knots[nurbs.domain.second], N);
    T t = nurbs.knots[nurbs.domain.first];
    int i = 0;
    for(int j = 0; j < nurbs.N_segments; ++j) {
        while (t <= nurbs.knots[nurbs.ks[j]+1] && i < N) {
            points[i] = nurbs.coefs[j].second.At(t);
            ++i; t = ts[i];
        }
    }
    plt::xlabel("t", {{"fontsize", "large"}});
    plt::ylabel("Q(t)", {{"fontsize", "large"}});
    plot_curve<T>(ts, points, path, "Знаменатель Q(t) NURBS кривой");
}






template <typename T>
void test_nurbs_integrals(int p, int n ) {
    std::vector<Point<T>> cv = generate_points_sorted<T>(n);
    std::vector<T> weights1 = linspace(T(1.0), T(2.0), cv.size());
    std::vector<T> weights2 = generate_vector<T>(n, 0.5, 1.0);
    std::vector<T> knots = generate_clamped_random_knot_vector<T>(n, p);

    NURBS2<T> nurbs1(p, 1.0, 2.0, cv); // w U -> линейная часть
    NURBS2<T> nurbs2(p, knots, weights1, cv); // w -> 
    NURBS2<T> nurbs3(p, weights2, cv); // U - 
    NURBS2<T> nurbs4(p, knots, weights2, cv); // - - 


    nurbs1.save_coefs("1/");
    nurbs2.save_coefs("2/");
    nurbs3.save_coefs("3/");
    nurbs4.save_coefs("4/");
        
    std::vector<Point<T>> points_1 = nurbs1.get_points(1000);
    std::vector<Point<T>> points_2 = nurbs2.get_points(1000);
    std::vector<Point<T>> points_3 = nurbs3.get_points(1000);
    std::vector<Point<T>> points_4 = nurbs4.get_points(1000);
    
    plot_den(nurbs1, "1/den.png");
    plot_den(nurbs2, "2/den.png");
    plot_den(nurbs3, "3/den.png");
    plot_den(nurbs4, "4/den.png");

    plot_den_curvature(nurbs1, "1/den_curve.png");
    plot_den_curvature(nurbs2, "2/den_curve.png");
    plot_den_curvature(nurbs3, "3/den_curve.png");
    plot_den_curvature(nurbs4, "4/den_curve.png");


 	plt::xlabel("X", {{"fontsize", "large"}});
    plt::ylabel("Y", {{"fontsize", "large"}});
    plot_curve(cv, points_1, "1/plot.png", "NURBS Кривая", "p =" + std::to_string(p));

 	plt::xlabel("X", {{"fontsize", "large"}});
    plt::ylabel("Y", {{"fontsize", "large"}});
    plot_curve(cv, points_2, "2/plot.png", "NURBS Кривая", "p =" + std::to_string(p));    

 	plt::xlabel("X", {{"fontsize", "large"}});
    plt::ylabel("Y", {{"fontsize", "large"}});    
	plot_curve(cv, points_3, "3/plot.png", "NURBS Кривая", "p =" + std::to_string(p));

 	plt::xlabel("X", {{"fontsize", "large"}});
    plt::ylabel("Y", {{"fontsize", "large"}});
    plot_curve(cv, points_4, "4/plot.png", "NURBS Кривая", "p =" + std::to_string(p));
    PLOT_END(); 
}


int main(int argc, char const *argv[]) {
    int p = std::atoi(argv[1]);
    int n = std::atoi(argv[2]);
    test_nurbs_integrals<long double>(p, n);    
    return 0;
}
