#pragma once
#include "../include/Point.h"
#include "../include/Poly.h"
#include "../include/tools.h"
#include "../include/integral.h"


/*
Модифицированный алгоритм де Бура для построения (получение коэффициентов) рациональных функций (P_x/Q(t), P_y/Q(t)), 
определяющих k-ый сегмент NURBS кривой. Посегментная параметризация 
Параметры:
    k1, k2 - индексы начала и конца интервала изменения параметра t.
    knots - вектор узлов;
    points - вектор контрольных точек (CV) на плоскости;
    p - степень кривой.
*/

template<typename T>
Rational<T> de_boor_nurbs2(int k1, int k2,
                          std::vector<T>& knots,
                          std::vector<T>& weights,
                          std::vector<Point<T>>& points,
                          int p) {

    size_t dim = points[0].dim;
    std::vector<Point<Poly<T>>> d((size_t)p+1, Point<Poly<T>>(dim));
    std::vector<Poly<T>> d_den((size_t)p+1);

    for (int i = 0; i < p+1; ++i){
        d_den[i][0] = weights[i + k1 - p];
        for (int j = 0; j < dim; ++j)
            d[i][j][0] = points[i + k1 - p][j] * weights[i + k1 - p];
    }

    T B,D, den;

    Poly<T> poly_up({knots[k2] - knots[k1], knots[k1]});
    Poly<T> poly_up_neg({knots[k1] - knots[k2], - knots[k1]});
 
    
    for (int r = 1; r < p+1; ++r) {
        for (int j=p ; j > r-1; --j) {
            den = knots[j+1+k1-r] - knots[j+k1-p];
            B = - knots[j+k1-p];
            D = knots[j+1+k1-r];
            Poly<T> temp1_den = d_den[j] * poly_up + B * d_den[j];
            Poly<T> temp2_den = d_den[j-1] * poly_up_neg + D * d_den[j-1];
            d_den[j] = (temp1_den + temp2_den) / den;
            for (size_t i = 0; i < dim; ++i){
                Poly<T> temp1 = d[j][i] * poly_up + B * d[j][i] ;
                Poly<T> temp2 = d[j-1][i] * poly_up_neg + D * d[j-1][i];
                d[j][i] = (temp1 + temp2) / den;
            }
        }
    }
    return {d[p], d_den[p]};
}


/*
    Класс NURBS кривой, в конструкторах которого строится аналитические выражение(многочлены) С(t) = (P_x/Q(t), P_y/Q(t))
    для каждого сегмента кривой.
    Посегметная параметризация (Параметризация [0, 1]). 
    Параметры:
        p - степень кривой;
        knots - вектор узлов;
        cv - вектор контрольных точек; 
        weights - веса.
*/
template<typename T>
struct NURBS2{
    int p;
    int dim;
    int N_segments; // ks.size() - 1

    std::vector<T> knots;
    std::vector<T> weights;
    std::vector<Point<T>> cv;
    std::pair<int,int> domain;
    std::vector<int> ks; // indexes of knots (consist 0.0 and 1.0)
    std::vector<Rational<T>> coefs;

public:

    // Конструктор для построение NURBS кривой.
    // Условие: неравномерный вектор узлов (NON-uniform KNOT VECTOR)
    // неравномерный вектор весов (NON-uniform weights)
    NURBS2(int p_,
          std::vector<T>& knots_,
          std::vector<T>& weights_,
          std::vector<Point<T>>& cv_) {
        p = p_;
        knots = knots_;
        weights = weights_;
        cv = cv_;
        dim = cv[0].dim;
        domain = {p, knots.size() - p - 1};
        ks = create_intervals(domain,knots);
        N_segments = ks.size() - 1;
        create_coefs();
        // polishing();
    }


    // Конструктор для построение NURBS кривой. 
    // Условие: равномерный вектор узлов (UNIFORM KNOT VECTOR)
    // неравномерный вектор весов (NON-uniform weights)
    NURBS2(int p_,
          std::vector<T>& weights_,
          std::vector<Point<T>>& cv_) {
        p = p_;
        cv = cv_;
        weights = weights_;
        knots = create_knots<T>(cv.size(), p);
        dim = cv[0].dim;
        domain = {p, knots.size() - p - 1};
        ks = create_intervals(domain,knots);
        N_segments = ks.size() - 1;
        create_coefs();
        // polishing();
    }

    // Конструктор для построение NURBS кривой. 
    // Условие: равномерный вектор узлов (UNIFORM KNOT VECTOR)
    // равномерный вектор весов (UNIFORM weights)
    NURBS2(int p_,  
            T w_start, 
            T w_end,
          std::vector<Point<T>>& cv_) {
        p = p_;
        cv = cv_;
        weights = linspace<T>(w_start, w_end, cv.size());
        knots = create_knots<T>(cv.size(), p);
        dim = cv[0].dim;
        domain = {p, knots.size() - p - 1};
        ks = create_intervals(domain,knots);
        N_segments = ks.size() - 1;
        create_coefs();
        // polishing_uniform();
    }

    // Конструктор для построение NURBS кривой. 
    // Условие: неравномерный вектор узлов (NON UNIFORM KNOT VECTOR)
    // равномерный вектор весов (UNIFORM weights)
    NURBS2(int p_,   
            std::vector<T>& knots_,
            T w_start, 
            T w_end,
          std::vector<Point<T>>& cv_) {
        p = p_;
        cv = cv_;
        weights = linspace<T>(w_start, w_end, cv.size());
        knots = knots_;
        dim = cv[0].dim;
        domain = {p, knots.size() - p - 1};
        ks = create_intervals(domain,knots);
        N_segments = ks.size() - 1;
        create_coefs();
        // polishing();
    }
    
    void create_coefs(){  
        for (int i = 0; i < N_segments; ++i) 
            coefs.push_back(de_boor_nurbs2(ks[i], ks[i+1], knots,weights, cv, p));
    }

    // Метод удаление старших коэффициентов при выполнения условия A < eps
    void polishing(T eps_ = 1e-11) {
        for (auto& frac: coefs){
            frac.second.update_real(eps_);
            for (int d=0; d < dim; ++d){
                frac.first[d].update_real(eps_);
            }
        }
    }

    // Нормализация знаменателя Q(t)
    void polishing2() {
        for (auto& frac: coefs){
            T temp = frac.second.normalize();
            for (int d=0; d < dim; ++d){
               frac.first[d] /= temp;
            }
        }
    } 
    
    void polishing_uniform(T eps_ = 1e-11) {
        int len = coefs.size();
        if (len <= 2*(p-1)) {
            polishing();
        } else {
            for (int i = 0; i < p-1; ++i){
                coefs[i].second.update_real(eps_);
                coefs[len - i - 1].second.update_real(eps_);
                for (int d=0; d < dim; ++d){
                    coefs[i].first[d].update_real(eps_);
                    coefs[len-i-1].first[d].update_real(eps_);
                }
            }
            for (int i = p-1; i < len-p+1; ++i) {
                coefs[i].second.resize(1);
                for (int d=0; d < dim; ++d){
                    coefs[i].first[d].update_real(eps_);
                }
            }
        } 
    }


    // Метод получения N точек на кривой
    std::vector<Point<T>> get_points(int N) {
        std::vector<T> ts = linspace(T(0.0), T(1.0), N / N_segments);        
        std::vector<Point<T>> points(ts.size() * N_segments, Point<T>(dim));
        
        for(int j = 0; j < N_segments; ++j){
            for (int i = 0; i < ts.size(); ++i)             
                for (int d = 0; d < dim; ++d) {
                    points[j * ts.size() + i][d] = coefs[j].first[d].At(ts[i]) / coefs[j].second.At(ts[i]) ;
                }
            }
        return  points;
    }

    // Метод нахождения производной по икс NURBS кривой.
    // Результат:
    // N - штук значений dy/dx
    std::vector<T> get_dy_dx_points(int N) {
        std::vector<T> points(N, 0.0);
        std::vector<T> ts = linspace(knots[domain.first],
                                     knots[domain.second], N);
        T t = knots[domain.first];
        int i = 0;
        for (int j = 0; j < N_segments; ++j) {
            Poly<T> p_x = coefs[j].first[0];
            Poly<T> p_y = coefs[j].first[1];
            Poly<T> den = coefs[j].second;

            while (t <= knots[ks[j+1]] && i < N) {
                points[i] = (p_y.der().At(t) * den.At(t) - p_y.At(t) *den.der().At(t)) /
                            (p_x.der().At(t) * den.At(t) - p_x.At(t) *den.der().At(t));
                ++i;
                t = ts[i];
            }
        }
        return points;
    }

    // Сохранение коэффициентов знаменателя для каждого сегмента в файл. 
    void save_denominators(std::string path = "") {
        std::ofstream out_den(path);
        for (const Rational<T>& coef : coefs){
            out_den << coef.second << "\n";
        }
    }
    // Сохранение коэффициентов многочленов для каждого сегмента в файл. 
    void save_coefs(std::string path_dir = "") {
        std::string path_x="coefs2_num_x.out";
        std::string path_y="coefs2_num_y.out";
        std::string den = "coefs2_den.out";

        std::ofstream out_x(path_dir + path_x);
        std::ofstream out_y(path_dir + path_y);
        std::ofstream out_den(path_dir + den);

        for (int i = 0; i < N_segments; ++i) {
            out_x << coefs[i].first[0] << "\n";
            out_y << coefs[i].first[1] << "\n";
            out_den << "[" << knots[ks[i]] << ":" << knots[ks[i+1]] << "]" << coefs[i].second << "\n";
        }
    }

    // Численный интеграл. Квадратура Гаусса-Лежандра
    T numerical_integral() {
        static const std::vector<long double> ps = POINTS;
        static const std::vector<long double> ws = WEIGHTS;
        int n = ps.size();
        T sum = 0.0;        
        T A = 0.0;
        T B = 1.0;            
            
        for (int i = 0; i < N_segments; ++i) {
    
            Poly<T> p_x = coefs[i].first[0];
            Poly<T> p_y = coefs[i].first[1];
            Poly<T> den = coefs[i].second;
            T sum_temp = 0.0;
            for (int j = 0; j < n; ++j) {
                T x = T(B-A)*ps[j]/2 + T(A+B)/2;
                // T x = ps[j] / 2 + T(0.5);
                
                T den_quadratic = std::pow(den.At(x), 2);
                sum_temp += ws[j] * p_y.At(x) * p_x.der().At(x) / den_quadratic;
                sum_temp -= ws[j] * p_y.At(x) * den.der().At(x)*p_x.At(x) / den_quadratic / den.At(x); //
            }
            sum += (B-A)/2.0 *sum_temp;
        }
        return sum;
    }

    // Аналитический интеграл (Алгоритм 2)
    // type -- метод факторизации знаменателя Q(t)
    //      1 - с помощью вспомогательной матрицы
    //      2 - метод Лаггера 
    std::complex<T> analytic_integral1(int type = 1, T t0=0.0, T t1=1.0) {
        polishing2();
        std::complex<T> sum = 0.0;
    
        std::complex<T> left_sum = 0.0;
        std::complex<T> right_sum = 0.0;
        
        T from = 0.0;
        T to = 1.0;
        
        for (int i = 0; i < N_segments; ++i) {
            std::complex<T> val = 0.0;
            Poly<T> p_x, p_x_der, p_y, den, num;
            Poly<T> temp;

            den = coefs[i].second;            
            p_x = coefs[i].first[0];
            p_y = coefs[i].first[1];
            p_x_der = coefs[i].first[0].der();
            
            std::vector<std::pair<int, std::complex<T>>> roots = den.solve(type);
            int n_roots = roots.size();
            std::vector<std::complex<T>> vec_root(n_roots);
            for (int j = 0; j < n_roots; ++j) {
                vec_root[j] = roots[j].second;
            }
            std::complex<T> temp_sum = 0.0;
                  
            auto div1 = p_y / den;
            auto div2 = p_x_der / den;
            
            Poly<T> P_12 = div1[0] * div2[0];
            temp_sum +=  P_12.integral(from, to) ;
            num = div1[0] * div2[1] + div2[0] * div1[1];
            auto div3 = num / den;    
            temp_sum += div3[0].integral(from, to) ;
            temp_sum += integral(div3[1], vec_root, from, to) ;
            temp_sum += integral(div1[1], div2[1], vec_root, from, to) ;
            left_sum += temp_sum;

            temp_sum = 0.0;
            div2 =  p_x / den;
            P_12 = div1[0] * div2[0];
            Poly<std::complex<T>> P_12_complex  = make_complex_poly<T>(P_12);
            num = div1[0] * div2[1] + div2[0] * div1[1];
            auto div4 = num / den;
            Poly<std::complex<T>> P4 = make_complex_poly<T>(div4[0]);
            for (int j = 0; j < n_roots; ++j) {
                Poly<std::complex<T>> divider(std::vector<std::complex<T>>{1.0, -vec_root[j]});
                auto div3_complex = P_12_complex / divider;    
                temp_sum +=  integral_type_1(div3_complex[1][0], vec_root[j], 1, from, to) ;
                temp_sum += div3_complex[0].integral(from, to) ;
                auto div5 = P4 / divider;
                temp_sum += div5[0].integral(from, to) ;
                temp_sum += integral_type_1(div5[1][0], vec_root[j], 1, from, to) ;
                temp_sum += integral(div4[1], vec_root, j, from, to) ;         
                temp_sum += integral(div1[1], div2[1], vec_root, j, from, to) ;   
            }
            right_sum += temp_sum;
        }
        sum = left_sum - right_sum;
        return sum;
    }

    //Аналитический интеграл (Алгоритм 1)
    // type -- метод факторизации знаменателя Q(t)
    //      1 - с помощью вспомогательной матрицы
    //      2 - метод Лаггера 
    std::complex<T> analytic_integral2(int type = 1) {
        T from, to;
        std::complex<T> sum = 0.0;
        std::complex<T> left_sum = 0.0;
        std::complex<T> right_sum = 0.0;

        std::vector<std::pair<int, std::complex<T>>> roots;
        std::complex<T> val = 0.0;
        Poly<T> den, num1, num2;
        Poly<T> temp;
        for (int i=0; i < N_segments; ++i) {
            from = knots[ks[i]];
            to = knots[ks[i+1]];

            den = coefs[i].second;
            num1 = coefs[i].first[1] * (coefs[i].first[0].der() * den);
            num2 = - coefs[i].first[1] *  den.der() * coefs[i].first[0];

            std::complex<T> a0 = den.normalize();

            roots = den.solve(type);
            den *= den * den;
            a0 *= a0 * a0;

            for (auto &root: roots) root.first *= 3;

            std::complex<T> b1 = num1.normalize();
            std::complex<T> b2 = num2.normalize();

            if (num1.deg() >= den.deg()) {
                auto div = num1 / den;
                left_sum +=  b1 / a0 * div[0].integral(from, to);
                num1 = div[1];
            }

            if (num2.deg() >= den.deg()) {
                auto div = num2 / den;
                right_sum +=  b2 / a0 * div[0].integral(from, to);
                num2 = div[1];
            }


            std::vector<std::complex<T>> res1 = frac_decomp_matrix(num1, den, roots);
            std::vector<std::complex<T>> res2 = frac_decomp_matrix(num2, den, roots);

            int j = 0;
            for (auto &root: roots) {
                for (int k = 1; k <= root.first; ++k) {
                    left_sum  +=   b1 / a0 * integral_type_1(res1[j], root.second, k, from, to);
                    right_sum +=   b2 / a0 * integral_type_1(res2[j], root.second, k, from, to);
                    ++j;
                }
            }
        }
        sum = left_sum + right_sum;
        return sum;
    }
};
