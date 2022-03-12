#pragma once
#include "../include/Point.h"
#include "../include/Poly.h"
#include "../include/tools.h"
#include "../include/integral.h"

template<typename T>
Rational<T> de_boor_nurbs(int k,
                          std::vector<T>& knots,
                          std::vector<T>& weights,
                          std::vector<Point<T>>& points,
                          int p) {

    size_t dim = points[0].dim;
    std::vector<Point<Poly<T>>> d((size_t)p+1, Point<Poly<T>>(dim));
    std::vector<Poly<T>> d_den((size_t)p+1);

    for (int i = 0; i < p+1; ++i){
        d_den[i][0] = weights[i + k - p];
        for (int j = 0; j < dim; ++j)
            d[i][j][0] = points[i + k - p][j] * weights[i + k - p];
    }

    T B,D, den;
    for (int r = 1; r < p+1; ++r) {
        for (int j=p ; j > r-1; --j) {
            den = knots[j+1+k-r] - knots[j+k-p];
            B = - knots[j+k-p];
            D = knots[j+1+k-r];
            Poly<T> temp1_den = d_den[j] * Poly<T>{1,0} + B * d_den[j];
            Poly<T> temp2_den = d_den[j-1] * Poly<T>{-1,0} + D * d_den[j-1];
            d_den[j] = (temp1_den + temp2_den) / den;
            for (size_t i = 0; i < dim; ++i){
                Poly<T> temp1 = d[j][i] * Poly<T>{1,0} + B * d[j][i] ;
                Poly<T> temp2 = d[j-1][i] * Poly<T>{-1,0} + D * d[j-1][i];
                d[j][i] = (temp1 + temp2) / den;
            }
        }
    }
    return {d[p], d_den[p]};
}


template <typename T>
std::vector<std::complex<T>> frac_decomp_matrix(Poly<T> num, Poly<T> den,
                                                std::vector<std::pair<int, std::complex<T>>>  roots) {
    // Ax = b
    size_t n = roots.size();
    std::vector<std::complex<T>> res;

    size_t num_deg = num.deg();
    size_t den_deg = den.deg();

    Poly<std::complex<T>> Q(den.deg(), std::complex<T>(0.0));
    Poly<std::complex<T>> P(num.deg(), std::complex<T>(0.0));

    for (int i = 0; i <= num_deg; ++i)
        P[i] = std::complex<T>(num[i]);

    for (int i = 0; i <= den_deg;++i)
        Q[i] = std::complex<T>(den[i]);

    Eigen::VectorX<std::complex<T>> b(Q.deg());
    b.setZero();
    for(int i =0; i <= P.deg(); ++i)
        b(i) = P[P.deg()-i];

    Eigen::MatrixX<std::complex<T>> matrix(Q.deg(), Q.deg());
    matrix.setZero();

    int ii = 0;
    for (int i = 0; i < n; ++i) {
        Poly<std::complex<T>> Q_temp = Q;
        Poly<std::complex<T>> x_a {1, -roots[i].second};

        for (int j=0; j < roots[i].first; ++j) {
            Q_temp /= x_a;
            size_t deg = Q_temp.deg();
            for (int k=0; k <= deg; ++k) {
                matrix(k, ii) = Q_temp[deg-k];
            }
            ++ii;
        }
    }
    //std::cout << matrix << "\n";
    Eigen::VectorX<std::complex<T>> x = matrix.colPivHouseholderQr().solve(b);

    //std::cout << x << "\n";
    for(int i =0; i < x.size(); ++i)
        res.push_back(x(i));

    return res;
}

template<typename T>
using Rational = std::pair<Point<Poly<T>>, Poly<T>>;


template<typename T>
struct NURBS{
    int p;
    int dim;

    std::vector<T> knots;
    std::vector<T> weights;
    std::vector<Point<T>> cv;
    std::pair<int,int> domain;
    std::vector<int> ks;
    std::vector<Rational<T>> coefs;

public:
    NURBS(int p_,
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
        create_coefs();
        polishing();
    }

    NURBS(int p_,
          std::vector<T>& weights_,
          std::vector<Point<T>>& cv_) {
        p = p_;
        cv = cv_;
        weights = weights_;
        knots = create_knots<T>(cv.size(), p);
        dim = cv[0].dim;
        domain = {p, knots.size() - p - 1};
        ks = create_intervals(domain,knots);
        create_coefs();
        polishing_uniform();
    }

    void create_coefs(){
        for (int k : ks) coefs.push_back(de_boor_nurbs( k, knots,weights, cv, p));
    }

    void polishing() {
        for (auto& frac: coefs){
            frac.second.update_real();
            for (int d=0; d < dim; ++d){
                frac.first[d].update_real();
            }
        }
    }

    void polishing_uniform() {
        int len = coefs.size();
        if (len <= 2*(p-1)) {
            polishing();
        } else {
            for (int i = 0; i < p-1; ++i){
                coefs[i].second.update_real();
                coefs[len - i - 1].second.update_real();

                for (int d=0; d < dim; ++d){
                    coefs[i].first[d].update_real();
                    coefs[len-i-1].first[d].update_real();
                }
            }
            for (int i = p-1; i < len-p+1; ++i) {
                coefs[i].second.resize(1);
                for (int d=0; d < dim; ++d){
                    coefs[i].first[d].update_real();
                }
            }
        } 
    }


    std::vector<Point<T>> get_points(int N) {
        std::vector<Point<T>> points(N, Point<T>(dim));
        std::vector<T> ts = linspace(knots[domain.first], knots[domain.second], N);
        T t = knots[domain.first];
        int i = 0;
        int j = 0;
        for(int k : ks){
            while (t <= knots[k+1] && i < N){
                for (int d = 0; d < dim; ++d){
                    points[i][d] = coefs[j].first[d].At(t) / coefs[j].second.At(t) ;
                }
                ++i; t = ts[i];
            }
            ++j;
        }
        return  points;
    }

    std::vector<T> get_dy_dx_points(int N) {
        std::vector<T> points(N, 0.0);
        std::vector<T> ts = linspace(knots[domain.first],
                                     knots[domain.second], N);
        T t = knots[domain.first];
        int i = 0;
        int j = 0;
        for (int k: ks) {
            Poly<T> p_x = coefs[j].first[0];
            Poly<T> p_y = coefs[j].first[1];
            Poly<T> den = coefs[j].second;

            while (t <= knots[k + 1] && i < N) {
                points[i] = (p_y.der().At(t) * den.At(t) - p_y.At(t) *den.der().At(t)) /
                            (p_x.der().At(t) * den.At(t) - p_x.At(t) *den.der().At(t));

                ++i;
                t = ts[i];
            }
            ++j;
        }
        return points;
    }

    void save_coefs(std::string path_dir = "") {
        std::string path_x="coefs_num_x.out";
        std::string path_y="coefs_num_y.out";
        std::string den = "coefs_den.out";

        std::ofstream out_x(path_dir + path_x);
        std::ofstream out_y(path_dir + path_y);
        std::ofstream out_den(path_dir + den);

        for (const Rational<T>& coef : coefs){
            out_x << coef.first[0] << "\n";
            out_y << coef.first[1] << "\n";
            out_den << coef.second << "\n";
        }
    }

    T numerical_integral() {
        auto ps = POINTS;
        auto ws = WEIGHTS;
        int n = ps.size();
        T A, B;
        T sum = 0.0;
        int i = 0;
        for(auto k : ks) {
            A = knots[k];
            B = knots[k+1];
            Poly<T> p_x = coefs[i].first[0];
            Poly<T> p_y = coefs[i].first[1];
            Poly<T> den = coefs[i].second;
            T sum_temp = 0.0;
            for (int j = 0; j < n; ++j) {
                T x = T(B-A)*ps[j]/2 + T(A+B)/2;
                T den_quadratic = std::pow(den.At(x),2);
                sum_temp += ws[j] * p_y.At(x) * p_x.der().At(x) / den_quadratic;
                sum_temp -= ws[j] * p_y.At(x) * den.der().At(x)*p_x.At(x) / den_quadratic / den.At(x); //
            }
            sum += (B-A)/2.0 *sum_temp;
            ++i;
        }
        return sum;
    }

    std::complex<T> analytic_integral1(int type = 1) {
        std::complex<T> sum = 0.0;
        size_t n = ks.size();
        std::complex<T> left_sum = 0.0;
        std::complex<T> right_sum = 0.0;
        
        for (int i = 0; i<n; ++i) {
            std::complex<T> val = 0.0;
            Poly<T> p_x, p_x_der, p_y, den, num;
            Poly<T> temp;
            T from = knots[ks[i]];
            T to = knots[ks[i] + 1];
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
            
            std::complex<T> norm_den = den.normalize();
            
            norm_den *= norm_den;         
            auto div1 = p_y / den;
            auto div2 = p_x_der / den;
            
            Poly<T> P_12 = div1[0] * div2[0];
            temp_sum +=  P_12.integral(from, to) / norm_den;
            num = div1[0] * div2[1] + div2[0] * div1[1];
            auto div3 = num / den;    
            temp_sum += div3[0].integral(from, to) / norm_den;
            temp_sum += integral(div3[1], vec_root, from, to) / norm_den;
            temp_sum += integral(div1[1], div2[1], vec_root, from, to) / norm_den;
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
                temp_sum +=  integral_type_1(div3_complex[1][0], vec_root[j], 1, from, to) / norm_den;
                temp_sum += div3_complex[0].integral(from, to) / norm_den;
                auto div5 = P4 / divider;
                temp_sum += div5[0].integral(from, to) / norm_den;
                temp_sum += integral_type_1(div5[1][0], vec_root[j], 1, from, to) / norm_den;
                temp_sum += integral(div4[1], vec_root, j, from, to) / norm_den;         
                temp_sum += integral(div1[1], div2[1], vec_root, j, from, to) / norm_den;   
            }
            right_sum += temp_sum;
        }
        sum = left_sum - right_sum;
        return sum;
    }


    //using second method
    std::complex<T> analytic_integral2(int type = 1) {
        T from, to;
        std::complex<T> sum = 0.0;
        std::complex<T> left_sum = 0.0;
        std::complex<T> right_sum = 0.0;

        size_t n = ks.size();
        std::vector<std::pair<int, std::complex<T>>> roots;
        std::complex<T> val = 0.0;
        Poly<T> den, num1, num2;
        Poly<T> temp;
        for (int i=0; i<n; ++i) {
            from = knots[ks[i]];
            to = knots[ks[i] + 1];

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
