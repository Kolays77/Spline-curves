#pragma once
#include "../include/Point.h"
#include "../include/Poly.h"
#include "../include/tools.h"
#include "../include/integral.h"


template<typename T>
using Rational = std::pair<Point<Poly<T>>, Poly<T>>;

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
struct NURBS{
    int p;
    int dim;

    std::vector<T> knots;
    std::vector<T> weights;
    std::vector<Point<T>> cv;
    std::pair<int,int> domain;
    std::vector<int> ks;
    std::vector<Rational<T>> coefs;
    std::string type;

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
        type = "default";
        create_coefs();

    }

    NURBS(int p_,
          std::vector<T>& weights_,
          std::vector<Point<T>>& cv_,
          std::string type_="uniform2") {
        p = p_;
        if (type_ == "uniform2") {
            knots = create_knots<T>(cv_.size(), p);
            type = type_;
        }
        weights = weights_;
        cv = cv_;
        dim = cv[0].dim;
        domain = {p, knots.size() - p - 1};
        ks = create_intervals(domain,knots);
        create_coefs();
    }

    void create_coefs(){
        for (int k : ks) 
            coefs.push_back(de_boor_nurbs( k, knots,weights, cv, p));

        if (type == "uniform2") {
            int len = coefs.size();
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
        } else {
            for (auto& frac: coefs){
                frac.second.update_real();
                for (int d=0; d < dim; ++d){
                    frac.first[d].update_real();
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

    T integral() {
        //auto ps  = get_integrate_points(); //  load_vector<T>("integrate_points.txt");
        //auto ws = get_integrate_weights(); // load_vector<T>("integrate_weights.txt");
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
                T den_cubic = std::pow(den.At(x),3);
                sum_temp += ws[j]* p_y.At(x) * (p_x.der().At(x)*den.At(x) - den.der().At(x)*p_x.At(x)) / den_cubic;
            }
            sum += (B-A)/2.0 *sum_temp;
            ++i;
        }
        return sum;
    }


    //using second method
    std::complex<T> analytic_integral(int type = 1) {
        T from, to;
        std::complex<T> sum = 0.0;
        size_t n = ks.size();
        std::vector<std::pair<int, std::complex<T>>> roots;
        std::vector<std::pair<int, std::complex<T>>> roots_num;

        Poly<T> p_x, p_y, den, NUM;
        Poly<T> temp;
        for (int i=0; i<n; ++i) {
            from = knots[ks[i]];
            to = knots[ks[i] + 1];
            p_x = coefs[i].first[0];
            p_y = coefs[i].first[1];
            den = coefs[i].second;
            NUM = p_y * (p_x.der() * den - den.der() * p_x);

            NUM.update_real();

            std::complex<T> a0 = den.normalize();
            roots = den.solve(type);

            den *= den * den;
            a0 *= a0 * a0;

            for (auto &root: roots) {
                root.first *= 3;
            }

            std::complex<T> b0 = NUM.normalize();

            if (NUM.deg() >= den.deg()) {
                auto div = NUM / den;
                sum +=  b0 / a0 * div[0].integral(from, to) ;
                NUM = div[1];
            }

            std::vector<std::complex<T>> res = frac_decomp_matrix(NUM, den, roots);
            int j = 0;
            for (const auto &root: roots) {
                for (int k = 1; k <= root.first; ++k) {
                    sum +=  b0 / a0 * Integral_first_type(res[j], root.second, k, from, to);
                    ++j;
                }
            }
        }
        return sum;
    }

};



