#pragma once

#include <vector>
#include <cmath>
#include <iostream>
#include <initializer_list>
#include <iomanip>
#include <utility>
#include <complex>

#include "tools.h"
#include "Eigen/Dense"
#include "solver.h"

template < typename Iterator >
std::ostream & PrintIter(std::ostream & os, Iterator begin, Iterator end) {
    for (auto it = begin; it != end; os << * it++)
        if (it != begin) os << ", ";
    return os;
}

template < class T >
std::ostream & PrintVec(std::ostream & os,
                        const std::vector<T> & vec) {
    os << "[";
    return PrintIter(os, vec.begin(), vec.end()) << "]";
}


// Класс многочлен.
// ПРИМЕР: x^3 + 2x^2 + x + 5 := [1, 2, 1, 5]
template < typename T = long double >
class Poly {
    std::vector<T> coef_;
    size_t deg_;


public:
    T EPS = 1e-15;
    Poly();
    Poly(const Poly & p);
    Poly(const std::initializer_list<T>& coefs);
    explicit Poly(size_t deg, T value);
    explicit Poly(std::vector<T> coef);
    explicit Poly(T value);
    size_t deg() const;
    const std::vector<T> & GetCoef() const;
    void Clear();
    Poly der() const;
    Poly der(int deg) const;
    Poly antider() const;
    T integral(T from, T to) const;
    void update();
    void update_real(T eps_ = 1e-10);
    void resize(int deg);
    T normalize();
    T eval(const T& x) const ;
    std::complex<T> At_complex(const std::complex<T> & x) const;
    T At(const T & x) const;
    T At2(const T & x) const; //fma only for double
    std::vector<T> At(const std::vector<T>& xs) const;
    T operator[](size_t index) const;
    T & operator[](size_t index);
    Poly & operator = (const Poly & rhs);
    Poly & operator *= (const Poly & rhs);
    Poly & operator *= (T rhs);
    Poly & operator += (const Poly & rhs);
    Poly & operator += (const T & rhs);
    Poly & operator -= (const Poly & rhs);
    Poly<T> & operator /= (const Poly<T>& poly);
    Poly & operator /= (T value);
    std::vector<std::pair<int, std::complex<T>>> solve(int type = 1);

};

template < typename T >
Poly<T> operator - (T value, Poly<T> p) { 
    p *= T(-1.0);
    p[p.deg()] += value;
	return p;
}

template < typename T >
Poly<T> operator - (Poly<T> p, T value) { 
	p[p.deg()] -= value;
	return p;
}


template < typename T >
Poly<T> operator + (Poly<T> p, T value) { 
	p[p.deg()] += value;
	return p;
}

template < typename T >
Poly<T> operator + (T value, Poly<T> p) { 
	p[p.deg()] += value;
	return p;
}

template < typename T >
Poly<T> operator + (const Poly<T> & p1,
                    const Poly<T> & p2);

template < typename T >
Poly<T> operator * (const Poly<T> & p1,
                    const Poly<T> & p2);

template < typename T >
Poly<T> operator * (const Poly<T> & p, T value);

template < typename T >
Poly<T> operator / (const Poly<T> & p,
                    const T & value);

template < typename T >
std::vector<Poly<T>> operator /(const Poly<T> & lhs,
                                const Poly<T> & rhs);


template < typename T >
std::ostream & operator << (std::ostream & out,
                            const Poly<T> & p);

template < typename T >
bool operator == (const Poly<T> & p1,
                  const Poly<T> & p2);

/////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////

// Конструкторы

template < typename T >
Poly<T> ::Poly(): coef_({0}), deg_(0) {}

template < typename T >
Poly<T> ::Poly(const Poly<T> & p): coef_(p.GetCoef()), deg_(p.deg()) {}

template < typename T >
Poly<T> ::Poly(std::vector<T> coef) {
    coef_ = std::move(coef);
    deg_ = (coef_.size() - 1);
}

template < typename T >
Poly<T> ::Poly(T value) {
    coef_ = std::vector<T> (1, value);
    deg_ = 0;
}

template < typename T >
Poly<T> ::Poly(size_t deg, T value) {
    coef_ = std::vector<T> (deg+1, value);
    deg_ = deg;
}

template < typename T >
T Poly<T>::normalize() {
    T a0 = coef_[0];
    if (a0 != 0.)
        for (T& c : coef_) {
            c /= a0;
        }
    return a0;
}

template < typename T >
size_t Poly<T> ::deg() const {
    return deg_;
}

template < typename T >
const std::vector<T> & Poly<T> ::GetCoef() const {
    return coef_;
}

template < typename T >
void Poly<T> ::Clear() {
    coef_.clear();
    coef_.push_back((T) 0);
    deg_ = 0;
}

template <typename T>
T Poly<T>::integral(T from, T to) const {
    Poly<T> antider = this->antider();
    return antider.At(to) - antider.At(from);
}


// Evaluating antiderivative F(x) : F'(x) = f(x)
// Вычисление первообразной F(x) : F'(x) = f(x)

template <typename T>
Poly<T> Poly<T>::antider() const {
    // x^2 [1, 0, 0] -> x^/3 [1/3, 0, 0, 0]
    if (this->deg() == 0)
        return Poly<T>(std::vector<T>({coef_[0], 0}));

    std::vector<T> tmp(deg_+2, 0);
    for (size_t i = 0; i <= deg_; ++i)
        tmp[i] = coef_[i] / T(deg_-i+1);

    return Poly<T>(tmp);

}

// Evaluating derivative
template <typename T>
Poly<T> Poly<T>::der() const {
    /// from x + x^2 [1 1 0] -> 2x + 1 [2 1]
    if (this->deg() == 0)
        return Poly<T>({0});

    if (this->deg() == 1)
        return Poly<T>({coef_[0]});

    std::vector<T> tmp(deg_, 0);
    for (size_t i = 0; i < deg_; ++i)
        tmp[i] = coef_[i] * T(deg_-i);

    return Poly<T>(tmp);
}

//TODO


// Help functions for calculating a polynomial at the point

template < typename T >
T two_sum(T & t, T a, T b) {
    T s = a + b;
    T bs = s - a;
    T as = s - bs;
    t = (b - bs) + (a - as);
    return s;
}

template < typename T >
T two_prod(T & t, T a, T b) {
    T p = a * b;
    t = std::fma(a, b, -p); // t = a*b-p
    return p;
}

template < typename T >
inline T Poly<T> ::At2(const T & x) const{
    T s = 0.0, c = 0.0, p, pi, t;
    for (const T & a: coef_) {
        p = two_prod(pi, s, x);
        s = two_sum(t, p, a);
        c *= x;
        c += pi + t;
    }
    return s + c;
}

template < typename T >
std::vector<T> Poly<T>::At(const std::vector<T>& xs) const{
    int n = xs.size();
    std::vector<T> res(xs.size());
    for (int i = 0; i < n;++i) {
        res[i] = this->At(xs[i]);
    }
    return res;
}

template < typename T >
inline T Poly<T> ::At(const T & x) const{
    T b = coef_[0];
    for (int d = 1; d <= deg_; ++d) {
        b = coef_[d] + b * x;
    }
    return b;
}

template < typename T >
std::complex<T> Poly<T> ::At_complex(const std::complex<T> & x) const{
    Poly<std::complex<T>> poly_complex = make_complex_poly(*this);
    std::complex<T> b = poly_complex[0];
    for (int d = 1; d <= deg_; ++d) {
        b = poly_complex[d] + b * x;
    }
    return b;
}

template < typename T >
T Poly<T> ::operator[](size_t index) const {
    return coef_[index];
}

template < typename T >
T & Poly<T> ::operator[](size_t index) {
    return coef_[index];
}

template < typename T >
Poly<T> & Poly<T> ::operator = (const Poly<T> & rhs) {
    coef_ = rhs.GetCoef();
    deg_ = rhs.deg();
    return *this;
}

template < typename T >
Poly<T> operator * (const Poly<T> & p1,
                    const Poly<T> & p2) {
    size_t lhs_deg = p1.deg();
    size_t rhs_deg = p2.deg();
    Poly<T> res(rhs_deg + lhs_deg, 0);
    for (size_t i = 0; i <= lhs_deg; ++i) {
        for (size_t j = 0; j <= rhs_deg; ++j) {
            res[i + j] += p1[i] * p2[j];
        }
    }
    return res;
}

template < typename T >
Poly<T> & Poly<T> ::operator *= (const Poly<T> & rhs) {
    *this = *this * rhs;
    return *this;
}

template < typename T >
Poly<T> & Poly<T> ::operator *= (T rhs) {
    for (size_t i = 0; i <= deg_; ++i) {
        this -> coef_[i] *= rhs;
    }
    return *this;
}

template < typename T >
Poly<T> & Poly<T> ::operator += (const T & rhs) {
    coef_[deg_] += rhs;
    return *this;
}

template < typename T >
Poly<T> & Poly<T> ::operator += (const Poly<T> & rhs) {
    Poly temp = *this + rhs;
    *this = temp;
    return *this;
}

template < typename T >
Poly<T> & Poly<T> ::operator /= (T value) {
    for (size_t i = 0; i <= deg_; ++i) {
        coef_[i] = coef_[i] / value;
    }
    return *this;
}

template < typename T >
Poly<T> & Poly<T> ::operator /= (const Poly<T>& poly) {
    *this = (*this / poly)[0];
    return *this;
}

template<typename T>
std::vector<std::pair<int, std::complex<T>>> solve_quadratic(const Poly<T>& poly){
    T a = poly[0];
    T b = poly[1];
    T c = poly[2];
    T D = b*b-4*a*c;
    if ( std::abs(D) < poly.EPS)
        return {{2, std::complex<T>(-b/2/a)}};
    
    if (D > 0) 
        return {{1, std::complex<T>((-b + std::sqrt(D))/2.0/a)},
                {1, std::complex<T>((-b - std::sqrt(D))/2.0/a)}};

    return {{1, std::complex<T>(-b/2.0/a, std::sqrt(-D)/2.0/a)}, 
            {1, std::complex<T>(-b/2.0/a, -std::sqrt(-D)/2.0/a)}};
}

template<typename T>
std::vector<std::pair<int, std::complex<T>>> solve_cubic2(const Poly<T>& poly){
    // solve cubic poly
    T a = poly[1] / poly[0];
    T b = poly[2] / poly[0];
    T c = poly[3] / poly[0];
    T q = (a * a - 3.0 * b);
    T r = (2.0 * a * a * a - 9.0 * a * b + 27.0 * c);
    T Q = q / 9.0;
    T R = r / 54.0;
    T Q3 = Q * Q * Q;
    T R2 = R * R;

    T CR2 = 729 * r * r;
    T CQ3 = 2916 * q * q * q;

    T sqrtQ = std::sqrt(Q);
    
    if ( std::abs(R) < 1e-12 && std::abs(Q) < 1e-12) {
        std::complex<T> root(-a / 3.0);
        return {{3, root}};

    } else if (std::abs(CR2 - CQ3) < 1e-12) {
        T sqrtQ = std::sqrt(Q); 
        if (R > 0) {
            std::complex<T> root1(-2 * sqrtQ  - a / 3);
            std::complex<T> root2(sqrtQ - a / 3);
            return {{1, root1}, {2,root2}};
        } else {
            std::complex<T> root1(- sqrtQ  - a / 3);
            std::complex<T> root2(2 * sqrtQ - a / 3);
            return {{2, root1}, {1, root2}};
        }
    } else if (R2 < Q3) {
        T sgnR = (R >= 0 ? 1.0 : -1.0);
        T ratio = sgnR * std::sqrt(R2 / Q3);
        T theta = std::acos(ratio);
        T norm = -2 * std::sqrt(Q);
        T root1 = norm * std::cos(theta / 3) - a / 3;
        T root2 = norm * std::cos((theta + 2.0 * M_PI) / 3) - a / 3;
        T root3 = norm * std::cos((theta - 2.0 * M_PI) / 3) - a / 3;
        return {{1, std::complex<T>(root1)}, {1, std::complex<T>(root2)}, {1, std::complex<T>(root3)}};
    } else {
        T sgnR = (R >= 0 ? 1.0 : -1.0);
        T A = -sgnR * std::pow ( std::abs(R) + std::sqrt(R2 - Q3), 1.0/3.0);
        T B = Q / A ;
        T root = (A + B - a / 3.0);
        Poly<T> quadratic = (poly / Poly<T>({1.0, -root}))[0];
        auto res = solve_quadratic(quadratic);
        res.push_back({1, std::complex<T>(root)});
        return res;
    }
}

template<typename T>
std::vector<std::pair<int, std::complex<T>>>  solve_cubic(Poly<T>& poly)  {
    std::vector<std::pair<int, std::complex<T>>> roots;
    T x[3];
    int type = solveP3(x, poly[1]/poly[0], poly[2]/poly[0],  poly[3]/poly[0]);
    switch (type) {
    case 1:
        roots.push_back({1, x[0]});
        roots.push_back({1, std::complex<T>(x[1], x[2])});
        roots.push_back({1, std::complex<T>(x[1], -x[2])});
        break;
    default:
        roots.push_back({1, x[0]});
        roots.push_back({1, x[1]});
        roots.push_back({1, x[2]});
        break;
    }
    return roots;
}

template<typename T>
std::vector<std::pair<int, std::complex<T>>>  solve_quartic(Poly<T>& poly)  {
    std::vector<std::pair<int, std::complex<T>>> roots(4);
    std::vector<std::complex<T>> vec_root = solveP4(poly[1]/poly[0], 
                                                    poly[2]/poly[0], 
                                                    poly[3]/poly[0], 
                                                    poly[4]/poly[0]);
    for (int i = 0; i < 4; ++i) {
        roots[i] = {1, vec_root[i]};
    }
    return roots;
}


template<typename T>
std::vector<std::pair<int, std::complex<T>>>  Poly<T>::solve(int type)  {
    switch (deg_) {
    case 1:
        return {{1, std::complex<T>(-coef_[1]/coef_[0])}};
    case 2:
        return solve_quadratic(*this);
    // case 3:
    //     return solve_cubic(*this);
    // case 4:
    //     return solve_quartic(*this);
    default:
        return solve_numerical(*this, type);
        break;
    }
    return {{1, 0}};
}



template<typename T>
T upper_bound(const Poly<T>& poly) {
    T upper_bound = std::abs(poly[1] / poly[0]);
    for(int i = 2; i <= poly.deg(); ++i)
        if (std::abs(poly[i] / poly[0]) > upper_bound)
            upper_bound = std::abs(poly[i] / poly[0]);

    return 1 + upper_bound;
}

template<typename T>
T lower_bound(const Poly<T>& poly) {
    return  1 / upper_bound(poly);
}


template<typename T>
T upper_bound(const Poly<std::complex<T>>& poly) {
    T upper_bound = std::abs(poly[1] / poly[0]);
    for(int i = 2; i <= poly.deg(); ++i)
        if (std::abs(poly[i] / poly[0]) > upper_bound)
            upper_bound = std::abs(poly[i] / poly[0]);

    return 1 + upper_bound;
}

template<typename T>
T lower_bound(const Poly<std::complex<T>>& poly) {
    return  1 / upper_bound(poly);
}

// NEED normalize before using !
// Optimization
template<typename T>
std::vector<std::pair<int, std::complex<T>>> solve_with_eigenvalues(Poly<T> poly) {
    if (poly[0] != 1.0)
        poly.normalize();

    std::vector<std::pair<int, std::complex<T>>> res;
    size_t n = poly.deg();
    Eigen::MatrixX<T> m(n,n);
    
    m.setZero();
    for(int i=1; i <= n; ++i){
        m(i-1,n-1) = -poly[n-i+1];
    }
    for(int i=0; i < n-1; ++i)
        m(i+1,i) = 1.0;
    
    auto roots =  m.eigenvalues();
    for (int i = 0; i < roots.size(); ++i) {
        res.push_back({1, std::abs(roots[i]) < poly.EPS ? 0.0 : roots[i]});
    }
    return res;
}

template<typename T>
std::vector<std::pair<int, std::complex<T>>> solve_laguerre(const Poly<T>& poly) {
    std::vector<std::pair<int, std::complex<T>>> res;
    int deg = poly.deg();
    T eps = 1e-12;
    int k = 0;
    while(std::abs(poly[deg-k]) < 1e-16) ++k;
    if (k != 0) res.push_back({k, 0.0});

    Poly<std::complex<T>> p(deg - k, 0);
    deg = p.deg();

    for (int i = 0; i <= deg; ++i) {
        p[i] = std::complex<T>{poly[i]};
    }

    std::complex<T> G, H, d1, d2, a;
    Poly<std::complex<T>> X{1,0};
    int m;
    int d = 0;

    while(d < deg){
        m = 1;
        std::complex<T> x = lower_bound(p);
        Poly<std::complex<T>> p1 = p.der();
        Poly<std::complex<T>> p11 = p1.der();
        for (k = 0; k < 1000; ++k) {
            if (std::abs(p.eval(x)) < eps) break;
            G = p1.eval(x) / p.eval(x);
            H = G * G - p11.eval(x) / p.eval(x);
            d1 = G + std::sqrt(std::complex<T>(deg - 1) * (std::complex<T>(deg) * H - G * G));
            d2 = G - std::sqrt(std::complex<T>(deg - 1) * (std::complex<T>(deg) * H - G * G));
            a = std::complex<T>(deg) / (std::abs(d1) > std::abs(d2) ? d1 : d2);
            if (std::abs(a) < eps) break;
            x -= a;
        }
        polish_complex(x, T(1e-14));

        if (std::abs(x.imag()) >= T(1e-9)) {
            ++d;
            X[1] = -std::conj(x);
            p /= X;
            res.push_back({m, std::conj(x)});         
        }      

        X[1] = -x;
        p /= X; ++d;
        res.push_back({m, x});
        X[1] = 0;
    }
    return res;
}

template<typename T>
std::vector<std::pair<int, std::complex<T>>> solve_numerical(const Poly<T>& poly, int type=1) {
    if (type == 1) {
        return solve_with_eigenvalues(poly);
    } else {
        return solve_laguerre(poly);
    }
}

template<typename T>
void Poly<T>::update() {
    size_t i = 0;
    while(coef_[i] == T(0) && i < deg_) {
        ++i;
    }
    if (i == deg_+1) {
        deg_ = 0;
        coef_ = std::vector<T>{0};
    } else {
        coef_ = std::vector<T>(coef_.begin() + i, coef_.end());
        deg_ = deg_ - i;
    }
}

template<typename T>
void Poly<T>::update_real(T eps_) {
    size_t i = 0;
    while(std::abs(coef_[i]) < eps_ && i < deg_) {
        ++i;
    }
    if (i == deg_+1) {
        deg_ = 0;
        coef_ = std::vector<T>{0};
    } else {
        coef_ = std::vector<T>(coef_.begin() + i, coef_.end());
        deg_ = deg_ - i;
    }
}


template <typename T>
void Poly<T>::resize(int deg) {
    int temp = deg_;
    deg_ = deg;
    coef_ = std::vector<T>(coef_.begin() + (temp - deg),  coef_.end());
}

template<typename T>
Poly<T>::Poly(const std::initializer_list<T> &coefs) {
    coef_ = coefs;
    deg_ = coef_.size() - 1;
}

template < typename T >
inline T Poly<T> ::eval(const T & x) const{
    T b = coef_[0];
    for (int d = 1; d <= deg_; ++d) {
        b = coef_[d] + b * x;
    }
    return b;
}

template<typename T>
Poly<T> &Poly<T>::operator-=(const Poly &rhs) {
    *this = *this - rhs;
    return *this;
}

template<typename T>
Poly<T> operator-(const Poly<T>& poly) {
    size_t deg = poly.deg();
    Poly<T> temp = poly;
    for (int i = 0; i <= deg; ++i){
        temp[i] = -poly[i];
    }
    return temp;
}


template < typename T >
Poly<T> plus (const Poly<T> & lhs,
              const Poly<T> & rhs) {
    size_t deg_lhs = lhs.deg();
    size_t deg_rhs = rhs.deg();

    Poly<T> temp = lhs;
    for (size_t i = 0; i <= deg_rhs; ++i){
        temp[deg_lhs -i ] += rhs[deg_rhs - i];
    }
    return temp;
}

template < typename T >
Poly<T> operator + (const Poly<T> & p1,
                    const Poly<T> & p2) {
    // [1,1,1,1] + [1,1] = [1,1,2,2]
    return p1.deg() >= p2.deg() ? plus(p1,p2) : plus(p2, p1);
}


template < typename T >
Poly<T> operator - (const Poly<T> & lhs,
                    const Poly<T> & rhs) {

    return lhs + (-rhs);
}


template < typename T >
Poly<T> operator * (const Poly<T> & p, T value) {
    if (value == 0)
        return Poly<T>({0});
    size_t deg = p.deg();
    Poly<T> ans(deg, 0);
    for (size_t i = 0; i <= deg; ++i) {
        ans[i] = p[i] * value;
    }
    return ans;
}

template < typename T >
Poly<T> operator * (T value, const Poly<T> & p) {
    return p * value;
}

template < typename T >
Poly<T> operator / (const Poly<T> & p,
                    const T & value) {
    size_t deg = p.deg();
    Poly<T> ans(deg, 0);
    T temp = 1.0 / value;
    for (size_t i = 0; i <= deg; ++i) {
        ans[i] = p[i] * temp;
    }
    return ans;
}

template<typename T>
std::vector<Poly<T>> operator/(const Poly<T>& lhs,
                                const Poly<T>& rhs) {
    // P(t) = rhs(t) * res(t) + R(t), degR < deg(rhs(t))
    // return [res(t), R(t)], where R - remainder poly
    size_t deg1 = lhs.deg();
    size_t deg2 = rhs.deg();

    if (deg1 < deg2) {
        return {Poly<T>{std::vector<T>(1, 0.0)}, lhs};
    }

    Poly<T> rem = lhs;
    std::vector<T> res;
    T a_i = 0;
    for (int i = 0; i <= deg1 - deg2  ; ++i) {
        a_i = rem[i];
        rem[i] = 0;
        res.push_back(a_i / rhs[0]);
        for (int j = 1; j <= deg2; ++j) {
            rem[i + j] -= rhs[j]*a_i / rhs[0];
        }
    }
    rem.update();
    return {Poly<T>{res}, rem};
}

template < typename T >
std::ostream & operator << (std::ostream & out,
                            const std::complex<T> & val) {
    return out << "(" << val.real() << " + (" << val.imag() << "i))";
}


template < typename T >
std::ostream & operator << (std::ostream & out,
                            const Poly<T> & p) {
    int deg = p.deg();
    
    if (deg == 0) {
        return out << p[0];
    }
    
    out << p[0] << " x^" << deg;
    
    for (int i = 1; i < deg; ++i) { 
        out << " + " << p[i] << "x^" << deg - i;
    }
    out << " + " << p[deg];
    return out;
}


template < typename T >
bool operator == (const Poly<T> & p1, Poly<T> & p2) {
    return p1.GetCoef() == p2.GetCoef();
}

template<typename T>
Poly<T> accumulate(std::vector<Poly<T>> polys){
    if (polys.size() == 0) return Poly<T>{0};

    Poly<T> res_poly = polys[0];
    for (int i = 1; i < polys.size(); ++i)
        res_poly  *= polys[i];
    return res_poly;
}

template<typename T>
Poly<std::complex<T>> make_complex_poly(const Poly<T>& poly) {
    int deg = poly.deg();
    Poly<std::complex<T>> p(deg, 0);
    for (int i = 0; i <= deg; ++i)
        p[i] = std::complex<T>{poly[i]};
    return p;
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

    Eigen::VectorX<std::complex<T>> x = matrix.colPivHouseholderQr().solve(b);
    for(int i =0; i < x.size(); ++i)
        res.push_back(x(i));

    return res;
}
