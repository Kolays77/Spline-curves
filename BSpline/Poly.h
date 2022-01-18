#pragma once
#include <vector>
#include <cmath>
#include <iostream>
#include <initializer_list>
#include <iomanip>
#include <utility>
#include <complex>
#include "Eigen/Core"

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

template < typename T = long double >
class Poly {
    // x^3 + 2x^2 + x + 5 := [1, 2, 1, 5]
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
    Poly antider() const;
    T integral(T from, T to) const;

    void update();
    void update_real();

    T normalize();

    T eval(const T& x) const ;

    T At(const T & x) const;
    std::vector<T> At(const std::vector<T> & points);

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
    std::vector<Poly<T>> solve();

};



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

// Constructors

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


// Evaluating antiderivative F'(x) = f(x)

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
inline T Poly<T> ::At(const T & x) const{
    T s = 0.0, c = 0.0, p, pi, t;
    for (const T & a: coef_) {
        p = two_prod(pi, s, x);
        s = two_sum(t, p, a);
        c *= x;
        c += pi + t;
    }
    return s + c;
}

template<typename T>
T Horner(const Poly<T>& poly, T x) {
    int deg = poly.deg();
    T b = poly[0];
    for (int d = 1; d <= deg; ++d) {
        b = poly[d] + b * x;
    }
    return b;
}

template < typename T >
std::vector<T> Poly<T> ::At(const std::vector<T> & points) {
    std::vector<T> values(points.size());
    const size_t len = points.size();
    for (size_t i = 0; i < len; ++i) {
        values[i] = Poly<T> ::At(points[i]);
    }
    return values;
}

template < typename T >
T Poly<T> ::operator[](size_t index) const {
    if (index > deg_)
        return 0;
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
std::vector<Poly<T>> solve_quadratic(const Poly<T>& poly){
    T b = poly[1] / poly[0];
    T c = poly[2] / poly[0];
    T D = b*b - 4*c;
    if (D < 0)
        return {Poly<T>({1, b, c})};

    return {Poly<T>({1, (b-sqrt(D))/2}),
            Poly<T>({1, (b+sqrt(D))/2})};
}

template<typename T>
std::vector<Poly<T>> Poly<T>::solve()  {
    if (deg_ == 1){
        // [2 1] -> [1 1/2]
        return { Poly<T>({ T(1.0), coef_[1]/coef_[0]}) };
    }
    if (deg_ == 2){
        return solve_quadratic(*this);
    }

    if (deg_ == 3){
        return solve_cubic(*this);
    }

    if (deg_ == 4){
        return solve_quartic(*this);
    }
    return { Poly<T>(std::vector<T>{0}) };
}

template<typename T>
inline T cubic_root1(T value){
    return value < 0 ? -pow(-value,1./3.) : pow(value, 1./3.);
}

template<typename T>
std::vector<Poly<T>> solve_cubic(const Poly<T>& poly) {
    std::cout << std::setprecision(16);
    T a = poly[0];
    T b = poly[1];
    T c = poly[2];
    T d = poly[3];
    T p = (3*a*c - b*b) / (3.*a*a);
    T q = (2*b*b*b - 9.*a*b*c + 27.*a*a*d) / (27.*a*a*a);
    T Q = (p/3.)*(p/3.)*(p/3.) + (q/2.)*(q/2.);

    if (Q  < 0 ) { //3 real unique root
        T y1 =   2. * sqrt(-p / 3.0) * cos(1/3.0 * acos(- q / (2*sqrt( - pow((p/3.),3)))));
        T y2 = - 2. * sqrt(-p / 3.0) * cos(1/3.0 * acos(- q / (2*sqrt( - pow((p/3.),3)))) + M_PI/3.);
        T y3 = - 2. * sqrt(-p / 3.0) * cos(1/3.0 * acos(- q / (2*sqrt( - pow((p/3.),3)))) - M_PI/3.);
        return  {Poly<T>{{1,-(y1- b/a/3.0)}},
                 Poly<T>{{1,-(y2- b/a/3.0)}},
                 Poly<T>{{1,-(y3- b/a/3.0)}}};
    }

    T alpha = cubic_root1(-q/2.0 + sqrt(Q));
    T beta = cubic_root1(-q/2.0 - sqrt(Q));
    T t = (alpha + beta);
    Poly<T> real_root({1,- (t - b/a/3.0)});

    if (Q > 0) { // 1 real root, 2 complex
        // z12 = v +- ui; v = - (alpha + beta)/2 ; u = (alpha - beta)/2 sqrt(3)
        T v = - (alpha + beta)/2 - b/a/3;
        T u = sqrt(3)*(alpha - beta)/2;
        Poly<T> temp2({1, -2*v, v*v + u*u});
        return {real_root, temp2};
    }

    if (Q < 1e-9) { // 3 real roots : (1 unique + 1 double) OR (1 triple)
        std::vector<Poly<T>> res = (poly / real_root)[0].solve();
        res.push_back(real_root);
        return res;
    }
    return {Poly<T>(std::vector<T>({0}))};
}

template<typename T>
std::vector<Poly<T>> solve_quartic(const Poly<T>& poly) {
    // method Ferrari
    // x^4 + ax^3 + bx^2 + cx + d = 0
    T a = poly[1] / poly[0];
    T b = poly[2] / poly[0];
    T c = poly[3] / poly[0];
    T d = poly[4] / poly[0];
    Poly<T> resolvent({1, -b, (a*c - 4.0*d), -a*a*d +4.0*b*d - c*c});
    std::vector<Poly<T>> resolvent_roots = resolvent.solve();
    T y1 = 0.0;
    T temp_y1 = 0.0;
    Poly<T> right_part;
    Poly<T> right_solve_temp;

    for (const Poly<T>& root : resolvent_roots){
        if (root.deg() == 1) {
            temp_y1 = -root[1];
            Poly<T> temp_right_poly({a*a/4.0 - b + temp_y1, temp_y1*a/2. - c, temp_y1*temp_y1/4. -d});
            std::vector<Poly<T>> temp_res = temp_right_poly.solve();
            if ( temp_res.size() == 2) {
                y1 = temp_y1;
                right_part = temp_res[0];
                break;
            }
        }
    }

    Poly<T> left_poly = Poly<T>(std::vector<T>({1., a/2., y1/2.}));

        std::cout << left_poly << " " << right_part << "\n";
        std::cout << " MINUS " << left_poly - right_part << "\n";
        std::cout << " PLUS " << left_poly + right_part << "\n";

    std::vector<Poly<T>> res1 = (left_poly - right_part).solve();
    std::vector<Poly<T>> res2 = (left_poly + right_part).solve();
    res1.insert(res1.end(), res2.begin(), res2.end());
    return res1;
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


template<typename T>
std::vector<std::pair<int, std::complex<T>>> solve_with_eigenvalues(const Poly<T>& poly){
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
        res.push_back({1, roots[i]});
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
        X[1] = -x; p /= X; ++d;
        res.push_back({m, x});
        X[1] = 0;
    }
    return res;
}

template<typename T>
std::vector<std::pair<int, std::complex<T>>> solve_numerical(const Poly<T>& poly, const std::string& type="eigenvalues") {
    if (type == "laguerre") {
        return solve_laguerre(poly);
    }
    return solve_with_eigenvalues(poly);

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
void Poly<T>::update_real() {
    size_t i = 0;
    while(std::abs(coef_[i]) < 1e-12 && i < deg_) {
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
Poly<T>::Poly(const std::initializer_list<T> &coefs) {
    coef_ = coefs;
    deg_ = coef_.size() - 1;
}

template<typename T>
T Poly<T>::eval(const T &x) const {
    if (deg_==0) {
        return coef_[0];
    }
    T sum = 0;
    T t = 1.0;
    for (auto it = coef_.rbegin(); it != coef_.rend(); ++it){
        sum += t * (*it);
        t *= x;
    }
    return sum;
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
    for (size_t i = 0; i <= deg; ++i) {
        ans[i] = p[i] / value;
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
        throw  std::logic_error("deg lhs poly < deg rhs poly. Operator /");
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
                            const Poly<T> & p) {
    return PrintVec(out, p.GetCoef());
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
