#pragma once
#include "Poly.h"

long double EPS = 1e-12;

// Для удобства представляется функции С(t) = (P_x/Q(t), P_y/Q(t))
template<typename T>
using Rational = std::pair<Point<Poly<T>>, Poly<T>>;

template<typename T>
T Abs(T value) {
    return value < 0 ? -value : value;
}

// ========== DECOMPOSITION METHODS ========== //

// P/(x-a1)(x-a2) => A1/(x-a1) + A2/(x-a2)
// degP = 0, 1
template<typename T>
std::vector<std::complex<T>> decomp(Poly<T>& P, std::complex<T>a1, std::complex<T> a2) {
    Poly<std::complex<T>> poly = make_complex_poly(P);
    std::complex<T> A1 = poly.eval(a1) / (a1 - a2); // A1
    std::complex<T> A2 = poly.eval(a2) / (a2 - a1); // A2
    return {A1, A2};
}

template<typename T>
std::vector<std::complex<T>> decomp(std::complex<T>a1, std::complex<T> a2) {
    std::complex<T> A1 = std::complex<T>(1.0) / (a1 - a2); // A1
    std::complex<T> A2 = std::complex<T>(1.0) / (a2 - a1); // A2
    return {A1, A2};
}


// P/((x-a1)(x-a2)...(x-a_m)) => A1/(x-a1) + A2/(x-a2) + ... + A_m/(x-a_m)
// degP = 0, 1, ..., a_m-1
// a_i -- DISTINCT ROOTS
template<typename T>
std::vector<std::complex<T>> decomp(Poly<T>& P, std::vector<std::complex<T>>& roots) {
    int n_roots = roots.size();
    std::vector<std::complex<T>> res(n_roots);
    Poly<std::complex<T>> poly = make_complex_poly(P);
    for (int i = 0; i < n_roots; ++i) {
        // FIND A_I
        std::complex<T> temp_den = 1;
        for (int k = 0; k < n_roots; ++k)
            temp_den = i != k ? temp_den * (roots[i] - roots[k]) : temp_den;
        res[i] = poly.eval(roots[i]) / temp_den;
    }
    return res;
}

// 1 / ((x-a1)(x-a2)...(x-a_m)) => A1/(x-a1) + A2/(x-a2) + ... + A_m/(x-a_m)
template<typename T>
std::vector<std::complex<T>> decomp(std::vector<std::complex<T>>& roots) {
    int n_roots = roots.size();
    std::complex<T> one = T(1.0);
    std::vector<std::complex<T>> res(n_roots);
    for (int i = 0; i < n_roots; ++i) {
        // FIND A_I
        std::complex<T> temp_den = 1;
        // Accumulate denominator
        for (int k = 0; k < n_roots; ++k)
            temp_den = i != k ? temp_den * (roots[i] - roots[k]) : temp_den;
        res[i] = one / temp_den;
    }
    return res;
}

// ========== INTEGRAL METHODS ========== //


// integral NUM / (x-a)
// integral NUM / (x-a)^n

template<typename T>
std::complex<T> integral_type_1(std::complex<T> NUM,
                                std::complex<T> a, int n,  T from,  T to) {
    //if (std::abs(den.eval(to)) < EPS || std::abs( den.eval(from)) < EPS)
    //    throw std::logic_error("Integral does not converge");
    std::complex<T> A(from);
    std::complex<T> B(to);

    if(n == 1)
        return NUM * ( std::log(B - a) - std::log(A - a)) ;
    else {
        return  NUM / std::complex<T>{n-1.} * (std::complex<T>{1.0}/std::pow(A-a, T(n)-1.) -
                                               std::complex<T>{1.0}/std::pow(B-a, T(n)-1.));
    }
}


template<typename T>
std::complex<T> integral_triple_roots(std::complex<T> compl_num,
                                      std::complex<T> a1,
                                      std::complex<T> a2,
                                      std::complex<T> a3,
                                      T from,
                                      T to) {

    std::complex<T> sum = 0.0;
    std::complex<T> A1 = compl_num / (a1 - a2) / ( a1 - a3);
    std::complex<T> A2 = compl_num / (a2 - a1) / ( a2 - a3);
    std::complex<T> A3 = compl_num / (a3 - a1) / ( a3 - a2);
    sum += integral_type_1(A1, a1, 1, from, to);
    sum += integral_type_1(A2, a2, 1, from, to);
    sum += integral_type_1(A3, a3, 1, from, to);

    return sum;
}


// integral A / ((x+a)(x+b)) = - 1 / (a-b) ln(| (x+a) / (x+b)|)
template<typename T>
std::complex<T> integral_type_7(std::complex<T> compl_num,
                                std::complex<T> a,
                                std::complex<T> b,  T from,  T to) {
    std::complex<T> A(from);
    std::complex<T> B(to);
    std::vector<std::complex<T>> res = decomp(a , b);
    return  compl_num * (integral_type_1(res[0], a, 1, from, to) + integral_type_1(res[1], b, 1, from, to));
}

// integral A / ((x-a)(x-b)^2)
template<typename T>
std::complex<T> integral_type_12(std::complex<T> compl_num,
                                 std::complex<T> a,
                                 std::complex<T> b,  T from,  T to) {
    std::complex<T> A(from);
    std::complex<T> B(to);
    std::vector<std::complex<T>> vec_coef = decomp(a, b);

    std::complex<T> int1 = vec_coef[0] * (integral_type_1(vec_coef[0], a, 1, from, to) + 
                                integral_type_1(vec_coef[1], b, 1, from, to));    
    std::complex<T> int2 = integral_type_1(vec_coef[1], b, 2, from, to);
    return compl_num * (int1 + int2);
}


template<typename T>
std::complex<T> integral(Poly<T>& poly,
                         std::vector<std::complex<T>>& roots, T from, T to) {

    int n_roots = roots.size();
    std::complex<T> res = 0.0;
    std::vector<std::complex<T>> vec_coef = decomp(poly, roots);
    for (int i = 0; i < n_roots; ++i) {
        res += integral_type_1(vec_coef[i], roots[i], 1, from, to);
    }
    return res;
}



// integral of (Poly1 / Q) * (Poly2 / Q)
template<typename T>
std::complex<T> integral(Poly<T>& poly1, Poly<T>& poly2,
                         std::vector<std::complex<T>>& vec_root, T from, T to) {

    int n_roots = vec_root.size();
    std::complex<T> sum = 0.0;

    std::vector<std::complex<T>> vec_coef1 = decomp(poly1, vec_root);
    std::vector<std::complex<T>> vec_coef2 = decomp(poly2, vec_root);
    
    for (int i = 0; i < n_roots; ++i) {
        for (int j = i; j < n_roots; ++j) {
            std::complex<T> val = 0.0;
            if (i == j) {
                sum += integral_type_1(vec_coef1[i] * vec_coef2[j], vec_root[i], 2, from, to);;
            } else {
                sum += integral_type_7(vec_coef1[i]*vec_coef2[j] + vec_coef1[j]*vec_coef2[i],
                                       vec_root[i],
                                       vec_root[j], from, to);
            }
        }
    }
    return sum;
}


// integral P/((x-a1)(x-a2) ... (x-a_m)) / (x-a_k)
template<typename T>
std::complex<T> integral(Poly<T>& poly,
                         std::vector<std::complex<T>>& roots, int k, T from, T to) {


    int n_roots = roots.size();
    std::complex<T> res = 0.0;
    std::vector<std::complex<T>> vec_coef = decomp(poly, roots);

    for (int i = 0; i < n_roots; ++i) {
        if (i == k) {
            res += integral_type_1(vec_coef[i], roots[i], 2, from, to);
        } else {
            res += integral_type_7(vec_coef[i], roots[i], roots[k], from, to);
        }
    }
    return res;
}

// integral P1*P2 / Q^2 / (x-a_k)
// where Q = (x-a1)(x-a2) ... (x-a_m)) / (x-a_k)
template<typename T>
std::complex<T> integral(Poly<T>& poly1, Poly<T>& poly2,
                         std::vector<std::complex<T>>& vec_root, int k, T from, T to) {


    int n_roots = vec_root.size();
    std::complex<T> sum = 0.0;

    std::vector<std::complex<T>> vec_coef1 = decomp(poly1, vec_root);
    std::vector<std::complex<T>> vec_coef2 = decomp(poly2, vec_root);

    for (int i = 0; i < n_roots; ++i) {
        for (int j = i; j < n_roots; ++j) {

            if (i == j && j == k) {
                sum += integral_type_1(vec_coef1[i] * vec_coef2[j], vec_root[i], 3, from, to);
            } else if (i == j) {
                sum += integral_type_12(vec_coef1[i]*vec_coef2[j],
                                        vec_root[k],
                                        vec_root[j], from, to);
            } else if (i == k) {
                sum += integral_type_12(vec_coef1[i]*vec_coef2[j]+vec_coef1[j]*vec_coef2[i], vec_root[j],
                                        vec_root[i], from, to);

            } else if (j == k) {
                sum += integral_type_12(vec_coef1[i]*vec_coef2[j]+vec_coef1[j]*vec_coef2[i], vec_root[i],
                                        vec_root[j], from, to);
            } else {
                sum += integral_triple_roots(vec_coef1[i]*vec_coef2[j]+vec_coef1[j]*vec_coef2[i],
                                             vec_root[i], vec_root[j], vec_root[k], from, to);
            }
        }
    }
    return sum;
}


// ========== NUMERICAL INTEGRAL METHODS ========== //
const std::vector<long double> POINTS({-0.9997137267734413,
                                       -0.9984919506395958,
                                        -0.9962951347331251,
                                       -0.9931249370374434,
                                       -0.9889843952429918,
                                       -0.983877540706057,
                                       -0.9778093584869183,
                                       -0.9707857757637064,
                                       -0.9628136542558156,
                                       -0.9539007829254917,
                                       -0.944055870136256,
                                       -0.9332885350430795,
                                       -0.921609298145334,
                                       -0.9090295709825297,
                                       -0.895561644970727,
                                       -0.8812186793850184,
                                       -0.8660146884971647,
                                       -0.8499645278795913,
                                       -0.8330838798884008,
                                       -0.8153892383391762,
                                       -0.7968978923903145,
                                       -0.7776279096494955,
                                       -0.7575981185197072,
                                       -0.7368280898020207,
                                       -0.7153381175730564,
                                       -0.693149199355802,
                                       -0.670283015603141,
                                       -0.6467619085141293,
                                       -0.6226088602037078,
                                       -0.5978474702471788,
                                       -0.5725019326213812,
                                       -0.5465970120650941,
                                       -0.520158019881763,
                                       -0.49321078920819095,
                                       -0.465781649773358,
                                       -0.4378974021720315,
                                       -0.40958529167830154,
                                       -0.38087298162462996,
                                       -0.3517885263724217,
                                       -0.32236034390052914,
                                       -0.292617188038472,
                                       -0.2625881203715035,
                                       -0.23230248184497396,
                                       -0.20178986409573602,
                                       -0.17108008053860327,
                                       -0.14020313723611397,
                                       -0.10918920358006111,
                                       -0.07806858281343663,
                                       -0.046871682421591634,
                                       -0.015628984421543084,
                                       0.015628984421543084,
                                       0.046871682421591634,
                                       0.07806858281343663,
                                       0.10918920358006111,
                                       0.14020313723611397,
                                       0.17108008053860327,
                                       0.20178986409573602,
                                       0.23230248184497396,
                                       0.2625881203715035,
                                       0.292617188038472,
                                       0.32236034390052914,
                                       0.3517885263724217,
                                       0.38087298162462996,
                                       0.40958529167830154,
                                       0.4378974021720315,
                                       0.465781649773358,
                                       0.49321078920819095,
                                       0.520158019881763,
                                       0.5465970120650941,
                                       0.5725019326213812,
                                       0.5978474702471788,
                                       0.6226088602037078,
                                       0.6467619085141293,
                                       0.670283015603141,
                                       0.693149199355802,
                                       0.7153381175730564,
                                       0.7368280898020207,
                                       0.7575981185197072,
                                       0.7776279096494955,
                                       0.7968978923903145,
                                       0.8153892383391762,
                                       0.8330838798884008,
                                       0.8499645278795913,
                                       0.8660146884971647,
                                       0.8812186793850184,
                                       0.895561644970727,
                                       0.9090295709825297,
                                       0.921609298145334,
                                       0.9332885350430795,
                                       0.944055870136256,
                                       0.9539007829254917,
                                       0.9628136542558156,
                                       0.9707857757637064,
                                       0.9778093584869183,
                                       0.983877540706057,
                                       0.9889843952429918,
                                       0.9931249370374434,
                                       0.9962951347331251,
                                       0.9984919506395958,
                                       0.9997137267734413
                                      });


const std::vector<long double> WEIGHTS({0.0007346344905008809,
                                            0.001709392653517807,
                                            0.002683925371554019,
                                            0.003655961201327216,
                                            0.004624450063421818,
                                            0.005588428003865117,
                                            0.00654694845084515,
                                            0.007499073255464816,
                                            0.008443871469668721,
                                            0.009380419653694542,
                                            0.01030780257486916,
                                            0.01122511402318622,
                                               0.012131457662979251,
                                               0.013025947892971715,
                                               0.01390771070371885,
                                               0.014775884527441474,
                                               0.015629621077546098,
                                               0.01646808617614516,
                                               0.017290460568323632,
                                               0.018095940722128407,
                                               0.018883739613374886,
                                               0.01965308749443545,
                                               0.020403232646209593,
                                               0.021133442112527594,
                                               0.02184300241624754,
                                               0.02253122025633626,
                                               0.02319742318525442,
                                               0.023840960265968263,
                                               0.024461202707957153,
                                               0.025057544481579718,
                                               0.025629402910208283,
                                               0.02617621923954582,
                                               0.02669745918357113,
                                               0.02719261344657694,
                                               0.027661198220792507,
                                               0.028102755659101357,
                                               0.028516854322395237,
                                               0.028903089601125278,
                                               0.029261084110638446,
                                               0.029590488059912694,
                                               0.02989097959333295,
                                               0.03016226510516929,
                                               0.030404079526454932,
                                               0.030616186583980524,
                                               0.03079837903115269,
                                               0.030950478850491105,
                                               0.03107233742756666,
                                               0.031163835696210035,
                                               0.03122488425484948,
                                               0.03125542345386349,
                                               0.03125542345386349,
                                               0.03122488425484948,
                                               0.031163835696210035,
                                               0.03107233742756666,
                                               0.030950478850491105,
                                               0.03079837903115269,
                                               0.030616186583980524,
                                               0.030404079526454932,
                                               0.03016226510516929,
                                               0.02989097959333295,
                                               0.029590488059912694,
                                               0.029261084110638446,
                                               0.028903089601125278,
                                               0.028516854322395237,
                                               0.028102755659101357,
                                               0.027661198220792507,
                                               0.02719261344657694,
                                               0.02669745918357113,
                                               0.02617621923954582,
                                               0.025629402910208283,
                                               0.025057544481579718,
                                               0.024461202707957153,
                                               0.023840960265968263,
                                               0.02319742318525442,
                                               0.02253122025633626,
                                               0.02184300241624754,
                                               0.021133442112527594,
                                               0.020403232646209593,
                                               0.01965308749443545,
                                               0.018883739613374886,
                                               0.018095940722128407,
                                               0.017290460568323632,
                                               0.01646808617614516,
                                               0.015629621077546098,
                                               0.014775884527441474,
                                               0.01390771070371885,
                                               0.013025947892971715,
                                               0.012131457662979251,
                                               0.01122511402318622,
                                               0.01030780257486916,
                                               0.009380419653694542,
                                               0.008443871469668721,
                                               0.007499073255464816,
                                               0.00654694845084515,
                                               0.005588428003865117,
                                               0.004624450063421818,
                                               0.003655961201327216,
                                               0.002683925371554019,
                                               0.001709392653517807,
                                               0.0007346344905008809
                                       });

template<typename T>
T numerical_integral(const Poly<T>& poly, T A, T B) {
    auto points = POINTS;
    auto weights = WEIGHTS;
    int n = points.size();
    T sum = 0.0;
    for (int i = 0; i < n; ++i) {
        sum += weights[i]*poly.At((B-A)*points[i]/2 + (A+B)/2);
    }
    return (B-A)/2.0 * sum;
}

// Length of the parametric curve
template<typename T>
T numerical_integral_length(const Poly<T>& x, const Poly<T>& y, T A, T B) {
    auto points = POINTS;
    auto weights = WEIGHTS;
    int n = points.size();
    T sum = 0.0;

    for (int i = 0; i < n; ++i) {
        T t = (B-A)*points[i]/2 + (A+B)/2;
        sum += weights[i] * std::sqrt( std::pow(x.der().At(t),2) +
                                       std::pow(y.der().At(t),2) );
    }
    return (B-A)/2.0 * sum;
}


// integral of a rational function f(t) = num(t)/den(t)
template<typename T>
T numerical_integral(const Poly<T>& num,
                     const Poly<T>& den, T A, T B) {
    auto points  = POINTS;
    auto weights = WEIGHTS;
    int n = points.size();
    T sum = 0.0;
    T x = 0.0;
    for (int i = 0; i < n; ++i) {
        x = T(B-A)*points[i]/2 + T(A+B)/2;
        sum += weights[i]* num.At(x) / den.At(x);
    }
    return (B-A)/2.0 * sum;
}
