#include "../include/Point.h"
#include <cassert>

template <typename T>
int find_span(int p, std::vector<T>& knots, T new_knot) {
    for (int i = 0; i < knots.size() - 1; ++i) {
        if (knots[i] <= new_knot && knots[i+1] > new_knot) {
            return i;
        }
    }
    return 0;
}

// interior knot may inserted only p times !
// TODO заменить пушбек на создание предварительной аллокацией
template <typename T>
int insert_knot(int p, std::vector<T>& knots, std::vector<Point<T>>& cv, std::vector<T>& w, T new_knot) {
    int n_span = find_span<T>(p, knots, new_knot);
    std::vector<Point<T>> Q;
    std::vector<T> R_w;
    for (int i = 0; i < cv.size(); ++i) cv[i] *= w[i];
    for (int i = 0; i <= n_span - p; ++i) {
        Q.push_back(cv[i]);
        R_w.push_back(w[i]);
    }
    for (int i = n_span - p + 1; i <= n_span; ++i) {
        if (std::abs(knots[i + p] - knots[i]) < 1e-18) throw "Division by zero condition!";
        
        T a  = (new_knot - knots[i]) / ( knots[i + p] - knots[i]);
        Q.push_back((T(1.0) - a) * cv[i - 1] + a * cv[i]);
        R_w.push_back((T(1.0) - a) * w[i - 1] + a * w[i]);
    }
    for (int i = n_span + 1; i <= cv.size(); ++i) {
        R_w.push_back(w[i-1]);
        Q.push_back(cv[i-1]);
    }
    cv = Q;
    w = R_w;
    for (int i = 0; i < cv.size(); ++i) cv[i] /= w[i];
    knots.insert(knots.begin() + n_span + 1, new_knot);
    return n_span;
}

template <typename T>
void normalize_knot_vector(std::vector<T>& knots) {
    if (std::abs(knots[0]) > 1e-10) {
        T first = knots[0];
        for (int i = 0; i < knots.size(); ++i) {
            knots[i] -= first;
        }
    } 

    if (std::abs(knots[knots.size() - 1] - T(1.0)) > 1e-10) 
        for (int i = 0; i < knots.size(); ++i) {
            knots[i] /= knots[knots.size() - 1];
        }
}

template <typename T>
void split_curve_new_knot(int p, std::vector<T>& knots, 
                        std::vector<Point<T>>& cv, 
                        std::vector<T>& w, 
                        T new_knot, 
                        std::vector<T>& knots1, 
                        std::vector<T>& knots2, 
                        std::vector<Point<T>>& cv_curve1,
                        std::vector<Point<T>>& cv_curve2,
                        std::vector<T>& weights_curve1,
                        std::vector<T>& weights_curve2) {
                            
    int n =  insert_knot<T>(p, knots, cv, w, new_knot);
    for (int i = 1; i < p; ++i) insert_knot<T>(p, knots, cv, w, new_knot);

    // TODO change push_back to []
    for (int i = 0; i < n + 1; ++i) knots1.push_back(knots[i]);
    for (int i = 0; i < p + 1; ++i) knots1.push_back(knots[n + 1]);
    for (int i = 0; i < p; ++i) knots2.push_back(knots[n+1]);
    for (int i = n+p; i < knots.size(); ++i) knots2.push_back(knots[i]) ;
    
    normalize_knot_vector(knots1);
    normalize_knot_vector(knots2);
    
    for (int i = 0; i < n + 1; ++i) {
        weights_curve1.push_back(w[i]);
        cv_curve1.push_back(cv[i]);
    } 

    for (int i = n; i < cv.size(); ++i) {
        weights_curve2.push_back(w[i]);
        cv_curve2.push_back(cv[i]);
    } 
}
