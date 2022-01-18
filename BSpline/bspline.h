#include <vector>
#include "Point.h"
#include "Poly.h"
#include "tools.h"

template<typename T>
Point<Poly<T>> de_boor(int k,
                     std::vector<T>& knots,
                     std::vector<Point<T>>& points,
                     int p) {

    size_t dim = points[0].dim;
    std::vector<Point<Poly<T>>> d((size_t)p+1, Point<Poly<T>>(dim));

    for (int i = 0; i < p+1; ++i)
        for (int j = 0; j < dim; ++j){
            d[i][j] = Poly<T>{points[i + k - p][j]};
        }

    T B,D, den;
    for (int r = 1; r < p+1; ++r) {
        for (int j=p ; j > r-1; --j) {
            den = knots[j+1+k-r] - knots[j+k-p];
            B = - knots[j+k-p];
            D = knots[j+1+k-r];
            for (size_t i = 0; i < dim; ++i){
                Poly<T> temp1 = d[j][i] * Poly<T>{1,0} + B * d[j][i] ;
                Poly<T> temp2 = d[j-1][i] * Poly<T>{-1,0} + D * d[j-1][i];
                d[j][i] = (temp1 + temp2) / den;
            }
        }
    }
    return d[p];
}

template<typename T>
std::vector<Point<T>> create_curve(int p,
                                   std::vector<T>& knots,
                                   std::vector<Point<T>>& cv,
                                   int N,
                                   std::vector<Point<Poly<T>>>& coefs_){

    int dim = cv[0].dim;
    std::pair<int, int> domain = {p, knots.size() - p - 1};
    std::vector<int> ks = create_intervals(domain, knots);
    std::vector<T> ts = linspace(knots[domain.first], knots[domain.second], N);
    T t = knots[domain.first];

    std::vector<Point<T>> points(N, Point<T>(dim));
    std::vector<Point<Poly<T>>> coefs;
    int i = 0;
    for(int k : ks){
        coefs.push_back(de_boor( k, knots,cv, p));
        while (t <= knots[k+1] && i < N){
            for (int d = 0; d < dim; ++d){
                //coefs[coefs.size()-1][d].update();
                points[i][d] = coefs[coefs.size()-1][d].At(t);
            }
            ++i;
            t = ts[i];
        }
    }
    coefs_ = coefs;
    return points;
}

template<typename T>
void save_coefs(const std::vector<Point<Poly<T>>> & coefs) {
    if (coefs[0].dim == 1) {
        std::ofstream out("coefs.out");
        for (const Point<Poly<T>>& point : coefs) {
            out << point[0] << "\n";
        }
    } else
        if (coefs[0].dim == 2) {
        std::ofstream out_x("coefs_x.out");
        std::ofstream out_y("coefs_y.out");
        for (const Point<Poly<T>>& point : coefs) {
            out_x << point[0] << "\n";
            out_y << point[1] << "\n";
        }

    } else
        if (coefs[0].dim == 3) {
        std::ofstream out_x("coefs_x.out");
        std::ofstream out_y("coefs_y.out");
        std::ofstream out_z("coefs_z.out");

        for (const Point<Poly<T>>& point : coefs) {
            out_x << point[0] << "\n";
            out_y << point[1] << "\n";
            out_z << point[2] << "\n";

        }
    }
}


template <typename T>
std::vector<Point<T>> test_bspline(int p){
    std::vector<Point<T>> cv = load_points<T>("points.in");
    std::vector<T> knots = create_knots<T>(cv.size(), p);
    std::vector<Point<Poly<T>>> coefs;

    std::vector<Point<T>> points = create_curve<T>(p, knots, cv, 1000, coefs);
    save_points("points.out", points);
    save_coefs(coefs);

    return points;
}