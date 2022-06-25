#include "../include/Point.h"
#include "../include/Poly.h"
#include "../include/tools.h"
#include "../include/integral.h"



//Модифицированный алгоритм де Бура для построения многочленов P_x, P_y, определяющих сегмент кривой под номером k.
//Параметры:
//    k - номер сегмента;
//    knots - вектор узлов;
//    points - вектор контрольных точек (CV) на плоскости;
//    p - степень кривой.

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



//    Класс B-сплайн кривой, в конструкторе которого строятся  многочлены для каждого сегмента кривой.
//    Параметры:
//       p - степень кривой;
//       knots - вектор узлов;
//       cv - вектор контрольных точек. 

template<typename T>
struct Bspline {
    int p; // Степень кривой.
    int dim; // Размерность пространства. Для плоскости dim=2.
    int N_segments; // Количество сегментов

    std::vector<T> knots;
    std::vector<Point<T>> cv;
    std::pair<int, int> domain; // интервал изменения параметра t для кривой.
    std::vector<int> ks; 
    std::vector<Point<Poly<T>>> coefs;

public:
    Bspline(int p_,
            std::vector<T> &knots_,
            std::vector<Point<T>> &cv_) {
        p = p_;
        knots = knots_;
        cv = cv_;
        dim = cv[0].dim;
        domain = {p, knots.size() - p - 1};
        ks = create_intervals(domain, knots);
        N_segments = ks.size() - 1;
        create_coefs();
    }

    void create_coefs() {
        for (int i = 0; i < N_segments; ++i) {
            coefs.push_back(de_boor(ks[i], knots, cv, p));
        }
        
        for (auto &coef: coefs) {
            for (int d = 0; d < dim; ++d) {
                coef[d].update_real();
            }
        }
    }

    // Метод получения N точек на кривой
    std::vector<Point<T>> get_points(int N) {
        std::vector<Point<T>> points(N, Point<T>(dim));
        std::vector<T> ts = linspace(knots[domain.first], knots[domain.second], N);
        T t = knots[domain.first];
        int i = 0;
        for(int j = 0; j < N_segments; ++j){
            while (t <= knots[ks[j+1]] && i < N){
                for (int d = 0; d < dim; ++d){
                    points[i][d] = coefs[j][d].At(t);
                }
                ++i; t = ts[i];
            }
        }
        return  points;
    }

    // Метод нахождения производной по икс B-сплайн кривой.
    // Результат:
    // N - штук значений dy/dx
    std::vector<T> get_dy_dx_points(int N) {
        std::vector<T> points(N, 0.0);
        std::vector<T> ts = linspace(knots[domain.first],
                                     knots[domain.second], N);
        T t = knots[domain.first];
        std::vector<Point<Poly<T>>> coefs_der;
        for (auto &c: coefs) {
            Point<Poly<T>> temp_point(dim);
            for (int d = 0; d < dim; ++d) {
                temp_point[d] = c[d].der();
            }
            coefs_der.push_back(temp_point);
        }

        int i = 0;
        int j = 0;
        for (int k: ks) {
            while (t <= knots[k + 1] && i < N) {
                points[i] = coefs_der[j][1].At(t) / coefs_der[j][0].At(t);
                ++i;
                t = ts[i];
            }
            ++j;
        }
        return points;
    }

    std::vector<Point<Poly<T>>> get_coefs(){
        return coefs;
    }

    std::vector<int> get_ks() { 
        return ks;
    }
    
    std::vector<T> get_knots() {
        return knots;
    }

    // Сохранение коэффициентов многочленов для каждого сегмента в файл. 
    void save_coefs() {
        if (coefs[0].dim == 1) {
            std::ofstream out("coefs.out");
            for (const Point<Poly<T>> &point: coefs) {
                out << point[0] << "\n";
            }
        } else if (coefs[0].dim == 2) {
            std::ofstream out_x("coefs_x.out");
            std::ofstream out_y("coefs_y.out");
            for (const Point<Poly<T>> &point: coefs) {
                out_x << point[0] << "\n";
                out_y << point[1] << "\n";
            }

        } else if (coefs[0].dim == 3) {
            std::ofstream out_x("coefs_x.out");
            std::ofstream out_y("coefs_y.out");
            std::ofstream out_z("coefs_z.out");

            for (const Point<Poly<T>> &point: coefs) {
                out_x << point[0] << "\n";
                out_y << point[1] << "\n";
                out_z << point[2] << "\n";

            }
        }
    }

    // Метод нахождения интеграла аналитически (непосредственно).
    T analytic_integral() {
        T sum = 0.0;
        for (int i = 0; i < N_segments; ++i) { 
            sum += (coefs[i][1] * coefs[i][0].der()).integral(knots[ks[i]], knots[ks[i+1]]);
        }
        return sum;
    }
    
    // Метод нахождения интеграла численно. Используется numerical_integral из include/integral.h
    T integral() {
        T sum = 0.0;
        for (int i = 0; i < N_segments; ++i) { 
            sum += numerical_integral(coefs[i][1] * coefs[i][0].der(), knots[ks[i]], knots[ks[i+1]]);
        }
        return sum;
    }
};
