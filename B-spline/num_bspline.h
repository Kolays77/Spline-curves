#include "../include/Point.h"
#include "../include/tools.h"


/*
Алгоритм де Бура для вычисления точки  кривой для k-го сегмента.
Параметры:
    t - параметр t \in [0, 1];
    k - номер сегмента;
    knots - вектор узлов;
    points - вектор констрольных точек на плоскости;
    p - степень кривой.
*/

template< typename T>
Point<T> num_de_boor(double t,
             int k,
             std::vector<T> & knots,
             std::vector<Point<T>>& points,
             int p){

    size_t dim = points[0].dim;
    std::vector<Point<T>> d(p+1);

    for (int i = 0; i < p+1; ++i)
        d[i] = points[i + k - p];

    T alpha, den;

    for (int r = 1; r < p+1; ++r) {
        for (int j=p ; j > r-1; --j) {
            den = knots[j+1+k-r] - knots[j+k-p];
            alpha = (t - knots[j+k-p]) / den;
            for (int i = 0; i < dim; ++i){
                d[j][i] = alpha*d[j][i] + (1-alpha)*d[j-1][i];
            }
        }
    }
    return d[p];
}

/*
Функция создания B-сплайн кривой численно.
Параметры:
    knots - вектор узлов;
    cv - вектор констрольных точек на плоскости;
    p - степень кривой;
    N - требуемое количество точек на выходе.
Результат:
    points - N точек на кривой. 
*/
template<typename T>
std::vector<Point<T>> create_curve(int p,
                                   std::vector<T>& knots,
                                   std::vector<Point<T>>& cv,
                                   int N){

    std::pair<int, int> domain = {p, knots.size() - p - 1};
    std::vector<int> ks = create_intervals(domain, knots);
    std::vector<T> ts = linspace(knots[domain.first], knots[domain.second], N);
    T t = knots[domain.first];
    std::vector<Point<T>> points(N);
    int i = 0;
    for(int k : ks)
        while (t <= knots[k+1] && i < N){
            points[i] = num_de_boor(t, k, knots, cv, p);
            ++i;
            t = ts[i];
        }
    return points;
}
