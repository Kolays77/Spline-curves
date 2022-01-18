//
// Created by reima on 30.12.2021.
//

#ifndef BSPLINE_TOOLS_H
#define BSPLINE_TOOLS_H

#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <iomanip>

#include "Point.h"

template<typename T>
std::vector<T> linspace(T start,
               T end,
               int num)
{
    std::vector<T> linspaced;
    if (num == 0)
        return linspaced;
    if (num == 1) {
        linspaced.push_back(start);
        return linspaced;
    }

    T delta = (end - start) / (num - 1);
    for(int i=0; i < num-1; ++i)
        linspaced.push_back(start + delta * (T)i);
    linspaced.push_back(end);
    return linspaced;
}


template<typename T>
std::vector<T> create_knots(int n, int p){
    std::vector<T> res(n + p + 1);
    for (int i = 0; i <= p; ++i) res[n + p - i] = n-p;
    for (int i = 1; i < n - p; ++i) res[p + i] = i;
    for (int i = 0; i < n + p + 1; ++i) res[i] /= n-p;
    return res;
}

template<typename T>
void print(const std::vector<T>& vec) {
    for (const T& v : vec) {
        std::cout << v << "\n";
    }
}

template<typename T>
std::vector<T> ones(size_t n) {
    return std::vector<T>(n, T(1.));
}

template<typename T>
std::vector<int> create_intervals(std::pair<int, int> domain, std::vector<T> knots){
    std::vector<int> intervals;
    for (int i = domain.first; i < domain.second; ++i)
        if (knots[i] != knots[i+1])
            intervals.push_back(i);
    return intervals;
}

template<typename T>
T compare_points(const std::vector<Point<T>>& points1,
                 const std::vector<Point<T>>& points2) {
    T res = 0.0;
    if (points1.size() != points2.size()) {
        throw std::logic_error("points1.size() != points2.size()");
    }
    for (int i = 0; i < points1.size(); ++i) {
        res += dist(points1[i], points2[i]);
    }
    return res;
}


template<typename T>
std::vector<T> load_vector(const std::string& file){
    std::vector<T> points;
    std::ifstream in(file);
    if (!in)
        throw std::logic_error("No here file " + file);

    T temp = 0.0;
    while(in >> temp)
        points.push_back(temp);
    return points;
}


#endif //BSPLINE_TOOLS_H
