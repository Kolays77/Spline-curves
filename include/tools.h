#pragma once
#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <iomanip>
#include <random>
#include <chrono>
#include <algorithm>

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
std::vector<T> generate_vector(int n, T t0= 0.0, T t1 = 0.0) {
    std::vector<T> res(n);
    std::mt19937_64 rng(std::chrono::steady_clock::now().time_since_epoch().count());
    for (int i = 0; i < n; ++i) {
        res[i] = std::uniform_real_distribution<T>(t0, t1)(rng);
    }
    return res;
}



template<typename T>
std::vector<T> generate_vector_sorted(int n, T t0= 0.0, T t1 = 0.0) {
    std::vector<T> res = generate_vector<T>(n, t0, t1);
    std::sort(res.begin(), res.end());
    return res;
}


template<typename T>
std::vector<Point<T>> generate_points(int n, T x0 = 0.0, T x1 = 1.0, T y0 = 0.0, T y1 = 1.0) {
    std::vector<T> x;
    std::vector<T> y;
    std::vector<Point<T>> res(n, Point<T>{std::vector<T>{0.0, 0.0}});

    x = generate_vector(n, x0, x1);
    y = generate_vector(n, y0, y1);
    
    for (int i = 0; i < n; ++i) {
        res[i][0] = x[i];
        res[i][1] = y[i];
    }
    return res;
}


template<typename T>
std::vector<Point<T>> generate_points_sorted(int n, T x0 = 0.0, T x1 = 1.0, T y0 = 0.0, T y1 = 1.0) {
    // SORTED BY X
    std::vector<Point<T>> res = generate_points(n, x0, x1, y0, y1);
    std::sort(res.begin(), res.end());
    return res;
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
std::vector<T> compare_points_vector(const std::vector<Point<T>>& points1,
                 const std::vector<Point<T>>& points2) {
    
    if (points1.size() != points2.size()) {
        throw std::logic_error("points1.size() != points2.size()");
    }
    int n = points1.size();
    std::vector<T> errors(n);
    for (int i = 0; i < n; ++i) {
        errors[i] = dist(points1[i], points2[i]);
    }
    return errors;
}


template<typename T>
void save_vector(const std::vector<T>& vec,
                           const std::string& file){

    std::ofstream out(file);
    out << std::setprecision(16);
    for (const T& v : vec) {
        out << v << "\n";
    }
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

void generate_n_file_points(int n, int len_cv, std::string dir) { 
    for (int i = 1; i <= n; ++i) {
        std::string path = dir + std::to_string(i) + ".in"; 
        std::ofstream out(path);
        std::vector<Point<long double>> points = generate_points_sorted<long double>(len_cv);
        save_points<long double>(path, points);
    }
}


template<typename T>
T curve_length(const std::vector<Point<T>>& points){
    size_t n = points.size();
    T sum = 0.0;
    for (size_t i = 0; i < n-1; ++i) {
        sum += dist(points[i], points[i+1]);
    }
    return sum;
};
