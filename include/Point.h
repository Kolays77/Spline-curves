#pragma once
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>

// TODO add precision to output

template<typename T>
struct Point{
    std::vector<T> arr;
    size_t dim{};

    Point();
    explicit Point(std::vector<T> vec);
    Point(std::initializer_list<T> list); // delete ?

    explicit Point(int dim_);
    T operator[](size_t index) const;
    T& operator[](size_t index);
    Point& operator = (const Point & rhs);
    Point& operator*=(const T& rhs);
    Point& operator/=(const T& rhs);
    Point operator/(T rhs);
    bool operator <(const Point& rhs);
};

template<typename T>
Point<T>::Point(std::initializer_list<T> list) { // delete ?
    dim = list.size();
    arr = std::vector<T>(list);
}

template<typename T>
Point<T>::Point(std::vector<T> vec) {
    dim = vec.size();
    arr = vec;
}


template<typename T>
T Point<T>::operator[](size_t index) const {
    // if (index > dim){
    //     throw std::logic_error(" index > dim \n");
    // }
    return arr[index];
}

template<typename T>
T &Point<T>::operator[](size_t index) {
    // if (index > dim){
    //     throw std::logic_error(" index > dim \n");
    // }
    return arr[index];
}

template<typename T>
Point<T>::Point() {
    arr = std::vector<T>();
    dim = 0;
}

template<typename T>
Point<T>::Point(int dim_) {
    dim = dim_;
    arr = std::vector<T>(dim);
}

template<typename T>
Point<T> &Point<T>::operator*=(const T &rhs) {
    for (T& p : arr){
        p *= rhs;
    }
    return *this;
}

template<typename T>
Point<T> &Point<T>::operator/=(const T &rhs) {
    for (T& p : arr){
        p /= rhs;
    }
    return *this;
}

template<typename T>
Point<T> & Point<T>::operator=(const Point &rhs) {
    dim = rhs.dim;
    arr = rhs.arr;
    return *this;
}

template<typename T>
Point<T>  operator + (const Point<T>& point1, const Point<T>& point2) {
    Point<T> temp = point1;
    for (int i = 0; i < temp.dim; ++i) {
        temp[i] += point2[i];
    }
    return temp;
}

template<typename T>
Point<T> operator * (const Point<T>& point, T value) {
    Point<T> temp = point;
    for (T& a : temp.arr){
        a *= value;
    }
    return temp;
}

template<typename T>
Point<T> operator * (T value, const Point<T>& point) {
    return point * value;
}

template<typename T>
Point<T> Point<T>::operator/(T rhs) {
    Point<T> point = *this;
    for (T& a : point.arr) {
        a /= rhs;
    }
    return point;
}


template < typename T >
std::ostream & operator << (std::ostream & out,
                            const Point<T> & p) {
    for (size_t d = 0; d < p.dim; ++d) {
        out << p.arr[d] << "\t";
    }
    return out;
}

template<typename T>
bool Point<T>::operator < (const Point<T>& rhs) {
    return this->arr[0] < rhs[0];
}


template<typename T>
void save_points(const std::string& file, const std::vector<Point<T>>& points){

    std::ofstream out(file);
    out << std::setprecision(16);
    for(const auto& point: points){
        out << point << "\n";
    }
    out.close();
}

template<typename T>
std::vector<Point<T>> load_points(const std::string& file){
    std::vector<Point<T>> points;
    std::ifstream in(file);
    std::string line;
    if (!in) {
        throw std::logic_error("No here file " + file);
    }
    T temp = 0.0;
    while(getline(in, line)){
        std::vector<T> point;
        std::istringstream ss(line);
        while(ss >> temp)
            point.push_back(temp);
        points.push_back(Point<T>(point));
    }
    return points;
}

template<typename T>
T dist(const Point<T>& p1, const Point<T>& p2){
    return sqrt(pow(p1[0] - p2[0], 2.0) + pow(p1[1] - p2[1], 2.0));
}