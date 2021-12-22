#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <iomanip>

typedef std::vector<long double> vec1d;
using namespace std;
 
vec1d linspace(long double start, 
                long double end, 
                long double num)
{
    vec1d linspaced;
    if (num == 0) 
        return linspaced; 
    if (num == 1) {
        linspaced.push_back(start);
        return linspaced;
    }

    long double delta = (end - start) / (num - 1);
    for(int i=0; i < num-1; ++i)
        linspaced.push_back(start + delta * i);
    linspaced.push_back(end); 
  return linspaced;
}


template<typename T>
void print_vec(T& v){
    for (auto vv : v){
        cout << vv << "\t";
    }
    cout << endl;
}

vec1d create_knots(int n, int p){
    vec1d res(n + p + 1);
    for (int i = 0; i <= p; ++i) res[n + p - i] = n-p;
    for (int i = 1; i < n - p; ++i) res[p + i] = i;
    //for (int i = 0; i < n + p + 1; ++i) res[i] /= n-p;
    return res;
}

std::vector<int> create_intervals(std::pair<int, int> domain, vec1d knots){
    vector<int> intervals;
    for (int i = domain.first; i < domain.second; ++i)
        if (knots[i] != knots[i+1])
            intervals.push_back(i);
    
    return intervals;
}

void print_vec1_file(const vec1d& vec, const string& file){
    std::ofstream out;
    out.open(file);
    out << setprecision(25) << fixed;
    for (auto v : vec)
        out << v << endl;
    out.close();
}

void print_coef_file(const vector<vec1d>&  coefs, const string& file){
    std::ofstream out;
    out.open(file);
    out << setprecision(25) << fixed;

    for (auto coef : coefs){
        for (auto c : coef){ 
            out << c << "\n";
        }
        out << "\n\n";
    }
    out.close();
}

vec1d load_vector(const string file){
    vec1d vec;
    std::ifstream in;
    in.open(file);
    if (!in) { 
        throw logic_error("No here file cv_x.in or cv_y.in");
    }
    long double temp = 0.0;
    while( in >> temp)  
        vec.push_back(temp);
    in.close();
    return vec;
}


long double poly_eval_naive(vec1d& coef, long double value){
    long double sum = 0.0;
    long double t = 1.0;
    for (auto it = coef.rbegin(); it != coef.rend(); ++it){
        sum += t * (*it);
        t *= value;
    }
    return sum;
}

long double two_sum (long double &t, long double a, long double b) {
  long double s = a+b;
  long double bs = s-a;
  long double as = s-bs;
  t = (b-bs) + (a-as);
  return s;
}

long double two_sum_res (long double a, long double b) {
  long double t = 0.0;
  long double s = two_sum(t,a,b);
  return s + t; 
  }

long double  sum_rump (const vec1d &X) {
  long double s=0.0, c=0.0, e;
  for (long double x: X) {
    s = two_sum (e, s, x);
    c += e;
  }
  return s+c;
}

long double  two_prod (long double &t, 
                      long double a, 
                      long double b) {
  long double p = a*b;
  t = fma (a, b, -p);  // t = a*b-p
  return p;
}

inline long double  graillat (const vec1d &A, long double x) {
  long double s=0.0, c=0.0, p, pi, t;

  for (auto a: A) {
    p = two_prod (pi, s, x);
    s = two_sum (t, p, a);
    c *= x;
    c += pi+t; 
  }
  return s+c;
}
