#include <vector>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <numeric>
#include <fstream>
#include <chrono>
#include <string>
#include <cmath>
#include <fstream>
using namespace std;

#include "tools.hpp"

template<typename T>
void print_vec_file(ofstream& out, T& v){
    for (auto vv : v){
        out << vv << "\t";
    }
    out << endl;
}

vec1d deBoor(int k, vec1d& knots, vec1d& points, int p){
    vector<vec1d> d(p+1);   
    for (int i = 0; i < p+1; ++i){
        d[i] = vec1d({points[i + k - p]});
    }

	long double A,B,C,D, den;
    for (int r = 1; r < p+1; ++r) {
        for (int j=p ; j > r-1; --j) {
            den = knots[j+1+k-r] - knots[j+k-p];
            A = 1.0;
            B = - knots[j+k-p];
        	C = -A;
            D = knots[j+1+k-r] ;
            
            size_t size = d[j].size() + 1;
            vec1d temp1(size, 0.0);
            vec1d temp2(size, 0.0);
            
            for (int s = 0; s < d[j].size(); ++s){
                temp1[s] = A*d[j][s];
                temp2[s] = C*d[j-1][s];  
            }

            for (int s = 1; s < d[j].size() + 1; ++s){
                temp1[s] = two_sum_res(temp1[s], B*d[j][s-1] );
                temp2[s] = two_sum_res(temp2[s], D*d[j-1][s-1] );
            }
            
            for (int s = 0; s < temp1.size(); ++s){
                temp1[s] = two_sum_res(temp1[s], temp2[s]) / den;
            }
            d[j] = temp1;
        }
    }
    return d[p];
}

int main(int argc, char *argv[]){
    if (argc==1){
        cerr << "Add params :  p, n" << endl;
        return 1;
    }

    int p = atoi(argv[1]);
    int N = 10000;
    
    vec1d xs = load_vector("cv_x.in");
    vec1d ys = load_vector("cv_y.in");

    int n = xs.size();
    vec1d knots = create_knots(n, p);
    pair<int, int> domain = {p, knots.size() - p - 1};

    vector<int> intervals = create_intervals(domain, knots);

    vec1d coef_x;
    vec1d coef_y;

    vec1d points_x(N-1);
    vec1d points_y(N-1);

    long double t0 = knots[domain.first];
    long double tn = knots[domain.second];

    int i = 0;
    long double t = t0;    
    vector<vec1d> coefs_x;
    vector<vec1d> coefs_y;

    vec1d ts = linspace(t0, tn, N);

    auto begin = std::chrono::steady_clock::now();
    for(int k : intervals){
        coef_x = deBoor(k, knots, xs, p);
        coef_y = deBoor(k, knots, ys, p);
        coefs_x.push_back(coef_x);
        coefs_y.push_back(coef_y);
        while (t < knots[k+1]){
            points_x[i] = graillat(coef_x, t);
            points_y[i] = graillat(coef_y, t);
            ++i;
            t = ts[i];
        }
    }
    auto end = std::chrono::steady_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "The time Bspline-prec: " << elapsed_ms.count() << " ms\n";
    
	print_vec1_file(points_x, "xs_prec.out");
    print_vec1_file(points_y, "ys_prec.out");
    print_coef_file(coefs_x, "coefs_x_prec.out");
    print_coef_file(coefs_y, "coefs_y_prec.out");
	return 0;
}