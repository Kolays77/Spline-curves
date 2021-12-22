#include <vector>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <numeric>
#include <fstream>
#include <chrono>
#include <string>
#include "tools.hpp"

using namespace std;

long double deBoor(double t, 
			int k, 
			vec1d& knots,
			vec1d& points, 
			int p){
 
    vec1d d(p+1);
    for (int i = 0; i < p+1; ++i)
        d[i] = points[i + k - p];
	
	long double alpha;
	long double den;

    for (int r = 1; r < p+1; ++r) {
        for (int j=p ; j > r-1; --j) {
			den = knots[j+1+k-r] - knots[j+k-p];
			alpha = (t - knots[j+k-p]) / den;
			d[j] = alpha*d[j] + (1-alpha)*d[j-1];
        }
    }
    return d[p];
}


int main(int argc, char *argv[]){
    int p = atoi(argv[1]);
	vec1d xs = load_vector("cv_x.in");
    vec1d ys = load_vector("cv_y.in");
    int n = xs.size();
	vec1d knots = create_knots(n, p);

    pair<int, int> domain = {p, knots.size() - p - 1};

    vector<int> intervals = create_intervals(domain, knots);
	
    long double t0 = knots[domain.first];
    long double tn = knots[domain.second];
	int N = 10000;
    
    vec1d ts = linspace(t0, tn, N);
	long double t = t0;

	vec1d points_x(N-1);
    vec1d points_y(N-1);
    int i = 0;
    auto begin = std::chrono::steady_clock::now();
    
    for(int k : intervals){
        while (t < knots[k+1]){
            points_x[i] = deBoor(t, k, knots, xs, p);
			points_y[i] = deBoor(t, k, knots, ys, p);
			++i;
            t = ts[i];
        }
    }
    auto end = std::chrono::steady_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "The time Numeric Bspline: " << elapsed_ms.count() << " ms\n";
    
	print_vec1_file(points_x, "num_xs.out");
    print_vec1_file(points_y, "num_ys.out");
	return 0;
}