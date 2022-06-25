#include "../../nurbs.h"
#include "../../nurbs2.h"
#include "../../../include/split_curve.h"
#include "../../../include/plot.h"

template <typename T>
void test_nurbs_integrals(int p, int n ) {
    std::vector<Point<T>> cv = generate_points_sorted<T>(n);        
    save_points("cv.out", cv);
    
    std::complex<T> int0_1;
    std::complex<T> int0_2;
    
    std::complex<T> int1;
    std::complex<T> int2;

    // FULL UNIFORM    
    NURBS<T>  nurbs1(p, 1.0, 2.0, cv);
    NURBS2<T> nurbs2(p, 1.0, 2.0, cv);
    
    nurbs1.save_coefs();
    nurbs2.save_coefs();
    
    std::vector<Point<T>> points_1 = nurbs1.get_points(1000);
    std::vector<Point<T>> points_2 = nurbs2.get_points(1000);
    
    int0_1 = std::complex<T>{nurbs1.numerical_integral()};
    int0_2 = std::complex<T>{nurbs2.numerical_integral()};
    
    int1 = nurbs1.analytic_integral1(1);
    int2 = nurbs2.analytic_integral1(1);
    
    std::cout << std::setprecision(15);
    std::cout << "Num: "        << int0_1 << "\n";
    std::cout << "Analytic: "    << int1 << "\n";
    std::cout << "Error: " << std::abs(int0_1 - int1) << "\n";
    std::cout << "========================\n";
    
    std::cout << "Num: "        << int0_2 << "\n";
    std::cout << "Analytic: "    << int2 << "\n";
    std::cout << "Error: " << std::abs(int0_2 - int2) << "\n";

    plot_curve(cv, points_1, "plot1.png", "NURBS1. p =" + std::to_string(p));
    plot_curve(cv, points_2, "plot2.png", "NURBS2. p =" + std::to_string(p));
    
    std::cout << "========================\n";
    std::vector<T> knots1;
    std::vector<T> knots2;

    std::vector<T> weights1;
    std::vector<T> weights2;

    std::vector<Point<T>> cv1;
    std::vector<Point<T>> cv2;
    int m = nurbs1.knots.size();

    //T new_knot = (nurbs1.knots[m / 2 + 1] + nurbs1.knots[m / 2]) / 2.0 ;
    T new_knot = 0.5;
    split_curve_new_knot<T>(p, nurbs1.knots, cv, nurbs1.weights, new_knot, knots1, knots2, cv1, cv2, weights1, weights2);          
    
    
    NURBS2<T> nurbs_part1(p, knots1, weights1, cv1);
    NURBS2<T> nurbs_part2(p, knots2, weights2, cv2);


    save_vector(nurbs_part1.cv, "part1/cv.out");
    save_vector(nurbs_part2.cv, "part2/cv.out");

    save_vector(knots1, "part1/knots.out");
    save_vector(knots2, "part2/knots.out");


    nurbs_part1.save_coefs("part1/");
    nurbs_part2.save_coefs("part2/");

    std::vector<Point<T>> points_part1 = nurbs_part1.get_points(1000);
    std::vector<Point<T>> points_part2 = nurbs_part2.get_points(1000);
    
    plot_curve_(cv1, points_part1, "curve 1");
    plot_curve(cv2,  points_part2, "plot_split.png", "NURBS by parts", "curve 2");

    PLOT_END(); 

    std::complex<T> int_0_parts = std::complex(nurbs_part1.numerical_integral() + nurbs_part2.numerical_integral());
    std::complex<T> int_1_parts = nurbs_part1.analytic_integral1(1) + nurbs_part2.analytic_integral1(1);
    
    std::cout << "Num: "        << int_0_parts << "\n";
    std::cout << "Analytic: "    << int_1_parts << "\n";
    std::cout << "Error: " << std::abs(int_0_parts - int_1_parts) << "\n";

}


int main(int argc, char const *argv[]) {
    int p, n;
    p = std::atoi(argv[1]);
    n = std::atoi(argv[2]);
    test_nurbs_integrals<long double>(p, n);    
    return 0;
}
