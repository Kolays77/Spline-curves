#include <vector>
#include <cmath>
#include <complex>

const long double twoPI = 2.0L*M_PI;
const long double onethird = 1./3L; 
const long double eps = 1e-13;

int solveP3(long double *x, long double a,long double b,long double c) {
	long double a2 = a*a;
  	long double a_3 = a*onethird;
  	long double q  = a_3*a_3 - b*onethird;
	long double r  = a_3*(a_3*a_3-b*0.5)+c*0.5;
	long double r2 = r*r;
	long double q3 = q*q*q;
	long double A, B;
    if (r2 < q3) {
        long double t=r/std::sqrt(q3);
        if(t<-1) t=-1;
        if(t> 1) t= 1;
        
		t=std::acos(t);
        q=-2*std::sqrt(q);
        x[0]=q*std::cos(t/3)-a_3;
        x[1]=q*std::cos((t+twoPI)*onethird)-a_3;
        x[2]=q*std::cos((t-twoPI)*onethird)-a_3;
        return 3;
    }
    else {
        A = - std::cbrt(std::abs(r)+std::sqrt(r2-q3));
        A = r < 0 ? -A : A;
		B = (std::abs(A) < eps ? 0 : q/A);
        x[0] =(A+B)-a_3;
        x[1] =-0.5*(A+B)-a_3;
        x[2] = 0.5*sqrt(3.)*(A-B);
        if(std::abs(x[2])<eps) { x[2]=x[1]; return 2; }
    return 1;
    }
    return 0;
}

// solve x^deg + num = 0
std::vector<std::complex<long double>> solve_sqrt(int deg, long double num) {
	std::complex<long double> r = std::pow(std::abs(num), 1.0/deg);
	long double phi = num < 0.0 ? 0.0 : M_PI;
	std::vector<std::complex<long double>> res;
	for (int i = 0 ; i < deg; ++i) {
		res.push_back(r * std::complex<long double>(std::cos((phi + twoPI * i) / deg), std::sin((phi + twoPI * i) / deg)));
	}
	return res;
}

std::vector<std::complex<long double>> solveP4(long double a, long double b, long double c, long double d) {
	if (std::abs(a) < eps && std::abs(b) < eps && std::abs(c) < eps) 
		return solve_sqrt(4, d);

	// y^3 − b*y^2 + (ac−4d)*y − a^2*d−c^2+4*b*d = 0
	long double a3 = -b;
	long double b3 =  a*c -4.*d;
	long double c3 = -a*a*d - c*c + 4.*b*d;

	long double x3[3];
	unsigned int iZeroes = solveP3(x3, a3, b3, c3);

	long double q1, q2, p1, p2, D, sqD, y;

	y = x3[0];
	if(iZeroes != 1) {
		if(std::abs(x3[1]) > std::abs(y)) y = x3[1];
		if(std::abs(x3[2]) > std::abs(y)) y = x3[2];
	}
	D = y*y - 4*d;
	if(std::abs(D) < eps) {
		q1 = q2 = y * 0.5;
		D = a*a - 4*(b-y);
		if(std::abs(D) < eps) 
			p1 = p2 = a * 0.5;
		else {
			sqD = std::sqrt(D);
			p1 = (a + sqD) * 0.5;
			p2 = (a - sqD) * 0.5;
		}
	}
	else {
		sqD = std::sqrt(D);
		q1 = (y + sqD) * 0.5;
		q2 = (y - sqD) * 0.5;
		p1 = (a*q1-c)/(q1-q2);
		p2 = (c-a*q2)/(q1-q2);
	}
    std::vector<std::complex<long double>> retval(4);
	D = p1*p1 - 4*q1;
	if(D < 0.0) {
		retval[0].real( -p1 * 0.5 );
		retval[0].imag( std::sqrt(-D) * 0.5 );
		retval[1] = std::conj(retval[0]);
	}
	else {
		sqD = sqrt(D);
		retval[0].real( (-p1 + sqD) * 0.5 );
		retval[1].real( (-p1 - sqD) * 0.5 );
	}
	D = p2*p2 - 4*q2;
	if(D < 0.0) {
		retval[2].real( -p2 * 0.5 );
		retval[2].imag( std::sqrt(-D) * 0.5 );
		retval[3] = std::conj(retval[2]);
	}
	else {
		sqD = std::sqrt(D);
		retval[2].real( (-p2 + sqD) * 0.5 );
		retval[3].real( (-p2 - sqD) * 0.5 );
	}
    return retval;
}

std::vector<std::complex<long double>> solveP4(std::vector<long double> poly) {
    return solveP4(poly[1]/poly[0], poly[2]/poly[0], 
                    poly[3]/poly[0], poly[4]/poly[0]);
}
