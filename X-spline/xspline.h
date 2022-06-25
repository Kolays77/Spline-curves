#pragma once

#include <stdexcept>
#include "../include/Poly.h"
#include "../include/tools.h"
#include "../include/plot.h"
#include "../include/integral.h"


#define EPS 1e-10


// Для всех методов : T -- RealType (float/double/long double)
// Вспомогательные аналитические (возвращающие многочлен) функции F | f | g | h.

template<typename T>
inline Poly<T> F(Poly<T> u, T p) {
    return u*u*u * (10.0 - p + (2.0*p - 15.0) * u + (6.0-p)*u*u);
}

template<typename T>
inline Poly<T> f(Poly<T> n, T d) {

    return F(n / d, 2.0*d*d);
}

template<typename T>
inline Poly<T> g(Poly<T> u, T q) {
    return (u*(q + u*(2.0*q+u*(8.0-12.0*q+u*(14.0 * q - 11.0 + u * (4.0 - 5.0 * q))))));
}

template<typename T>
inline Poly<T> h(Poly<T> u, T q) {
    return (u * (q + u * (2.0 * q + u*u * (-2.0 * q - u * q))));
}

// Вспомогательные численные функции  F | f | g | h.

template<typename T>
inline T F(T u, T p) {
    return u*u*u * (10.0 - p + (2.0*p - 15.0) * u + (6.0-p)*u*u);
}

template<typename T>
inline T f(T n, T d) {
    return F(n / d, 2.0*d*d);
}

template<typename T>
inline T g(T u, T q) {
    return (u*(q + u*(2.0*q+u*(8.0-12.0*q+u*(14.0 * q - 11.0 + u * (4.0 - 5.0 * q))))));
}

template<typename T>
inline T h(T u, T q) {
    return (u * (q + u * (2.0 * q + u*u * (-2.0 * q + ((-1) *u) * q))));
}

// Вспомогательные численные функции  A0, A1, A2, A3

template<typename T>
T A0(T t, T s_k_1) {
    if(s_k_1 < 0.0)  
        return h(-t, -s_k_1);
    
    if ( s_k_1 >= 0.0 && t < s_k_1)
        return f(t - s_k_1, -1.0 - s_k_1);
    
    return 0.0;
}

template<typename T>
T A1(T t, T s_k_2) {
    if (s_k_2 < 0.0) return g(1.0 - t, -s_k_2);
    return f(t - 1.0 - s_k_2, -1.0 - s_k_2);
}

template<typename T>
T A2(T t, T s_k_1) {
    if(s_k_1 < 0.0) return g(t, -s_k_1);
    return f(t+s_k_1, 1.0 + s_k_1);
}

template<typename T>
T A3(T t, T s_k_2) { 
    if(s_k_2 < 0.0)
        return h(t - 1.0, -s_k_2);
    if(s_k_2 >= 0.0 && t > 1.0 - s_k_2) 
        return f(t-1.0 + s_k_2, 1.0+s_k_2);
    return 0.0;
}

/*
    Функция получения N  точек на сегменте X-сплайн кривой численно.
    Параметры:
        cv_x_0, ... , cv_x_3 - x-составляющие четырём контрольным точкам P1, ..., P3
        cv_y_0, ... , cv_y_3 - y-составляющие четырём контрольным точкам P1, ..., P3
        N - требуемое количество точек на кривой C(t)
*/
template<typename T>
std::pair<std::vector<T>, std::vector<T>> get_points_from_segment(
    T cv_x_0,  T cv_x_1, T cv_x_2, T cv_x_3,
    T cv_y_0,  T cv_y_1, T cv_y_2, T cv_y_3, 
    T s_k_1, T s_k_2, int N = 1000) {    

    std::vector<T> ts = linspace<T>(0.0, 1.0, N);
    std::vector<T> x(N);
    std::vector<T> y(N);

    for (int i = 0; i < N; ++i) {
        T A0_temp = A0(ts[i], s_k_1);
        T A1_temp = A1(ts[i], s_k_2);
        T A2_temp = A2(ts[i], s_k_1);
        T A3_temp = A3(ts[i], s_k_2);
        T den = A0_temp + A1_temp + A2_temp + A3_temp;
        T temp_x = A0_temp*cv_x_0 + A1_temp*cv_x_1 + 
                   A2_temp*cv_x_2 + A3_temp*cv_x_3;
        T temp_y = A0_temp*cv_y_0 + A1_temp*cv_y_1 + A2_temp*cv_y_2 + A3_temp*cv_y_3;
        x[i] = temp_x / den; 
        y[i] = temp_y / den; 
    }
    return {x, y};
}

/*
    Функция получения N  точек на сегменте X-сплайн кривой численно.
    Параметры:
        cv_x - x-составляющие вектора контрольных точек P_i
        cv_y - y-составляющие вектора контрольных точек P_i
        N_per_segment - требуемое количество точек кривой C(t) для каждого сегмента.
*/
template<typename T>
std::pair<std::vector<T>, std::vector<T>> make_curve(std::vector<T> cv_x, 
                                std::vector<T> cv_y, 
                                std::vector<T> s, int N_per_segment = 1000) {
    
    int n = s.size();
    s[0] = 0.0;
    s[n-1] = 0.0;
    s.push_back(s[n-1]);
    s.insert(s.begin(), s[0]);
    cv_x.push_back(cv_x[n-1]);
    cv_x.insert(cv_x.begin(), cv_x[0]);
    cv_y.push_back(cv_y[n-1]);
    cv_y.insert(cv_y.begin(), cv_y[0]);
    int N = cv_x.size();
    std::vector<T> x;
    std::vector<T> y;
    for (int i = 0; i < N-3; ++i) {
        auto res = get_points_from_segment(cv_x[i], cv_x[i+1], cv_x[i+2], cv_x[i+3],
                                           cv_y[i], cv_y[i+1], cv_y[i+2], cv_y[i+3],
                                           s[i+1], s[i+2], N_per_segment);
        x.insert(x.end(), res.first.begin(), res.first.end());
        y.insert(y.end(), res.second.begin(), res.second.end());
    }    
    return {x, y};
}


/*
    Вспомогательный шаблонный класс SubSegment представляет рациональную функцию
    c(t) = (x(t)/den(t), y(t)/den(t)) при t в пределе [from, to]
    Параметры конструктора:
        den - многочлен Q(t) - знаменатель
        x - многочлен x(t)
        y - многочлен y(t)
        from, to - область определения
*/
template<typename T>
struct SubSegment{
    Poly<T> den;
    Poly<T> x;
    Poly<T> y;
    
    T from;
    T to;

    SubSegment(const Poly<T>& den_, const Poly<T>& x_, const Poly<T>& y_, T from_=0.0L, T to_=1.0L) {
        den = den_;
        x = x_;
        y = y_;
        from = from_;
        to = to_;
    }

    void print() {
        std::cout << "den :" << den << "\n";
        std::cout << "x :" << x << "\n";
        std::cout << "y :" << y << "\n";
        std::cout << "from :" << from << "\n";
        std::cout << "to :" << to << "\n";
    }
};


template<typename T> 
struct Xspline {
	// Структура : [segm1 : [[SubSegment<T>[0, a]], [SubSegment<T>[a, 1.0]]], ..., segm_n : ... ]
    std::vector<std::vector<SubSegment<T>>> vv_segments;
	std::vector<T> cv_x;
 	std::vector<T> cv_y;
    std::vector<T> s;
    
    // Число сегментов : N (number of cv) - 1
    int N_segments;

    Xspline(const std::vector<T>& cv_x_, 
            const std::vector<T>& cv_y_, 
            const std::vector<T>& s_) {
		
        cv_x = cv_x_;
		cv_y = cv_y_;
		s = s_;   
        N_segments = cv_x.size() - 1;
		prepare_data();
		create_curve(); 
        polishing_segments();	
    }

    // Подготовка данных
    // Фиктивное добавление первой и последней контрольной точки : для определения первого и последнего сегмента
	void prepare_data() {
		int n = s.size();
	    s[0] = 0.0;
	    s[n-1] = 0.0;
	    s.push_back(s[n-1]);
	    s.insert(s.begin(), s[0]);
	    cv_x.push_back(cv_x[n-1]);
	    cv_x.insert(cv_x.begin(), cv_x[0]);
	    cv_y.push_back(cv_y[n-1]);
	    cv_y.insert(cv_y.begin(), cv_y[0]);
	}

    void polishing_segments(T eps_ = 1e-11) {
        for (auto& segm : vv_segments) {
            for (auto& subsegm : segm) {
                subsegm.den.update_real(eps_);
                subsegm.x.update_real(eps_);
                subsegm.y.update_real(eps_);

                T temp = subsegm.den.normalize();
                subsegm.x /= temp;
                subsegm.y /= temp;
            }
        }

    } 
    
    // Метод создания всех рациональных функций сегментов кривой
	void create_curve() {
		int N = cv_x.size();
		for(int i = 0; i < N - 3; ++i) {
			auto res = make_segment(cv_x[i], cv_x[i+1], cv_x[i+2], cv_x[i+3],
			                        cv_y[i], cv_y[i+1], cv_y[i+2], cv_y[i+3],
			                        s[i+1], s[i+2]);
			vv_segments.push_back(res);	
		}
	}

    // Метод создания одной рациональной функции сегмента X-spline кривой
	std::vector<SubSegment<T>> make_segment(
	    T cv_x_0,  T cv_x_1, T cv_x_2, T cv_x_3,
	    T cv_y_0,  T cv_y_1, T cv_y_2, T cv_y_3, 
	    T s_k_1, T s_k_2) {

		bool fl_A0 = false;
		bool fl_A3 = false;

        bool zero_A0 = false;
		bool zero_A3 = false;


        Poly<T> zero_poly({0.0});

		Poly<T> A0, A1, A2, A3;
		std::vector<SubSegment<T>> res;

		if (s_k_1 < 0.0) { 
			A0 = h(Poly<T>{-1.0, 0.0}, -s_k_1);
			A2 = g(Poly<T>{1.0, 0.0},  -s_k_1);
		} else if(std::abs(s_k_1) < EPS) {
			zero_A0 = true;
			A0 = zero_poly;
			A2 = f(Poly<T>{1.0, s_k_1}, 1.0 + s_k_1);
		} else {
			fl_A0 = true;
			A2 = f(Poly<T>{1.0,   s_k_1},  1.0 + s_k_1);  // [0, 1]
        	A0 = f(Poly<T>{1.0, - s_k_1}, -1.0 - s_k_1); // [0, s_k_1]
		}

		if (s_k_2 < 0) {
			A1 = g(Poly<T>{-1.0, 1.0},  -s_k_2);
			A3 = h(Poly<T>{1.0, -1.0},  -s_k_2);
		} else if(std::abs(s_k_2) < EPS) { 
			zero_A3 = true;
			A3 = zero_poly;
			A1 = f(Poly<T>{1.0, -1.0-s_k_2}, -1.0 - s_k_2);
		} else {
			fl_A3 = true;
			A3 = f(Poly<T>{1.0, -1.0+s_k_2},  1.0 + s_k_2);
			A1 = f(Poly<T>{1.0, -1.0-s_k_2}, -1.0 - s_k_2);
		}

		T from = 0.0;
		T to = 1.0;

		Poly<T> den_all = A0 + A1 + A2 + A3;
		Poly<T> x_all = A0*cv_x_0 + A1*cv_x_1 + A2*cv_x_2 + A3*cv_x_3;
        Poly<T> y_all = A0*cv_y_0 + A1*cv_y_1 + A2*cv_y_2 + A3*cv_y_3;


		Poly<T> den = A1 + A2;
		Poly<T> x =  A1*cv_x_1 + A2*cv_x_2;
        Poly<T> y =  A1*cv_y_1 + A2*cv_y_2;



		if (zero_A0 && zero_A3) {
			return {SubSegment<T>(den, x, y, from, to)};	
        }
        
        if (zero_A3) {
            if (fl_A0 && std::abs(s_k_1-1) > EPS) {
                    return {SubSegment<T>(den_all - A3, x_all - A3*cv_x_3, y_all - A3*cv_x_3, 0.0, s_k_1), 
                            SubSegment<T>(den, x, y, s_k_1, 1.0)};
            }
            else {
                return {SubSegment<T>(den_all - A3, x_all - A3*cv_x_3, y_all - A3*cv_y_3, from, to)};
            }
        }

        if (zero_A0) {
            if (fl_A3 && std::abs(s_k_2-1) > EPS) {
                    return {SubSegment<T>(den, x, y, 0.0, 1.0-s_k_2), 
                            SubSegment<T>(den_all - A0, x_all - A0*cv_x_0, y_all - A0*cv_y_0, 1.0-s_k_2, 1.0)};
                }
            else {
                return {SubSegment<T>(den_all - A0, x_all - A0*cv_x_0, y_all - A0*cv_y_0, from, to)};        
            }
        }

        // if (zero_A3) {
        //     if (fl_A0) {
        //         if (std::abs(s_k_1-1) > EPS) {
        //             return {SubSegment<T>(den_all - A3, x_all - A3*cv_x_3, y_all - A3*cv_x_3, 0.0, s_k_1), 
        //                     SubSegment<T>(den, x, y, s_k_1, 1.0)};
        //         }
        //         else { // [0, 1]
        //             return {SubSegment<T>(den_all - A3, x_all - A3*cv_x_3, y_all - A3*cv_y_3, 0.0, s_k_1)};
        //         }
        //     } else {
        //         return {SubSegment<T>(den_all - A3, x_all - A3*cv_x_3, y_all - A3*cv_y_3, from, to)};
        //     }
        // }

        // if (zero_A0) {
        //     if (fl_A3) {
        //         if(std::abs(s_k_2-1) > EPS) {
        //             return {SubSegment<T>(den, x, y, 0.0, 1.0-s_k_2), 
        //                     SubSegment<T>(den_all - A0, x_all - A0*cv_x_0, y_all - A0*cv_y_0, 1.0-s_k_2, 1.0)};
        //         } else { // [0, 1]
        //             return {SubSegment<T>(den_all - A0, x_all - A0*cv_x_0, y_all - A0*cv_y_0, from, to)};
        //         }
        //     }
        //     else {
        //         return {SubSegment<T>(den_all - A0, x_all - A0*cv_x_0, y_all - A0*cv_y_0, from, to)};        
        //     }
        // }

        if(!fl_A0 && !fl_A3) {
            return {SubSegment<T>(den_all, x_all, y_all, from, to)};
        }
        
        if (fl_A0 && fl_A3) {
            // [0.0, s_k_1 == 1 - s_k_2 , 1.0]
            if (std::abs(s_k_1 + s_k_2 - 1.0) < EPS) {
                return {SubSegment<T>(den_all - A3, x_all - A3*cv_x_3, y_all - A3*cv_y_3, from, s_k_1), 
                        SubSegment<T>(den_all - A0, x_all - A0*cv_x_0, y_all - A0*cv_y_0, T(1.0)-s_k_2, to)};
            }
            else  {
                if (s_k_1 < 1-s_k_2) { // [0, s_k_1, 1-s_k_2, 1.0]
                    return {SubSegment<T>(den_all - A3, x_all - A3*cv_x_3, y_all - A3*cv_y_3, from, s_k_1), 
                            SubSegment<T>(den, x, y, s_k_1, 1.0-s_k_2), 
                            SubSegment<T>(den_all - A0, x_all - A0*cv_x_0, y_all - A0*cv_y_0, T(1.0)-s_k_2, to)}; 

                } else {
                    //пересечение [0, 1-s_k_2, s_k_1, 1.0]
                    std::vector<SubSegment<T>> res;
                    if (std::abs(s_k_2 - 1.0) > EPS) 
                        res.push_back(SubSegment<T>(den_all - A3, x_all - A3*cv_x_3, y_all - A3*cv_y_3, from, T(1.0)-s_k_2));
                    res.push_back(SubSegment<T>(den_all, x_all, y_all, 1-s_k_2, s_k_1));
                    if (std::abs(s_k_1 - 1.0) > EPS)
                        res.push_back(SubSegment<T>(den_all - A0, x_all - A0*cv_x_0, y_all - A0*cv_y_0, s_k_1, to));
                    return res;
                }
            }
        }      
        if (fl_A0) {
            if (std::abs(s_k_1-1.0) > EPS) {
                return {SubSegment<T>(den_all, x_all, y_all, from, s_k_1), 
                        SubSegment<T>(den_all - A0, x_all - A0*cv_x_0, y_all - A0*cv_y_0, s_k_1, to)};
            } else { 
                return {SubSegment<T>(den_all, x_all, y_all, from, s_k_1)};
            }
        }
            
        if (fl_A3) {
            if (std::abs(s_k_2-1.0) > EPS) {
                return {SubSegment<T>(den_all - A3, x_all - A3*cv_x_3, y_all - A3*cv_y_3, from, T(1.0)-s_k_2),
                        SubSegment<T>(den_all, x_all, y_all, T(1.0)-s_k_2, to)};
            }
            else {
                return {SubSegment<T>(den_all, x_all, y_all, from, to)};
            }
        }

		return {SubSegment<T>(den_all, x_all, y_all, from, to)};
	}

    // Получение N  точек на j  сегменте
    std::pair<std::vector<T>, std::vector<T>> get_points_per_segment(int j, int N = 1000) {
        if (j < 0 || j > N_segments - 1)
            throw std::invalid_argument("Wrong number(" + std::to_string(j) + 
                ") segments. Right : [0, " + std::to_string(N_segments - 1) + "].");

        std::vector<T> x_segm;
        std::vector<T> y_segm;
        int n_subsegm = vv_segments[j].size();
        int N_acc = 0;
        int NN = 0;
        int k = 0;
        
        for (auto subsegm : vv_segments[j]) {
            if (k == (n_subsegm - 1)) NN = N - N_acc;
            else NN = int((subsegm.to - subsegm.from) * N); 

            N_acc += NN;
            k += 1;

            std::vector<T> ts = linspace<T>(subsegm.from, subsegm.to, NN );
            std::vector<T> xs_temp(NN);
            std::vector<T> ys_temp(NN);
            for (int i = 0; i < NN; ++i) {
                xs_temp[i] = subsegm.x.At(ts[i]) / subsegm.den.At(ts[i]);
                ys_temp[i] = subsegm.y.At(ts[i]) / subsegm.den.At(ts[i]);
            }
            x_segm.insert(x_segm.end(), xs_temp.begin(), xs_temp.end());
            y_segm.insert(y_segm.end(), ys_temp.begin(), ys_temp.end());
        }
        return {x_segm, y_segm};
    }

    // Получение N_per_segment точек для каждого сегмента кривой.
    std::pair<std::vector<T>, std::vector<T>> get_points(int N_per_segment = 1000) {
        std::vector<T> xs;
        std::vector<T> ys;
        for (int i = 0; i < N_segments; ++i) {
            auto res_x_y_segm = get_points_per_segment(i, N_per_segment);
            xs.insert(xs.end(), res_x_y_segm.first.begin(), res_x_y_segm.first.end());
            ys.insert(ys.end(), res_x_y_segm.second.begin(), res_x_y_segm.second.end());
        }
        return {xs, ys};
    }

    // Печать в консоль коэффициентов рациональных функций сегментов
    void print_segments() {
        int i = 0;
        for (auto segm : vv_segments) {
            std::cout << "============_" + std::to_string(i) + "_============\n";
            for (auto subsegm : segm) {
                subsegm.print();
                std::cout << "\n";
            }
            ++i;
        }
    }
    

    // Метод отрисовки X сплайн кривой по-сегментно.
    // path -- путь к сохраненному файлу формата .png
    void plot_segments(std::string path="segments.png") {
        std::cout << "Number of segments : " << N_segments << "\n";
        
        plt::xlabel("X", {{"fontsize", "large"}});
        plt::ylabel("Y", {{"fontsize", "large"}});
        
        for (int i = 0; i < N_segments - 1; ++i) {
            auto res = get_points_per_segment(i);
            plot_curve_(res.first, res.second);
        }
        auto res = get_points_per_segment(N_segments-1);
        plot_curve(cv_x, cv_y, res.first, res.second, 
                    path);
    }

    // сохранение вычисленных  N точек кривой в файлы. 
    void save_points(int N = 1000) {
        auto res = get_points(N);
        save_vector(res.first,  "x.txt");
        save_vector(res.second, "y.txt");
    }
    
    //  Нахождение интеграла численно по квадратуре Гаусса-Лежандра.
    std::complex<T> numerical_integral() {
        static const std::vector<long double> ps = POINTS;
        static const std::vector<long double> ws = WEIGHTS;
        int n = ps.size();
        T A, B;
        T sum = 0.0;
        
        for (auto segm : vv_segments) {
            for (auto subsegm : segm) {
                T A = subsegm.from;
                T B = subsegm.to;
                T sum_temp = 0.0;
                
                for (int j = 0; j < n; ++j) {
                    T x = T(B-A)*ps[j]/2 + T(A+B)/2;
                    T den_quadratic = std::pow(subsegm.den.At(x),2);
                    sum_temp += ws[j] * subsegm.y.At(x) * subsegm.x.der().At(x) / den_quadratic;
                    sum_temp -= ws[j] * subsegm.y.At(x) * subsegm.den.der().At(x)*subsegm.x.At(x) / den_quadratic / subsegm.den.At(x); //
                }
                sum += (B-A)/2.0 *sum_temp;
            }
        }
        return sum;
    }

    // Нахождение интеграла по Алгоритму 2.
    std::complex<T> analytical_integral(int type = 1) {
        std::complex<T> sum = 0.0;
        for (auto segm : vv_segments) {
            for (auto subsegm : segm) {
                
                T from =  subsegm.from;
                T to = subsegm.to;

                Poly<T> p_x, p_x_der, p_y, den, num;
                
                den = subsegm.den;            
                p_x = subsegm.x;
                p_y = subsegm.y;
                p_x_der = subsegm.x.der();

                std::vector<std::pair<int, std::complex<T>>> roots = den.solve(type);
                int n_roots = roots.size();
                std::vector<std::complex<T>> vec_root(n_roots);
                for (int j = 0; j < n_roots; ++j) {
                    vec_root[j] = roots[j].second;
                }
                std::complex<T> temp_sum = 0.0;
                
                auto div1 = p_y / den;
                auto div2 = p_x_der / den;
                
                Poly<T> P_12 = div1[0] * div2[0];
                temp_sum +=  P_12.integral(from, to) ;
                num = div1[0] * div2[1] + div2[0] * div1[1];
                auto div3 = num / den;    
                temp_sum += div3[0].integral(from, to) ;
                temp_sum += integral(div3[1], vec_root, from, to) ;
                temp_sum += integral(div1[1], div2[1], vec_root, from, to) ;
                sum += temp_sum;


                temp_sum = 0.0;
                div2 =  p_x / den;
                P_12 = div1[0] * div2[0];
                Poly<std::complex<T>> P_12_complex  = make_complex_poly<T>(P_12);
                num = div1[0] * div2[1] + div2[0] * div1[1];
                auto div4 = num / den;
                Poly<std::complex<T>> P4 = make_complex_poly<T>(div4[0]);
                for (int j = 0; j < n_roots; ++j) {
                    Poly<std::complex<T>> divider(std::vector<std::complex<T>>{1.0, -vec_root[j]});
                    auto div3_complex = P_12_complex / divider;    
                    temp_sum +=  integral_type_1(div3_complex[1][0], vec_root[j], 1, from, to) ;
                    temp_sum += div3_complex[0].integral(from, to) ;
                    auto div5 = P4 / divider;
                    temp_sum += div5[0].integral(from, to) ;
                    temp_sum += integral_type_1(div5[1][0], vec_root[j], 1, from, to) ;
                    temp_sum += integral(div4[1], vec_root, j, from, to) ;         
                    temp_sum += integral(div1[1], div2[1], vec_root, j, from, to) ;   
                }    
                sum -= temp_sum;
            }
        }
        return sum;
    }
};