#pragma once
#include "matplotlibcpp.h"

#include "Point.h"

//Флаги компиляции: -I/usr/include/python3.10  -lpython3.10 

namespace plt = matplotlibcpp;


// Вспомогательные функции 


// After calling the rendering functions, you need to call the function to kill python_interpreter
//https://stackoverflow.com/questions/67533541/py-finalize-resulting-in-segmentation-fault-for-python-3-9-but-not-for-python/67577360#67577360
void PLOT_END() {
    plt::detail::_interpreter::kill();
}

// Настройка параметров рисунка.
void set_figure(std::string title="") {
    plt::title(title);
    std::map<std::string, std::string> keywords_fig;
    keywords_fig["savefig.bbox"] = "tight"; // узкие рамки.
    plt::rcparams(keywords_fig);
}


// Функции отрисовки графиков с сохранением рисунка в файл: 

template<typename T>
void plot_y_line_(T value, std::vector<int> xs) {
    std::vector<T> ys(xs.size(), value);
    plt::semilogy(xs, ys, "b--");
}

template<typename T>
void plot_errors_( std::vector<T> errors, std::string legend="") {
    int n = errors.size();
    std::vector<T> xs(n);
    for (int i = 1; i <= n; ++i) 
        xs[i-1] = T(i) / n;

    if (legend != "") {
        plt::named_semilogy(legend, xs, errors);
    } else {
        plt::semilogy( xs, errors);
    }
}

    
template<typename T>
void plot_errors_(std::vector<int> xs, 
                std::vector<T> errors, std::string legend="") {
    if (legend != "") {
        plt::named_semilogy(legend, xs, errors);
    } else {
        plt::semilogy( xs, errors);
    }
}

template<typename T>
void plot_errors_(std::vector<T> xs, 
                std::vector<T> errors, std::string legend="") {
    if (legend != "") {
        plt::named_semilogy(legend, xs, errors);
    } else {
        plt::semilogy(xs, errors);
    }
}

template<typename T>
void plot_curve_( std::vector<T>& x, 
                std::vector<T>& y,
                const std::string& legend="") {
    
    std::map<std::string, std::string> keywords_points;
    if (legend == "") {
        plt::plot(x, y);
    } else {
        keywords_points["label"] = legend;
        keywords_points["linewidth"] = "2";
        plt::plot(x, y, keywords_points);
    }
}


template<typename T>
void plot_curve_( std::vector<T>& cv_x,
                std::vector<T>& cv_y,
                std::vector<T>& x, 
                std::vector<T>& y, 
                const std::string& legend="") {
    
    std::map<std::string, std::string> keywords_cv;
    keywords_cv["color"] = "green";
    keywords_cv["label"] = "Контрольные т.";

    plt::scatter(cv_x, cv_y, 20.0, keywords_cv);
    plot_curve_(x, y, legend);
}

template<typename T>
void plot_curve_(std::vector<Point<T>>& cv, 
            std::vector<Point<T>>& points, 
            const std::string& legend="") {
    
    size_t len = points.size();
    size_t len_cv = cv.size();
    
    std::vector<T> x(len);
    std::vector<T> y(len);
    
    std::vector<T> cv_x(len_cv);
    std::vector<T> cv_y(len_cv);
    
    for (int i=0; i < len; ++i) {
        x[i] = points[i][0];
        y[i] = points[i][1];
    }

    for (int i=0; i < len_cv; ++i) {
        cv_x[i] = cv[i][0];
        cv_y[i] = cv[i][1];
        
    }

    plot_curve_(cv_x, cv_y, x, y, legend);
}


template<typename T>
void plot_curve_(std::vector<Point<T>>& points,
                const std::string& legend="") {
    size_t len = points.size();

    std::vector<T> x(len);
    std::vector<T> y(len);

    for (int i=0; i < len; ++i) {
        x[i] = points[i][0];
        y[i] = points[i][1];
    }

    plot_curve_(x, y, legend);
}


// Функции отрисовки графиков с сохранением рисунка в файл:

template<typename T>
void plot_errors(std::vector<T> xs, 
                std::vector<T> errors,  
                const std::string& out="plot.png",
                const std::string& title="", 
                const std::string& legend="") {
    set_figure(title); 
    plot_errors_(xs, errors, legend);
    plt::legend();
    plt::save(out, 300);
    plt::close();  
}

template<typename T>
void plot_errors(std::vector<int> xs, 
                std::vector<T> errors,  
                const std::string& out="plot.png",
                const std::string& title="", 
                const std::string& legend="") {
    set_figure(title); 
    plot_errors_(xs, errors, legend);
    plt::legend();
    plt::save(out, 300);
    plt::close();  
}

template<typename T>
void plot_errors_scatter(std::vector<int> xs, 
                std::vector<T> errors,  
                const std::string& out="plot.png",
                const std::string& title="", 
                const std::string& legend="") {
    set_figure(title); 
    plot_errors_(xs, errors, legend);
    plt::legend();
    plt::save(out, 300);
    plt::close();  
}

template<typename T>
void plot_errors( std::vector<T> errors,  
                const std::string& out="plot.png",
                const std::string& title="",
                const std::string& legend="") {

    set_figure(title); 
    int n = errors.size();
    std::vector<T> xs(n);
    for (int i = 1; i <= n; ++i) 
        xs[i-1] = T(i) / n;
    
    plot_errors_(xs, errors, legend);
    plt::legend();
    plt::save(out, 300);
    plt::close();  
}

template<typename T>
void plot_curve( std::vector<T>& x, 
                std::vector<T>& y, 
                const std::string& out="plot.png",
                const std::string& title="",
                const std::string& legend="") {
    set_figure(title);                
    plot_curve_(x, y, legend);
    plt::legend();
    plt::save(out, 300);
    plt::close();
}

template<typename T>
void plot_curve( std::vector<T>& cv_x,
                std::vector<T>& cv_y,
                std::vector<T>& x, 
                std::vector<T>& y, 
                const std::string& out="plot.png",
                const std::string& title="",
                const std::string& legend="") {
    set_figure(title);                
    plot_curve_(cv_x, cv_y, x, y, legend);
    plt::legend();
    plt::save(out, 300);
    plt::close();
}

template<typename T>
void plot_curve(std::vector<Point<T>>& cv, 
            std::vector<Point<T>>& points, 
            const std::string& out="plot.png",
            const std::string& title="",
            const std::string& legend="") {
    set_figure(title); 
    plot_curve_(cv, points, legend);
    plt::legend();
    plt::save(out, 300);
    plt::close();
}

template<typename T>
void plot_curve(std::vector<Point<T>>& points, 
                const std::string& out="plot.png",
                const std::string& title="",
                const std::string& legend="") {
    set_figure(title);           
    plot_curve_(points, legend);
    std::map<std::string, std::string> m;
    m["fontsize"] = "medium";
    plt::legend(m);
    plt::save(out, 300);
    plt::close();
}