//
// Created by a_d_e on 17.04.2021.
//
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

#define         HIGH_PRECISION    0.1
#define         LOW_PRECISION     1.0
#define         ZOOM_PRECISION    5.0
#define         MAX_SPLINE_STEP   0.01

class XSpline {
public:
    XSpline(const double *control_x, const double *control_y, int count, const double *shape) {
        this->count = count;
        this->control_x = new double[count + 2];
        this->control_y = new double[count + 2];
        this->shape = new double[count + 2];
        this->control_x[0] = control_x[0];
        this->control_y[0] = control_y[0];
        this->control_x[count + 1] = control_x[count - 1];
        this->control_y[count + 1] = control_y[count - 1];
        this->shape[0] = 0;
        this->shape[1] = 0;
        this->shape[count] = 0;

        for (int i = 0; i < this->count; ++i) {
            this->control_x[i + 1] = control_x[i];
            this->control_y[i + 1] = control_y[i];
            this->shape[i + 2] = shape[i];
        }
    }

    void compute(double precision, const string &filename) {
        out.open(filename);
        out << setprecision(6) << fixed;
        compute_open_spline(precision);
        out.close();
    }

    virtual ~XSpline() {
        delete[] control_x;
        delete[] control_y;
        delete[] shape;
    }

private:
    double *control_x;
    double *control_y;
    double *shape;
    int count;
    ofstream out;

    static double f_blend(double numerator, double denominator) {
        double p = 2 * denominator * denominator;
        double u = numerator / denominator;
        double u2 = u * u;

        return (u * u2 * (10 - p + (2 * p - 15) * u + (6 - p) * u2));
    }

    static double g_blend(double u, double q) {/* p equals 2 */
        return (u * (q + u * (2 * q + u * (8 - 12 * q + u * (14 * q - 11 + u * (4 - 5 * q))))));
    }

    static double h_blend(double u, double q) {
        double u2 = u * u;
        return (u * (q + u * (2 * q + u2 * (-2 * q - u * q))));
    }

    void negative_s1_influence(double t, int s1, double *A0, double *A2) {
        *A0 = h_blend(-t, -shape[s1]);
        *A2 = g_blend(t, -shape[s1]);
    }

    void negative_s2_influence(double t, int s2, double *A1, double *A3) {
        *A1 = g_blend(1 - t, -shape[s2]);
        *A3 = h_blend(t - 1, -shape[s2]);
    }

    void positive_s1_influence(int k, double t, int s1, double *A0, double *A2) {
        double Tk;

        Tk = k + 1 + shape[s1];
        *A0 = (t + k + 1 < Tk) ? f_blend(t + k + 1 - Tk, k - Tk) : 0.0;

        Tk = k + 1 - shape[s1];
        *A2 = f_blend(t + k + 1 - Tk, k + 2 - Tk);
    }

    void positive_s2_influence(int k, double t, int s2, double *A1, double *A3) {
        double Tk;

        Tk = k + 2 + shape[s2];
        *A1 = f_blend(t + k + 1 - Tk, k + 1 - Tk);

        Tk = k + 2 - shape[s2];
        *A3 = (t + k + 1 > Tk) ? f_blend(t + k + 1 - Tk, k + 3 - Tk) : 0.0;
    }

    void add_point(double x, double y) {
#pragma omp critical
        out << "x: " << x << ", y: " << y << endl;
    }

    void point_adding(const double *blend, int p0, int p1, int p2, int p3) {
        double weights_sum;

        weights_sum = blend[0] + blend[1] + blend[2] + blend[3];
        add_point((control_x[p0] * blend[0] + control_x[p1] * blend[1] + control_x[p2] * blend[2] +
                   control_x[p3] * blend[3]) / weights_sum,
                  (control_y[p0] * blend[0] + control_y[p1] * blend[1] + control_y[p2] * blend[2] +
                   control_y[p3] * blend[3]) / weights_sum);
    }

    void point_computing(const double *blend, int p0, int p1, int p2, int p3, double *x, double *y) {
        double weights_sum;

        weights_sum = blend[0] + blend[1] + blend[2] + blend[3];

        *x = (control_x[p0] * blend[0] + control_x[p1] * blend[1] + control_x[p2] * blend[2] +
              control_x[p3] * blend[3]) / weights_sum;
        *y = (control_y[p0] * blend[0] + control_y[p1] * blend[1] + control_y[p2] * blend[2] +
              control_y[p3] * blend[3]) / weights_sum;
    }

    double step_computing(int k, int p0, int p1, int p2, int p3, int s1, int s2, double precision) {
        double blend[4];
        double x_start, y_start, x_end, y_end, x_mid, y_mid, x_length, y_length, start_to_end_dist, number_of_steps;
        double step, angle_cos, scale_prod, xv1, xv2, yv1, yv2, sides_length_prod;

        /* This function computes the step used to draw the segment (p1, p2)
           (xv1, yv1) : coordinates of the vector from middle to origin
           (xv2, yv2) : coordinates of the vector from middle to extremity */

        if ((shape[s1] == 0) and (shape[s2] == 0))
            return (1.0);              /* only one step in case of linear segment */

        /* compute coordinates of the origin */
        if (shape[s1] > 0) {
            if (shape[s2] < 0) {
                positive_s1_influence(k, 0.0, s1, &blend[0], &blend[2]);
                negative_s2_influence(0.0, s2, &blend[1], &blend[3]);
            } else {
                positive_s1_influence(k, 0.0, s1, &blend[0], &blend[2]);
                positive_s2_influence(k, 0.0, s2, &blend[1], &blend[3]);
            }
            point_computing(blend, p0, p1, p2, p3, &x_start, &y_start);
        } else {
            x_start = control_x[p1];
            y_start = control_y[p1];
        }

        /* compute coordinates  of the extremity */
        if (s2 > 0) {
            if (s1 < 0) {
                negative_s1_influence(1.0, s1, &blend[0], &blend[2]);
                positive_s2_influence(k, 1.0, s2, &blend[1], &blend[3]);
            } else {
                positive_s1_influence(k, 1.0, s1, &blend[0], &blend[2]);
                positive_s2_influence(k, 1.0, s2, &blend[1], &blend[3]);
            }
            point_computing(blend, p0, p1, p2, p3, &x_end, &y_end);
        } else {
            x_end = control_x[p2];
            y_end = control_y[p2];
        }

        /* compute coordinates  of the middle */
        if (s2 > 0) {
            if (s1 < 0) {
                negative_s1_influence(0.5, s1, &blend[0], &blend[2]);
                positive_s2_influence(k, 0.5, s2, &blend[1], &blend[3]);
            } else {
                positive_s1_influence(k, 0.5, s1, &blend[0], &blend[2]);
                positive_s2_influence(k, 0.5, s2, &blend[1], &blend[3]);
            }
        } else if (s1 < 0) {
            negative_s1_influence(0.5, s1, &blend[0], &blend[2]);
            negative_s2_influence(0.5, s2, &blend[1], &blend[3]);
        } else {
            positive_s1_influence(k, 0.5, s1, &blend[0], &blend[2]);
            negative_s2_influence(0.5, s2, &blend[1], &blend[3]);
        }
        point_computing(blend, p0, p1, p2, p3, &x_mid, &y_mid);

        xv1 = x_start - x_mid;
        yv1 = y_start - y_mid;
        xv2 = x_end - x_mid;
        yv2 = y_end - y_mid;

        scale_prod = xv1 * xv2 + yv1 * yv2;

        sides_length_prod = sqrt((xv1 * xv1 + yv1 * yv1) * (xv2 * xv2 + yv2 * yv2));

        if (sides_length_prod == 0.0)
            angle_cos = 0.0;
        else
            angle_cos = scale_prod / sides_length_prod;

        x_length = x_end - x_start;
        y_length = y_end - y_start;

        start_to_end_dist = sqrt((double) x_length * (double) x_length + (double) y_length * (double) y_length);

        number_of_steps = sqrt(start_to_end_dist) / 2;

        number_of_steps += (int) ((1 + angle_cos) * 10);

        if (number_of_steps == 0)
            step = 1;
        else
            step = precision / number_of_steps;

        if ((step > MAX_SPLINE_STEP) || (step == 0))
            step = MAX_SPLINE_STEP;
        return step;
    }

    void spline_segment_computing(double step, int k, int p0, int p1, int p2, int p3, int s1, int s2) {
        double blend[4];

        if (s1 < 0) {
            if (s2 < 0) {
                for (double t = 0.0; t < 1; t += step) {
                    negative_s1_influence(t, s1, &blend[0], &blend[2]);
                    negative_s2_influence(t, s2, &blend[1], &blend[3]);

                    point_adding(blend, p0, p1, p2, p3);
                }
            } else {
                for (double t = 0.0; t < 1; t += step) {
                    negative_s1_influence(t, s1, &blend[0], &blend[2]);
                    positive_s2_influence(k, t, s2, &blend[1], &blend[3]);

                    point_adding(blend, p0, p1, p2, p3);
                }
            }
        } else if (s2 < 0) {
            for (double t = 0.0; t < 1; t += step) {
                positive_s1_influence(k, t, s1, &blend[0], &blend[2]);
                negative_s2_influence(t, s2, &blend[1], &blend[3]);

                point_adding(blend, p0, p1, p2, p3);
            }
        } else {
            for (double t = 0.0; t < 1; t += step) {
                positive_s1_influence(k, t, s1, &blend[0], &blend[2]);
                positive_s2_influence(k, t, s2, &blend[1], &blend[3]);

                point_adding(blend, p0, p1, p2, p3);
            }
        }
    }

    void compute_open_spline(double precision) {
#pragma omp parallel for default(none) shared(precision)
        for (int i = 0; i < count; i++) {
            double step = step_computing(i, i, i + 1, i + 2, i + 3, i + 1, i + 2, precision);
            spline_segment_computing(step, i, i, i + 1, i + 2, i + 3, i + 1, i + 2);
        }
    }
};

int main() {
    double x[] = {1, 3, 5, 7, 9};
    double y[] = {4, 6, 4, 6, 4};
    double s[] = {1, 1, 1, 1};
    XSpline spline(x, y, 5, s);
    spline.compute(0.01, "result.txt");
    return 0;
}