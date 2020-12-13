#define _USE_MATH_DEFINES

#include "matplotlibcpp.h"

#include <iostream>
#include <cmath>
#include <functional>

struct StewartData
{
    double L1, L2, L3, gamma, p1, p2, p3, x1, x2, y2;
    StewartData(const double &fl1,
                const double &fl2,
                const double &fl3,
                const double &fgamma,
                const double &fp1,
                const double &fp2,
                const double &fp3,
                const double &fx1,
                const double &fx2,
                const double &fy2)
        : L1(fl1), L2(fl2), L3(fl3), gamma(fgamma),
          p1(fp1), p2(fp2), p3(fp3), x1(fx1), x2(fx2), y2(fy2) {}
};

inline double power(const double &x)
{
    return x * x;
}

double getSolution(double x0,
                   double x1,
                   std::function<double(const double &)> function,
                   const double &epsilon = 0.0000001,
                   const int &max_iteration = 100)
{
    double res = 0;
    for (int i = 0; i < max_iteration; i++)
    {
        res = x1 - (function(x1) * (x1 - x0)) / (function(x1) - function(x0));
        if (std::abs(function(res)) <= epsilon || std::abs(res - x1) <= epsilon)
        {
            break;
        }

        x0 = x1;
        x1 = res;
    }

    return res;
}

int main()
{
    const double SQRT2 = std::sqrt(2);
    const double PI = std::acos(-1);
    const double SQRT5 = std::sqrt(5);
    constexpr size_t NUMBER_OF_DATA = 2;

    namespace plt = matplotlibcpp;

    StewartData data[NUMBER_OF_DATA] = {StewartData(2, SQRT2, SQRT2, PI / 2, SQRT5, SQRT5, SQRT5, 4, 0, 4),
                                        StewartData(3, 3 * SQRT2, 3, PI / 4, 5, 5, 3, 5, 0, 6)};

    for (int i = 0; i < NUMBER_OF_DATA; i++)
    {
        std::cout << "--------- "
                  << "data " << i << " ---------" << std::endl;

        auto p = data[i];
        auto A2 = [p](const double &theta) -> double {
            return p.L3 * std::cos(theta) - p.x1;
        };
        auto B2 = [p](const double &theta) -> double {
            return p.L3 * std::sin(theta);
        };
        auto A3 = [p](const double &theta) -> double {
            return p.L2 * std::cos(theta + p.gamma) - p.x2;
        };
        auto B3 = [p](const double &theta) -> double {
            return p.L2 * std::sin(theta + p.gamma) - p.y2;
        };
        auto D = [p, A2, A3, B2, B3](const double &theta) -> double {
            return 2 * (A2(theta) * B3(theta) - B2(theta) * A3(theta));
        };
        auto N1 = [p, A2, A3, B2, B3](const double &theta) -> double {
            return B3(theta) * (power(p.p2) - power(p.p1) - power(A2(theta)) - power(B2(theta))) - B2(theta) * (power(p.p3) - power(p.p1) - power(A3(theta)) - power(B3(theta)));
        };
        auto N2 = [p, A2, A3, B2, B3](const double &theta) -> double {
            return -A3(theta) * (power(p.p2) - power(p.p1) - power(A2(theta)) - power(B2(theta))) + A2(theta) * (power(p.p3) - power(p.p1) - power(A3(theta)) - power(B3(theta)));
        };
        auto f = [N1, N2, D, p](const double &theta) -> double {
            return power(N1(theta)) + power(N2(theta)) - power(p.p1) * power(D(theta));
        };

        std::vector<double> curve_x, curve_y;
        for (double i = -1; i < 2.2; i += 0.01)
        {
            curve_x.push_back(i);
            curve_y.push_back(f(i));
        }

        plt::plot(curve_x, curve_y);

        if (i == 0)
        {
            double theta1 = getSolution(-1.0, -0.9, f, 1e-12),
                   theta2 = getSolution(1.0, 0.9, f, 1e-12);

            double x1 = N1(theta1) / D(theta1), y1 = N2(theta1) / D(theta1);
            double x2 = N1(theta2) / D(theta2), y2 = N2(theta2) / D(theta2);

            std::cout << "theta: " << theta1 << ", x: " << x1 << ", y: " << y1 << std::endl;
            std::cout << "theta: " << theta2 << ", x: " << x2 << ", y: " << y2 << std::endl;
        }
        else
        {
            double theta1 = getSolution(-1.0, -0.9, f, 1e-12),
                   theta2 = getSolution(-0.4, -0.5, f, 1e-12),
                   theta3 = getSolution(0.9, 1.0, f, 1e-12),
                   theta4 = getSolution(1.8, 1.9, f, 1e-12);

            double x1 = N1(theta1) / D(theta1), y1 = N2(theta1) / D(theta1);
            double x2 = N1(theta2) / D(theta2), y2 = N2(theta2) / D(theta2);
            double x3 = N1(theta3) / D(theta3), y3 = N2(theta3) / D(theta3);
            double x4 = N1(theta4) / D(theta4), y4 = N2(theta4) / D(theta4);

            std::cout << "theta: " << theta1 << ", x: " << x1 << ", y: " << y1 << std::endl;
            std::cout << "theta: " << theta2 << ", x: " << x2 << ", y: " << y2 << std::endl;
            std::cout << "theta: " << theta3 << ", x: " << x3 << ", y: " << y3 << std::endl;
            std::cout << "theta: " << theta4 << ", x: " << x4 << ", y: " << y4 << std::endl;
        }
    }

    std::vector<double> curve_x, curve_y;
    for (double i = -1; i < 2.2; i += 0.01)
    {
        curve_x.push_back(i);
        curve_y.push_back(0);
    }

    plt::plot(curve_x, curve_y);

    plt::save("curve.png");

    return 0;
}
