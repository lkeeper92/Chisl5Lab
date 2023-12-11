#include <iostream>
#include <cmath>
#include <string>

using namespace std;

double function(double x) {
    return 1/(pow(log10(x),2)+1);
}
double I1(double a, double b, int n) {
    double _In = 0;
    double h = (b - a) / n;
    for (double x = a; x < b; x += h) {
        _In += function(x);
    }
    return _In;
}
double I2(double a, double b, int n, int temp) {
    double _In = 0;
    double h = (b - a) / n;
    for (double x = a + temp * h; x < b; x += 2 * h) {
        _In += function(x);
    }
    return _In;
}
double trapIntegral(double a, double b, double E) {
    double e = 1, h;
    double In1 = 0, In2 = 0;
    int n = 2;
    while (e > 3 * E) {
        h = (b - a) / n;
        In1 = (h / 2) * (function(a) + 2 * I1(a, b, n) + function(b));
        In2 = (h / 4) * (function(a) + 2 * I1(a, b, n * 2) + function(b));
        e = abs(In2 - In1);
        n *= 2;
    }
    return In2;
}
double simpsonIntegral(double a, double b, double E) {
    double e = 1, h;
    double In1, In2 = 0;
    int n = 2;
    while (e > 15 * E) {
        h = (b - a) / n;
        In1 = (h / 3) * (function(a) + 4 * I2(a, b, n, 0) + 2 * I2(a, b + h, n, 1) + function(b));
        In2 = (h / 6) * (function(a) + 4 * I2(a, b, n * 2, 0) + 2 * I2(a, b + h, n * 2, 1) + function(b));
        e = abs(In2 - In1);
        n *= 2;
    }
    return In2;
}
double functionXY(double x, double y) {
    return (1 / ((x + y) * (x + y)));
}

double I(int i, double a, double hx, int j, double c, double hy) {
    return functionXY(a + i * hx, c + j * hy);
}
double simpsonIntegral(double a, double b, double c, double d, double E) {
    double e = 1, hx, hy;
    double In1 = 0, In2 = 0;
    double temp = 0;
    int n = 2, m = 2;
    while (e > E) {
        hx = (b - a) / (2 * n);
        hy = (d - c) / (2 * m);
        for (int i = 0; i <= n - 1; i++)
        {
            for (int j = 0; j <= m - 1; j++)
            {
                In1 += I(2 * i, a, hx, 2 * j, c, hy) + 4 * I(2 * i + 1, a, hx, 2 * j, c, hy) + I(2 * i + 2, a, hx, 2 * j, c, hy) + 4 * I(2 * i, a, hx, 2 * j + 1, c, hy)
                    + 16 * I(2 * i + 1, a, hx, 2 * j + 1, c, hy) + 4 * I(2 * i + 2, a, hx, 2 * j + 1, c, hy) + I(2 * i, a, hx, 2 * j + 2, c, hy) + 4 * I(2 * i + 1, a, hx, 2 * j + 2, c, hy)
                    + I(2 * i + 2, a, hx, 2 * j + 2, c, hy);
            }
            In2 += In1;
            In1 = 0;
        }
        In2 *= (hx * hy) / 9;
        e = abs(In2 - temp);
        temp = In2;
        n *= 2;
    }
    return In2;
}

int main() {
    int n;
    double e1 = 1e-4, e2 = 1e-5;
    double I1 = trapIntegral(1.0, 2.835, e1);
    double I2 = simpsonIntegral(1.0, 2.835, e1);
    double I3 = trapIntegral(1.0, 2.835, e2);
    double I4 = simpsonIntegral(1.0, 2.835, e2);
    cout << "Epsilon = 1e-4: \n Trapezoid method: " << I1 << "\n Simpson method: " << I2 << endl;
    cout << "=========================================\n";
    cout << "Epsilon = 1e-5: \n Trapezoid method: " << I3 << "\n Simpson method: " << I4 << endl;
    cout << "=========================================\n";
    cout <<"Epsilon = 1e-4: \nSimpson's method for function of two variables: " << simpsonIntegral(3.0, 4.0, 1.0, 2.0, e1)<<endl;
    cout << "=========================================\n";
    cout << "Epsilon = 1e-5: \nSimpson's method for function of two variables: " << simpsonIntegral(3.0, 4.0, 1.0, 2.0, e2) << endl;
    return 0;
}

