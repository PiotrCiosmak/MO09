#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

using namespace std;

constexpr double p{1.0}, q{}, r{-1.0 / 4.0}, s{};
constexpr double alfa{}, beta_{1.0}, gamma{-1.0}, theta{}, fi{}, psi{1.0};
constexpr double xBegin{}, xEnd{2.0};

double explication(const double &x);

void thomasAlgorithm(const vector<double> &L, const vector<double> &D, const vector<double> &U, const vector<double> &B, vector<double> &X, const int &n);

double maxError(const vector<double> &error, const int &n);

double numberMethod(const double &h, const int &n);

double conventionalMethod(const double &h, const int &n);

int main()
{
    double h, conventionalError, numberError;
    fstream bledyNumeryczne("bledyNumeryczne.txt", fstream::out);
    fstream bledyKonwencjonalne("bledyKonwencjonalne.txt", fstream::out);

    for (int N = 10; N < 20000; N += 10)
    {
        h = (xEnd - xBegin) / (N - 1);
        conventionalError = log10(conventionalMethod(h, N));
        numberError = log10(numberMethod(h, N));

        bledyNumeryczne << log10(h) << " " << numberError << endl;
        bledyKonwencjonalne << log10(h) << " " << conventionalError << endl;
    }

    bledyKonwencjonalne.close();
    bledyNumeryczne.close();
}

double explication(const double &x)
{
    return (exp(x / 2.0) - exp(2.0 - x / 2.0)) / (1.0 - exp(2.0));
}

void thomasAlgorithm(const vector<double> &L, const vector<double> &D, const vector<double> &U, const vector<double> &B, vector<double> &X, const int &n)
{
    vector<double> tmpB(n, 0);
    vector<double> tmpD(n, 0);

    tmpD.at(0) = D.at(0);
    tmpB.at(0) = B.at(0);

    for (int i = 1; i < n; i++)
        tmpD.at(i) = D.at(i) - L.at(i - 1) * (U.at(i - 1) / tmpD.at(i - 1));

    for (int i = 1; i < n; i++)
        tmpB.at(i) = B.at(i) - L.at(i - 1) * tmpB.at(i - 1) / tmpD.at(i - 1);

    X.at(n - 1) = tmpB.at(n - 1) / tmpD.at(n - 1);

    for (int i = n - 2; i >= 0; i--)
        X.at(i) = (tmpB.at(i) - U.at(i) * X.at(i + 1)) / tmpD.at(i);
}

double maxError(const vector<double> &error, const int &n)
{
    double max = fabs(error.at(0));

    for (int i = 0; i < n; ++i)
        if (fabs(error.at(i)) > max)
            max = fabs(error.at(i));

    return max;
}

double numberMethod(const double &h, const int &n)
{
    double x0{xBegin}, x1{xBegin};

    vector<double> B(n, 0);
    vector<double> D(n, 0);
    vector<double> L(n, 0);
    vector<double> U(n, 0);
    vector<double> X(n, 0);
    vector<double> error(n, 0);

    U.at(0) = alfa / h;
    D.at(0) = beta_ - alfa / h;
    B.at(0) = -gamma;

    for (int i = 1; i < n - 1; ++i)
    {
        L.at(i - 1) = U.at(i) = 1.0 / (h * h) + (-1.0 / 4.0) / 12.0;
        D.at(i) = (-2.0) / (h * h) + (-1.0 / 4.0) * 10.0 / 12.0;
        B.at(i) = -s;
    }

    L.at(n - 2) = -fi / h;
    D.at(n - 1) = -fi / h + psi;
    B.at(n - 1) = -theta;

    thomasAlgorithm(L, D, U, B, X, n);

    for (int i = 0; i < n; ++i, x0 += h)
        error.at(i) = fabs(X.at(i) - explication(x0));

    if (n == 20)
    {
        fstream NumeralResult("NumeralResult.txt", fstream::out);
        fstream exactResult("exactResult.txt", fstream::out);

        for (int i = 0; i < n; ++i, x1 += h)
        {
            NumeralResult << x1 << " " << X.at(i) << endl;
            exactResult << x1 << " " << explication(x1) << endl;
        }

        NumeralResult.close();
        exactResult.close();
    }

    return maxError(error, n);
}

double conventionalMethod(const double &h, const int &n)
{
    double x0{xBegin}, x1{xBegin};

    vector<double> B(n, 0);
    vector<double> D(n, 0);
    vector<double> L(n, 0);
    vector<double> U(n, 0);
    vector<double> X(n, 0);
    vector<double> error(n, 0);

    U.at(0) = alfa / h;
    D.at(0) = beta_ - alfa / h;
    B.at(0) = -gamma;

    for (int i = 1; i < n - 1; ++i)
    {
        L.at(i - 1) = U.at(i) = p / (h * h) - q / (2.0 * h);;
        D.at(i) = (-2.0 * p) / (h * h) + r;
        B.at(i) = -s;
    }

    L.at(n - 2) = -fi / h;
    D.at(n - 1) = -fi / h + psi;
    B.at(n - 1) = -theta;

    thomasAlgorithm(L, D, U, B, X, n);

    for (int i = 0; i < n; ++i, x0 += h)
        error.at(i) = fabs(X.at(i) - explication(x0));

    if (n == 20)
    {
        fstream ConventionalResult("ConventionalResult.txt", fstream::out);

        for (int i = 0; i < n; ++i, x1 += h)
        {
            error.at(i) = fabs(X.at(i) - explication(x1));
            ConventionalResult << x1 << " " << X.at(i) << endl;
        }

        ConventionalResult.close();
    }

    return maxError(error, n);
}