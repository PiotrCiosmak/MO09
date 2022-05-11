#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

using namespace std;

constexpr double alfa{}, beta_{1.0}, gamma{-1.0}, theta{};
constexpr double fi{}, psi{1.0};
constexpr double p{1.0}, q{}, r{1.0}, s{};

double explication(const double &x);

void
thomasAlgorithm(vector<double> &B, vector<double> &D, const vector<double> &L, const vector<double> &U,
                vector<double> &X, const int &n);

double calculateError(const vector<double> &X, const double &h, const int &n, double start);

void saveFile(const string &fileName, const vector<double> &X, const double &h, const int &n, double xStep);

double numberMethod(double h, int n, double x0, double x1);

double conventionalMethod(double h, int n, double x0, double x1);

int main()
{
    constexpr int start{10}, stop{1000}, step{10};
    constexpr double x0{}, x1{1.0};
    double h, numberError, conventionalError;

    try
    {
        fstream file("../blad.txt", ios::out);
        for (int i = start; i < stop; ++i)
        {
            h = (x1 - x0) / (i - 1.0);
            numberError = numberMethod(h, i, x0, x1);
            conventionalError = conventionalMethod(h, i, x0, x1);
            file << log10(h) << " " << log10(numberError) << " " << log10(conventionalError) << endl;
        }
        file.close();
    }
    catch (...)
    {
        cout << "Wystapil blad poczas zapisywania pliku blad.txt" << endl;
        exit(1);
    }

}

double explication(const double &x)
{
    return (exp(2.0 - 2.0 * x) - 4.0 * exp(4.0 - 2.0 * x) * 4.0 * exp(2.0 + 2.0 * x) - x * x * exp(4)) /
           (4.0 - 4.0 * exp(4));
}

void
thomasAlgorithm(vector<double> &B, vector<double> &D, const vector<double> &L, const vector<double> &U,
                vector<double> &X, const int &n)
{
    for (int i = 0; i < n; ++i)
    {
        D.at(i) = D.at(i) - (L.at(i - 1) * U.at(i - 1) / D.at(i - 1));
        B.at(i) = B.at(i) - (L.at(i - 1) * D.at(i - 1) / B.at(i - 1));
    }
    X.at(n - 1) = 1.0 / D.at(n - 1) * B.at(n - 1);

    for (int i = n - 2; i >= 0; --i)
        X.at(i) = 1.0 / D.at(i) * (B.at(i) - U.at(i) * X.at(i + 1));

}

double calculateError(const vector<double> &X, const double &h, const int &n, double start)
{
    vector<double> errors(n, 0);
    double maxError{fabs(X.at(0)) - explication(start)};

    for (int i = 0; i < n; ++i, start += h)
    {
        errors.at(i) = fabs(X.at(i)) - explication(start);
        if (maxError < errors.at(i))
            maxError = errors.at(i);
    }
    return maxError;
}

void saveFile(const string &fileName, const vector<double> &X, const double &h, const int &n, double xStep)
{
    try
    {
        fstream file(fileName.c_str(), ios::out);
        for (int i = 0; i < n; ++i, xStep += h)
            file << xStep << " " << explication(xStep) << " " << X.at(i) << endl;
        file.close();
    }
    catch (...)
    {
        cout << "Wystapil blad poczas zapisywania pliku " + fileName << endl;
        exit(1);
    }

}

double numberMethod(double h, int n, double x0, double x1)
{
    vector<double> B(n, 0);
    vector<double> D(n, 0);
    vector<double> L(n - 1, 0);
    vector<double> U(n - 1, 0);
    vector<double> X(n, 0);

    B.at(0) = -gamma;
    B.at(n - 1) = -theta;
    D.at(0) = beta_ - (alfa / h);
    D.at(n - 1) = psi + (fi / h);
    L.at(n - 2) = -fi / h;
    U.at(0) = alfa / h;

    double hPow{h * h};
    for (int i = 1; i < n - 1; ++i)
    {
        B.at(i) = (x0 + i * h - h) / 12.0 + ((10.0 / 12.0) * (x0 + i * h)) + (x0 + i * h + h) / 12.0;
        D.at(i) = ((-2.0 * p) / hPow) + (r * (10.0 / 12.0));
        L.at(i - 1) = U.at(i) = (p / hPow) + (r / 12.0);
    }
    thomasAlgorithm(B, D, L, U, X, n);
    saveFile("../numberMethod.txt", X, h, n, x0);

    return calculateError(X, h, n, x0);
}

double conventionalMethod(double h, int n, double x0, double x1)
{
    vector<double> B(n, 0);
    vector<double> D(n, 0);
    vector<double> L(n - 1, 0);
    vector<double> U(n - 1, 0);
    vector<double> X(n, 0);

    B.at(0) = -gamma;
    B.at(n - 1) = -theta;
    D.at(0) = beta_ - (alfa / h);
    D.at(n - 1) = psi + (fi / h);
    L.at(n - 2) = -fi / h;
    U.at(0) = alfa / h;

    double hPow{h * h};
    for (int i = 1; i < n - 1; ++i)
    {
        B.at(i) = i * h + x0;
        D.at(i) = r + (-2.0 * p) / hPow;
        L.at(i - 1) = U.at(i) = (p / hPow) - (q / (2.0 * h));
    }
    thomasAlgorithm(B, D, L, U, X, n);
    saveFile("../conventionalMethod.txt", X, h, n, x0);

    return calculateError(X, h, n, x0);
}