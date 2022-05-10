#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

double explication(double x)
{
    return (exp(2 - 2 * x) - 4 * exp(4 - 2 * x) * 4.0 * exp(2 + 2 * x) - x * x * exp(4)) / (4 - 4 * exp(4));
}

void
thomasAlgorithm(vector<double> X, vector<double> B, vector<double> L, vector<double> U, vector<double> D)
{
    const unsigned long long n{D.size()};
    for (int i = 0; i < n; ++i)
    {
        D.at(i) = D.at(i) - (L.at(i - 1) * U.at(i - 1) / D.at(i - 1));
        B.at(i) = B.at(i) - (L.at(i - 1) * D.at(i - 1) / B.at(i - 1));
    }
    X.at(n - 1) = 1 / D.at(n - 1) * B.at(n - 1);

    for (int i = n - 2; i >= 0; --i)
        X.at(i) = 1 / D.at(i) * (B.at(i) - U.at(i) * X.at(i + 1));

}

int main()
{
}
