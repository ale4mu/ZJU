#ifndef Polynomial_HPP
#define Polynomial_HPP
using namespace std;
#include "Function.hpp"
#include <vector>
#include <cmath>

template <typename T>
class Polynomial : public Function<T>
{

private:
    vector<T> a; //多项式系数

public:
    Polynomial(vector<T> _a) : a(_a) {}

    T operator()(const T &x) const override 
    {
        T result = 0;
        T power_x = 1; 
        for (const auto &val : a)
        {
            result += val * power_x;
            power_x *= x;
        }
        return result;
    }

    ~Polynomial() override = default;

    void print(const T &x) //打印多项式和结果
    {
        cout << "x=" << x << ", f(x) = ";
        for (int i = 0; i < a.size(); i++)
        {
            cout << a[i] << "*x^" << i;
            if (i != a.size() - 1)
            {
                cout << " + ";
            }
        }
        cout << " = " << (*this)(x) << endl;
    }
};

#endif