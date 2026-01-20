#ifndef LinearFunction_hpp
#define LinearFunction_hpp
#include "Function.hpp"
#include <complex>
#include <vector>
#include <stdexcept>
#include <iostream>
#include<algorithm>
#include <gmpxx.h>

// 默认支持标量和复数,主模板
template <typename T>
class LinearFunction : public Function<T>
{
private:
    T w;
    T b;

public:
    LinearFunction(T w, T b) : w(w), b(b) {}

    T operator()(const T &x) const override
    {
        return w * x + b;
    }

    ~LinearFunction() override = default;

    void print(const T &x) const
    {
        std::cout << w << "*" << x << "+" << b << "=" << (*this)(x) << std::endl;
    }

    void print_complex(const T &x) const //打印复数
    {
        T result = (*this)(x);
        std::cout << "(" << w.real() << "+" << w.imag() << "i)*(" << x.real() << "+" << x.imag() << "i)+(" << b.real() << "+" << b.imag() << "i)=" << result.real() << "+" << result.imag() << "i" << std::endl;
    }

    void print_GMP(const T& x)const //打印GMP
    {
        mpf_class result = (*this)(x);
        gmp_printf("%.15Ff * %.15Ff + %.15Ff = %.15Ff\n",w.get_mpf_t(),x.get_mpf_t(),b.get_mpf_t(),result.get_mpf_t());
    }
};

// 支持向量（输入和参数均为向量），特化模板
template <typename T>
class LinearFunction<std::vector<T>> : public Function<std::vector<T>>
{
private:
    std::vector<T> w; // w 和 b 也是向量
    std::vector<T> b;

public:
    LinearFunction(std::vector<T> w, std::vector<T> b) : w(w), b(b)
    {
        if (w.size() != b.size())
        {
            throw std::invalid_argument("w and b must have the same size!");
        }
    }

    std::vector<T> operator()(const std::vector<T> &x) const override
    {
        if (x.size() != w.size())
        {
            throw std::invalid_argument("x must have the same size as w and b!");
        }

        std::vector<T> result;
        for (size_t i = 0; i < x.size(); ++i)
        {
            result.push_back(w[i] * x[i] + b[i]);
        }
        return result;
    }

    ~LinearFunction() override = default;

    void print(const std::vector<T> &x) const //打印向量
    {
        std::cout << "(";
        for (auto it = w.begin(); it != w.end(); ++it)
        {
            std::cout << *it;
            if (next(it) != w.end())
            {
                std::cout << " , ";
            }
        }
        std::cout << ") * (";

        for (auto it = x.begin(); it != x.end(); ++it)
        {
            std::cout << *it;
            if (next(it) != x.end())
            {
                std::cout << " , ";
            }
        }
        std::cout << ") + (";

        for (auto it = b.begin(); it != b.end(); ++it)
        {
            std::cout << *it;
            if (next(it) != b.end())
            {
                std::cout << " , ";
            }
        }
        std::cout << ") = (";

        std::vector<T> result = (*this)(x);
        for (auto it = result.begin(); it != result.end(); ++it)
        {
            std::cout << *it;
            if (next(it) != result.end())
            {
                std::cout << " , ";
            }
        }
        std::cout << ")" << std::endl;
    }
};

// 支持标量参数对向量（w是标量，x 是向量）
template <typename T>
class LinearFunctionVec : public Function<std::vector<T>>
{
private:
    T w;
    std::vector<T> b;

public:
    LinearFunctionVec(T w, std::vector<T> b) : w(w), b(b) {}

    std::vector<T> operator()(const std::vector<T> &x) const override
    {
        std::vector<T> result;
        for (int i = 0; i < x.size(); i++)
        {
            result.push_back(w * x[i] + b[i]); 
        }
        return result;
    }

    ~LinearFunctionVec() override = default;

    void print(const std::vector<T> &x) const  //打印向量
    {
        std::cout << w << "*(";
         for (auto it = x.begin(); it != x.end(); ++it)
        {
            std::cout << *it;
            if (next(it) != x.end())
            {
                std::cout << " , ";
            }
        }
        std::cout << ") + (" ;

        for (int i = 0; i < b.size(); i++)
        {
            std::cout << b[i];
            if (i != b.size() - 1)
            {
                std::cout << " , ";
            }
        }
        std::cout << ") = (";

        std::vector<T> result = (*this)(x);
        for (auto it = result.begin(); it != result.end(); ++it)
        {
            std::cout << *it;
            if (next(it) != result.end())
            {
                std::cout << " , ";
            }
        }
        std::cout << ")" << std::endl;
    }
};

#endif