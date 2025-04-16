#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>
using namespace std;

class Fraction
{
private:
    int numerator;   // 分子
    int denominator; // 分母

    int gcd(int a, int b) const // 求最大公约数
    {
        return b == 0 ? a : gcd(b, a % b);
    }

    void reduce() // 约分
    {
        int g = gcd(abs(numerator), abs(denominator));
        numerator /= g;
        denominator /= g;
    }

public:
    Fraction() {} // 默认构造函数

    Fraction(int n, int d) : numerator(n), denominator(d)
    {
        if (d == 0)
            throw invalid_argument("Denominator cannot be 0!");
        reduce(); // 约分
    }

    Fraction(const Fraction &f) : Fraction(f.numerator, f.denominator) {} // 拷贝构造函数

    // 重载算术运算符
    Fraction operator+(const Fraction &f1) // 加法
    {
        int new_denominator = this->denominator * f1.denominator;
        int new_numerator = this->numerator * f1.denominator + f1.numerator * this->denominator;
        return Fraction(new_numerator, new_denominator);
    }

    Fraction operator-(const Fraction &f1) // 减法
    {
        int new_denominator = this->denominator * f1.denominator;
        int new_numerator = this->numerator * f1.denominator - f1.numerator * this->denominator;
        return Fraction(new_numerator, new_denominator);
    }

    Fraction operator*(const Fraction &f1) // 乘法
    {
        int new_denominator = this->denominator * f1.denominator;
        int new_numerator = this->numerator * f1.numerator;
        return Fraction(new_numerator, new_denominator);
    }

    Fraction operator/(const Fraction &f1) // 除法
    {
        if (f1.numerator == 0)
            throw invalid_argument("Cannot / 0!");
        int new_denominator = this->denominator * f1.numerator;
        int new_numerator = this->numerator * f1.denominator;
        return Fraction(new_numerator, new_denominator);
    }

    // 重载关系运算符
    bool operator==(const Fraction &f)
    {
        return this->denominator * f.numerator == this->numerator * f.denominator;
    }

    bool operator<=(const Fraction &f)
    {
        // return this->denominator * f.numerator >= this->numerator * f.denominator;
        Fraction temp = *this - f;
        return temp.numerator>= 0;
    }

    bool operator<(const Fraction &f)
    {
        // return this->denominator * f.numerator > this->numerator * f.denominator;
        Fraction temp = *this - f;
        return temp.numerator> 0;
    }

    bool operator>=(const Fraction &f)
    {
        // return this->denominator * f.numerator <= this->numerator * f.denominator;
        Fraction temp = *this - f;
        return temp.numerator<= 0;
    }

    bool operator>(const Fraction &f)
    {
        // return this->denominator * f.numerator < this->numerator * f.denominator;
        Fraction temp = *this - f;
        return temp.numerator < 0;
    }

    bool operator!=(const Fraction &f)
    {
        return !(*this == f);
    }

    explicit operator double() const // typecast to double
    {
        return static_cast<double>(numerator) / denominator;
    }

    explicit operator string() const // typecast to string
    {
        stringstream ss;
        ss << numerator << "/" << denominator;
        return ss.str();
    }

    friend ostream &operator<<(ostream &out, const Fraction &f) // inserter
    {
        out << f.numerator << "/" << f.denominator;
        return out;
    }

    friend istream &operator>>(istream &in, Fraction &f) // extractor
    {
        char c;
        in >> f.numerator >> c >> f.denominator;
        return in;
    }

    Fraction to_fraction(double a) // 有限小数化分数
    {
        int denominator = 1;
        while (a - (int)a > 0)
        {
            a *= 10;
            denominator *= 10;
        }
        return Fraction((int)a, denominator);
    }
};
