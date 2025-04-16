#include "fraction.h"

int main()
{
    Fraction f1(1, 2);
    Fraction f2(3, 4);

    cout << "f1: " << f1 << endl;
    cout << "f2: " << f2 << endl;

    cout << "f1 +  = 4/5 : " << f1 + f2 << endl;
    cout << "f1 - f2 = 1/4 : " << f2 - f1 << endl;
    cout << "f1 * f2 = 3/8 : " << f1 * f2 << endl;
    cout << "f1 / f2 = 2/3 : " << f1 / f2 << endl;

    cout << "test f1 == f2: " << (f1 == f2) << endl;
    cout << "test f1 <= f2: " << (f1 <= f2) << endl;
    cout << "test f1 < f2: " << (f1 < f2) << endl;
    cout << "test f1 >= f2: " << (f1 >= f2) << endl;
    cout << "test f1 > f2: " << (f1 > f2) << endl;
    cout << "test f1 != f2: " << (f1 != f2) << endl;

    double num = 0.85;
    Fraction f3 = Fraction().to_fraction(num);
    cout << "0.85 to fraction: " << f3 << endl;

    Fraction f4;
    string input = "5/8"; // 模拟用户输入
    istringstream iss(input);
    cin.rdbuf(iss.rdbuf()); // 重定向，将 cin 的输入缓冲区替换为 iss 的缓冲区
    cin >> f4;
    cout << "test cin 5/8: " << f4 << endl;

    return 0;
}