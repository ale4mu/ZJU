#include "Function.hpp"
#include "LinearFunction.hpp"
#include"Polynomial.hpp"
#include <iostream>
#include <complex>
#include <vector>


int main()
{
    std::cout << "test T=int" << std::endl;
    LinearFunction<int> f_int(2, 3);
    f_int.print(4);
    std::cout << "----------------------" << std::endl;

    std::cout << "test T=float" << std::endl;
    LinearFunction<float> f_float(2.5f, 3.5f);
    f_float.print(1.6f);
    std::cout << "----------------------" << std::endl;

    std::cout << "test T=double" << std::endl;
    LinearFunction<double> f_double(2.5, 3.5);
    f_double.print(1.6);
    std::cout << "----------------------" << std::endl;

    std::cout << "test T=complex" << std::endl;
    LinearFunction<std::complex<double>> f_complex({2.0, 2.0}, {3.0, 1.0});
    f_complex.print_complex({2.0, 2.0});
    std::cout << "----------------------" << std::endl;

    std::cout << "test T=vector(w is vector)" << std::endl;
    LinearFunction<std::vector<int>> f_vec1({1, 2,3}, {3, 4,5});
    f_vec1.print({1,1,1});
    std::cout << "----------------------" << std::endl;

    std::cout << "test T=vector(w is not vector)" << std::endl;
    LinearFunctionVec<int> f_vec2(2,{3,4,5});
    f_vec2.print({1,1,1});
    std::cout << "----------------------" << std::endl;

    std::cout << "test T=GMP" << std::endl;
    mpf_class w("1.234567890123456789", 100); // 100位精度
    mpf_class x("2.123566425613145611", 100); 
    mpf_class b("9.876543210987654321", 100);
    LinearFunction<mpf_class> f_GMP(w, b);
    f_GMP.print_GMP(x);// 10为小数
    std::cout << "----------------------" << std::endl;


    std::cout << "test Polynomial" << std::endl;
    std::vector<int> a = {5,2,3};
    Polynomial<int> p(a);
    p.print(2);

    return 0;
}