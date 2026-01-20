#include "polynomials.h"
#include <iostream>

using namespace poly;

/// @brief 对单个函数做正交多项式投影
/// @param name 函数名
/// @param f 函数
/// @param poly 正交多项式
/// @param N_values 正交多项式的最高阶数
/// @param x_grid 网格
/// @return 不同N值的误差
std::vector<double> test(
    const std::string &name,
    const std::function<double(double)> &f,
    const Polynomial &poly,
    const std::vector<int> &N_values,
    const std::vector<double> &x_grid)
{
    std::string filename = "../HW3/" + poly.name() + "_" + name + ".csv";
    std::ofstream outfile(filename);

    std::cout << poly.name() << ": " << name << std::endl;
    std::cout << "N\tL2_Error" << std::endl;

    std::map<int, std::vector<double>> results;
    std::vector<double> errors;
    for (int N : N_values)
    {
        auto coeffs = projection(f, poly, N);        // 计算投影系数
        results[N] = evaluate(coeffs, poly, x_grid); // 在网格上计算
        double error = l2_error(f, coeffs, poly);    // l2误差
        errors.push_back(error);
        std::cout << N << "\t" << std::scientific << std::setprecision(4) << error << std::endl;
    }

    outfile << "x,original";
    for (int N : N_values)
    {
        outfile << ",N=" << N;
    }
    outfile << "\n";

    for (size_t i = 0; i < x_grid.size(); ++i)
    {
        outfile << x_grid[i] << "," << f(x_grid[i]);
        for (int N : N_values)
        {
            outfile << "," << results[N][i];
        }
        outfile << "\n";
    }
    outfile.close();
    return errors;
}

void run_test(
    const std::string &filename,
    const Polynomial &poly,
    const std::vector<double> &x_grid,
    const std::map<std::string, std::function<double(double)>> &functions,
    const std::vector<int> &N_values)
{
    std::ofstream outfile(filename);
    outfile << "Function,N,Error\n";

    for (const auto &pair : functions)
    {
        const std::string &func_name = pair.first;
        const std::function<double(double)> &func = pair.second;

        std::vector<double> errors = test(func_name, func, poly, N_values, x_grid);

        for (size_t i = 0; i < N_values.size(); ++i)
        {
            outfile << func_name << "," << N_values[i] << "," << errors[i] << "\n";
        }
    }

    outfile.close();
}

int main()
{
    auto f1 = [](double x)
    { return std::pow(std::abs(std::sin(pi * x)), 3); };
    auto f2 = [](double x)
    { return std::abs(x); };
    auto f3 = [](double x)
    { return std::cos(pi * x); };
    auto f4 = [](double x)
    { return (x > 0.0) ? 1.0 : ((x < 0.0) ? -1.0 : 0.0); };

    std::map<std::string, std::function<double(double)>> functions = {
        {"f1", f1},
        {"f2", f2},
        {"f3", f3},
        {"f4", f4}};

    std::vector<int> N_values = {2, 4, 8, 16, 32};

    std::vector<double> legendre_res = create_grid(-1.0, 1.0, 800); // [-1, 1]
    std::vector<double> hermite_res = create_grid(-2.0, 2.0, 800);  // [-2,2]

    // 实例化多项式对象
    Legendre<200> legendre_poly;
    Hermite<200> hermite_poly;

    run_test("../HW3/Legendre_errors.csv", legendre_poly, legendre_res, functions, N_values);
    run_test("../HW3/Hermite_errors.csv", hermite_poly, hermite_res, functions, N_values);

    return 0;
}