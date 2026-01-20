#include "Galerkin.h"
#include "solution.h"
#include "polynomials.h"
#include <iostream>
#include <fstream>
#include <random>
#include <functional>
#include <Eigen/Dense>

// 辅助函数，给文件名加后缀：filename.csv -> filename_suffix.csv
std::string make_filename(const std::string &base, const std::string &suffix)
{
    if (base.empty())
        return "";
    std::string name = base;
    std::string ext = "";
    size_t dot_pos = base.find_last_of('.');
    if (dot_pos != std::string::npos)
    {
        name = base.substr(0, dot_pos);
        ext = base.substr(dot_pos);
    }
    return name + "_" + suffix + ext;
}

// 辅助函数，保存误差 (N, Mean_Error, Var_Error)
void save_errors_csv(const std::string &filename,
                     const std::vector<int> &N,
                     const std::vector<double> &mean_errs,
                     const std::vector<double> &var_errs)
{
    if (filename.empty())
        return;
    std::string final_name = make_filename(filename, "errors");

    std::ofstream file(final_name);
    if (file.is_open())
    {
        file << "N,Mean_Error,Var_Error\n";
        for (size_t i = 0; i < N.size(); ++i)
        {
            file << N[i] << ","
                 << std::scientific << std::setprecision(6)
                 << mean_errs[i] << "," << var_errs[i] << "\n";
        }
    }
}

// 辅助函数，保存数值解和准确解的对比
void save_comparison_csv(const std::string &filename,
                         const std::string &grid_col_name,
                         const std::vector<double> &grid,
                         const std::vector<double> &exact,
                         const std::vector<int> &N,
                         const std::vector<std::vector<double>> &results)
{
    if (filename.empty())
        return;

    std::ofstream file(filename);
    if (file.is_open())
    {
        file << grid_col_name << ",Exact";
        for (int n : N)
            file << ",N=" << n;
        file << "\n";

        for (size_t i = 0; i < grid.size(); ++i)
        {
            file << std::fixed << std::setprecision(6) << grid[i] << "," << exact[i];
            for (size_t j = 0; j < N.size(); ++j)
            {
                if (i < results[j].size())
                    file << "," << results[j][i];
                else
                    file << ",NaN";
            }
            file << "\n";
        }
    }
}

// 辅助函数，保存N阶gpc展开的系数
void save_coefficients_csv(const std::string &base_filename,
                           int n,
                           const std::vector<double> &grid,
                           const std::vector<std::vector<double>> &coeffs)
{
    if (base_filename.empty())
        return;
    std::string final_name = make_filename(base_filename, "coeffs_N" + std::to_string(n));

    std::ofstream file(final_name);
    if (file.is_open())
    {
        file << "Grid"; // t或者x
        for (int k = 0; k <= n; ++k)
            file << ",y" << k;
        file << "\n";

        for (size_t i = 0; i < grid.size(); ++i)
        {
            file << std::fixed << std::setprecision(6) << grid[i];
            if (i < coeffs.size())
            {
                for (double val : coeffs[i])
                    file << "," << val;
            }
            file << "\n";
        }
    }
}

/// @brief 使用Galerkin + Runge Kutta方法求解ODE
/// @param N Hermite正交多项式的最高阶数
/// @param ode 要求解的ODE
/// @param t_grid ODE的时间间隔点
/// @param filename 保存的文件名
void solve1(const std::vector<int> &N, const solution::ODE &ode, const std::vector<double> &t_grid, const std::string &filename = "")
{
    std::vector<std::vector<double>> results;
    std::vector<double> mean;
    std::vector<double> var;
    int NN = t_grid.size() - 1;

    for (int n : N)
    {
        Galerkin::ODEResult s = Galerkin::solve_ode(n, ode, {0}, NN); // alpha = 0时的数值解
        std::vector<std::vector<double>> res = s.solutions;
        std::vector<double> alpha0;
        for (auto &r : res)
        {
            alpha0.push_back(r[0]);
        }
        results.push_back(alpha0);
        auto error = s.errors;
        mean.push_back(error.first);                                // 期望误差
        var.push_back(error.second);                                // 方差误差
        save_coefficients_csv(filename, n, t_grid, s.coefficients); // 保存每个N的gpc展开系数
    }

    std::vector<double> exact1 = solution::solution_ode_series(ode, NN); // 少一个初始点
    // std::cout << exact1.size() << std::endl;
    std::vector<double> exact;
    exact.reserve(NN + 1);
    exact.push_back(ode.solution_ode(t_grid[0]));
    exact.insert(exact.end(), exact1.begin(), exact1.end());
    save_comparison_csv(filename, "Time", t_grid, exact, N, results); // 保存alpha = 0的数值解与准确解
    save_errors_csv(filename, N, mean, var);                          // 保存误差
}

void solve2(const std::vector<int> &N, const solution::PDE &pde, const std::vector<double> &x_grid, const std::string &filename = "")
{
    std::vector<double> mean;
    std::vector<double> var;
    double T = pde.get_T();
    int NN = x_grid.size() - 1;
    for (int n : N)
    {
        auto s = Galerkin::solve_pde(n, pde, {1}, NN);
        std::vector<std::vector<double>> res = s.second;
        auto error = s.first;
        mean.push_back(error.first); // 期望误差
        var.push_back(error.second); // 方差误差
    }
    save_errors_csv(filename, N, mean, var);
}

void solve3(int N, const solution::PDE &pde,const std::vector<int>& NN)
{
    std::cout << std::left << std::setw(15) << "grids"
              << std::left << std::setw(20) << "mean_error"
              << std::left << std::setw(20) << "mean_order"
              << std::left << std::setw(20) << "var_error"
              << std::left << std::setw(20) << "var_order" << std::endl;
    double past1 = 0.0, past2 = 0.0;
    for (int i = 0; i < NN.size(); ++i)
    {
        int n = NN[i];
        std::vector<double> x_grid = poly::create_grid(0, 2 * poly::pi,n);
        int n_past = 0;
        if (i != 0)
        {
            n_past = NN[i - 1];
        }
        auto s = Galerkin::solve_pde(N, pde, {1}, n);
        std::pair<double, double> errors = s.first;
        double error1 = errors.first;
        double error2 = errors.second;
        double order1 = solution::order(past1, error1, n_past, n);
        double order2 = solution::order(past2, error2, n_past, n);
        std::cout << std::left << std::setw(15) << n
                  << std::left << std::setw(20) << error1
                  << std::left << std::setw(20) << order1
                  << std::left << std::setw(20) << error2
                  << std::left << std::setw(20) << order2 << std::endl;
        past1 = error1;
        past2 = error2;
    }
}

int main()
{
    Eigen::setNbThreads(1);
    std::vector<double> t_grid1 = poly::create_grid(0.0, 1.0, 4096);
    solution::ODE ode(0.0, 1.0, 0.0, 1.0, 1.0);
    std::cout << "=========== ODE Q1 ===========" << std::endl;
    solve1({0, 1, 2, 3, 4}, ode, t_grid1, "../HW5/ode_1.csv");
    std::cout << "=========== ODE Q2 ===========" << std::endl;
    solve1({1, 2, 3, 4, 5, 6, 7}, ode, t_grid1, "../HW5/ode_2.csv");

    std::vector<double> x_grid2 = poly::create_grid(0, 2 * poly::pi, 8192);
    solution::PDE pde(1.0, 1.0, 5.0, 2 * poly::pi, 0.8);
    std::cout << "=========== PDE ===========" << std::endl;
    solve2({0, 1, 2, 3, 4}, pde, x_grid2, "../HW5/pde_1.csv");

    std::cout << "=========== PDE  N=4===========" << std::endl;
    solve3(4,pde,{32,64,128,256});
}
