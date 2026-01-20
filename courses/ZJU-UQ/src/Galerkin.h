#ifndef GALERKIN_H
#define GALERKIN_H

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <Eigen/Dense>
#include "solution.h"
#include "polynomials.h"
#include <omp.h>
namespace Galerkin
{
    // 计算阶乘的辅助函数
    long long factorial(int n)
    {
        if (n < 0)
            return 0;
        long long res = 1.0;
        for (int i = 1; i <= n; ++i)
        {
            res *= i;
        }
        return res;
    }

    // 由于alpha服从标准正态分布，a0 = 0, a1 = 1
    // 这里的A实际上是笔记中的A的转置
    Eigen::MatrixXd A_build_ode(int N)
    {
        int size = N + 1;
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(size, size);
        if (N == 0)
            return A;
        for (int j = 0; j < size; ++j)
        {
            long long jj = factorial(j); // j!
            for (int k = 0; k < size; ++k)
            {
                long long gamma = factorial(k); // γ_k = k!
                if (gamma == 0)
                    continue;
                double e_1jk = 0.0;
                int sum_ijk = 1 + j + k; // i=1
                if (sum_ijk % 2 == 0)
                {
                    int s = sum_ijk / 2; // s = (1+j+k)/2
                    if (s >= 1 && s >= j && s >= k)
                    {
                        // 用笔记中的公式直接计算
                        e_1jk = static_cast<double>(jj * factorial(k)) / (factorial(s - 1) * factorial(s - j) * factorial(s - k));
                    }
                }
                A(j, k) = (1.0 / gamma) * e_1jk;
            }
        }

        return A;
    }

    // 这个函数用于计算每一步的 Runge-Kutta 迭代
    Eigen::VectorXd rk3_step(const Eigen::VectorXd &y, const Eigen::MatrixXd &A, double dt)
    {
        Eigen::VectorXd y1 = y + dt * (A * y);
        Eigen::VectorXd y2 = 0.75 * y + 0.25 * y1 + 0.25 * dt * (A * y1);
        Eigen::VectorXd y_next = (1.0 / 3.0) * y + (2.0 / 3.0) * y2 + (2.0 / 3.0) * dt * (A * y2);
        return y_next;
    }

    // 存储结果
    struct ODEResult
    {
        std::pair<double, double> errors;              // (均值误差, 方差误差)
        std::vector<std::vector<double>> solutions;    // gpc解
        std::vector<std::vector<double>> coefficients; // 投影系数V_hat
    };

    /// @param N 正交多项式最高阶数
    /// @param alpha_grid 正交多项式的网格点（实际上是alpha作为随机变量的不同取值）
    /// @param NN 划分ODE的时间间隔
    ODEResult solve_ode(int N, const solution::ODE ode, const std::vector<double> &alpha_grid, int NN)
    {
        int size = N + 1;
        Eigen::MatrixXd A_T = A_build_ode(N); // 矩阵A
        Eigen::MatrixXd A = (-1.0) *A_T.transpose(); // 这里与笔记中不同，应该要取负号
        // Hermite多项式（概率型），有200个积分节点
        poly::Hermite<200> hermite_poly(poly::Method::USE_LIBRARY, poly::HermiteType::PROBABILITY); 

        // 初始条件 V(0) = [1(beta), 0, ..., 0]^T
        Eigen::VectorXd V = Eigen::VectorXd::Zero(size);
        V(0) = ode.get_u0(); // 初值
        std::vector<std::vector<double>> solutions; // gpc解
        std::vector<std::vector<double>> coefficients; // gpc系数
        double t0 = ode.get_t0();
        double t1 = ode.get_t1();
        double dt = (t1 - t0) / NN; // 时间离散度
        std::vector<double> Mean;
        std::vector<double> Var;
        std::vector<double> exact_mean;
        std::vector<double> exact_var;
        // 记录初始条件
        {
            std::vector<double> V_hat(V.data(), V.data() + size); 
            Mean.push_back(V_hat[0]); // 数值解的期望
            Var.push_back(0.0); // 数值解的方差
            exact_mean.push_back(1.0); // 精确解的期望
            exact_var.push_back(0.0); // 精确解的方差
            solutions.push_back(poly::evaluate(V_hat, hermite_poly, alpha_grid));
            coefficients.push_back(V_hat);
        }
        int count = 0;

        // runge-kutta的迭代过程
        while (t0 < t1)
        {
            if (t0 + dt > t1)
            {
                dt = t1 - t0;
            }
            V = rk3_step(V, A, dt); // 用runge-kutta解出当前时间步的V
            std::vector<double> V_hat(V.data(), V.data() + size); // 类型转换
            Mean.push_back(V_hat[0]); // 记录均值（均值是gpc投影系数的第一个元素）
            double sum = 0.0;

            // 计算方差，除了第一个元素外，...求和
            for (int i = 1; i < size; ++i)
            {
                sum += V_hat[i] * V_hat[i] * factorial(i);
            }
            Var.push_back(sum);
            double mu_exact = std::exp(t0 * t0 / 2.0);                      // 精确解的数学期望是exp(t^2/2)
            double var_exact = std::exp(2.0 * t0 * t0) - std::exp(t0 * t0); // 精确解的方差是exp(2t^2) - exp(t^2)
            exact_mean.push_back(mu_exact);
            exact_var.push_back(var_exact);
            // 这里容易出错，求和与正交多项式的投影有一点不同
            // V的自变量是t，而正交多项式的自变量是随机变量alpha
            solutions.push_back(poly::evaluate(V_hat, hermite_poly, alpha_grid)); // 记录当前时间步的gpc解
            coefficients.push_back(V_hat); // 记录当前时间步的gpc系数
            t0 += dt;
            ++count;
        }
        // 计算绝对值误差
        double error_mean = solution::L_Chebyshev(Mean, exact_mean);
        double error_var = solution::L_Chebyshev(Var, exact_var);
        std::cout << "N = : " << N << std::endl;
        std::cout << "error_mean: " << error_mean << std::endl;
        std::cout << "error_var: " << error_var << std::endl;

        ODEResult result;
        result.errors = std::make_pair(error_mean, error_var);
        result.solutions = solutions;
        result.coefficients = coefficients;
        return result;
    }

    // 构建矩阵A
    Eigen::MatrixXd A_build_pde(int N)
    {
        int size = N + 1;
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(size, size);

        for (int i = 0; i < size; ++i)
        {

            A(i, i) = 1.0; // 注意归一化
            if (i - 1 >= 0)
            {
                A(i, i - 1) += 0.1 * (double)i / (2.0 * i - 1.0);
            }
            if (i + 1 < size)
            {
                A(i, i + 1) += 0.1 * (double)(i + 1) / (2.0 * i + 3.0);
            }
        }
        return A;
    }

    // scheme3
    std::vector<Eigen::VectorXd> scheme3_step(const std::vector<Eigen::VectorXd> &U_old, const Eigen::MatrixXd &A,
                                              double ratio)
    {
        int size = U_old.size(); // 空间网格数
        std::vector<Eigen::VectorXd> U_new(size);
#pragma omp parallel for schedule(static)
        for (int j = 0; j < size; ++j)
        {
            int j_next = (j + 1) % size;
            Eigen::VectorXd diff = U_old[j_next] - U_old[j];
            // U{j+1} = U{j} + (dt/dx) * A * (U{j+1} - U{j})
            U_new[j] = U_old[j] + ratio * A * diff;
        }
        return U_new;
    }

    /// @brief 用Galerkin + scheme3 求解PDE
    /// @param N 正交多项式的最高阶数
    /// @param pde 要求解的PDE
    /// @param alpha_grid 正交多项式的网格点（实际上是alpha作为随机变量的不同取值）
    /// @param NN 划分PDE的空间坐标数
    /// @return （（期望误差，方差误差），解）
    std::pair<std::pair<double, double>, std::vector<std::vector<double>>> solve_pde(int N, const solution::PDE pde, const std::vector<double> &alpha_grid, int NN)
    {
        int size = N + 1;
        Eigen::MatrixXd A = A_build_pde(N);
        poly::Legendre<200> legendre_poly(poly::Method::USE_LIBRARY); // legendre多项式
        const double a = pde.get_a();
        const double k = pde.get_k();
        const double T = pde.get_T();
        const double L = pde.get_L();
        const double cfl = pde.get_cfl();
        const double dx = L / NN;
        double dt = cfl * dx;
        std::vector<Eigen::VectorXd> U(NN, Eigen::VectorXd::Zero(size)); // 即笔记中的V
        std::vector<std::vector<double>> solutions; // 记录gpc解
        solutions.resize(NN);

        // 初始条件u(x,0) = sinx
        for (int j = 0; j < NN; ++j)
        {
            U[j](0) = std::sin(j * dx);
        }

        double current_t = 0.0;
        double ratio = cfl; // dt/dx
        int count = 0;
        // 先走完每个时间步
        while (current_t < T)
        {
            // 走完最后一个时间步
            if (current_t + dt > T)
            {
                dt = T - current_t;
                ratio = dt / dx;
            }
            U = scheme3_step(U, A, ratio);
            current_t += dt;
            ++count;
        }

        std::vector<double> Mean(NN); // 存储数值解的期望
        std::vector<double> Var(NN);
        std::vector<double> exact_mean(NN); // 存储精确解的期望
        std::vector<double> exact_var(NN);

// 只在最后一个时间步上计算物理空间的误差
#pragma omp parallel for schedule(static)
        for (int j = 0; j < NN; ++j)
        {
            //     if (j == 0)
            // {
            //     int num_threads = omp_get_num_threads();
            //     std::cout << "Number of threads: " << num_threads << std::endl;
            // }
            double current_x = j * dx; // 当前x
            std::vector<double> U_current(U[j].data(), U[j].data() + size);
            Mean[j] = U_current[0]; // 记录均值（第一个元素）
            double sum = 0.0;
            // 计算方差（从第二个元素开始）
            for (int i = 1; i < size; ++i)
            {
                sum += U_current[i] * U_current[i] * legendre_poly.norm(i) * 0.5; // 这里legendre多项式要*0.5
            }
            Var[j] = sum;
            // 准确解的期望和方差
            double mu_exact = std::sin(current_x + T) * std::sin(0.1 * T) / (0.1 * T); 
            double var_exact = 0.5 * (1 - std::cos(2 * (current_x + T)) * std::sin(0.2 * T) / (0.2 * T)) - mu_exact * mu_exact;
            exact_mean[j] = mu_exact;
            exact_var[j] = var_exact;
            // 实际上最后一题只需要误差，因此只需要gpc系数
            solutions[j] = poly::evaluate(U_current, legendre_poly, alpha_grid); // 计算gpc解
        }

        double error_mean = solution::L_Chebyshev(Mean, exact_mean);
        double error_var = solution::L_Chebyshev(Var, exact_var);
        // std::cout << "N = : " << N << std::endl;
        // std::cout << "error_mean: " << error_mean << std::endl;
        // std::cout << "error_var: " << error_var << std::endl;
        std::pair<double, double> error = std::make_pair(error_mean, error_var);
        return std::make_pair(error, solutions);
    }
}

#endif // GALERKIN_H