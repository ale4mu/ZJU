#include "solution.h"
#include "polynomials.h"
#include <string>
#include <fstream>
#include <iostream>
#include <omp.h>
using ExactFunc = std::function<std::pair<double, double>(double)>;

std::pair<double, double> process(
    int N,
    double step, // dt 或 dx
    const std::vector<double> &mean_y,
    const std::vector<double> &mean_y2,
    const std::string &filename,
    ExactFunc exact_func,
    int start_index = 0 // ODE是1
)
{
    std::ofstream outfile(filename);
    outfile << "Coord,Mean_Num,Mean_Exact,Var_Num,Var_Exact\n";
    outfile << std::fixed << std::setprecision(8);

    int _size = mean_y.size() - start_index;
    std::vector<double> numer_means(_size), exact_means(_size);
    std::vector<double> numer_vars(_size), exact_vars(_size);

    for (int i = 0; i < _size; ++i)
    {
        int idx = i + start_index;
        double coord = idx * step;

        // 数值解
        double mean = mean_y[idx];
        double var = mean_y2[idx] - mean * mean;

        // 期望解
        auto [exact_mean, exact_var] = exact_func(coord);
        numer_means[i] = mean;
        exact_means[i] = exact_mean;
        numer_vars[i] = var;
        exact_vars[i] = exact_var;

        // 写入文件
        outfile << coord << "," << mean << "," << exact_mean << ","
                << var << "," << exact_var << "\n";
    }

    outfile.close();
    double mean_error = solution::L_Chebyshev(numer_means, exact_means);
    double var_error = solution::L_Chebyshev(numer_vars, exact_vars);
    return {mean_error, var_error};
}

/// @brief 用随机配点法 + Runge-Kutta法求解ODE
/// @tparam M 配置点个数
/// @param N 时间离散度（离散时间步数）
/// @return （期望误差，方差误差）
template <int M>
std::pair<double, double> solve_ode_template(int N)
{
    poly::Hermite<M> hermite(poly::Method::USE_LIBRARY, poly::HermiteType::PROBABILITY);// Hermite多项式对象（概率型）
    const auto &nodes = hermite.get_points(); // 高斯积分节点
    const auto &weights = hermite.get_weights(); // 高斯积分权重
    int num_points = nodes.size(); // 节点个数

    std::vector<double> mean_y(N + 1, 0.0); // 存储数值解每个时间点的期望解
    std::vector<double> mean_y2(N + 1, 0.0); // 存储数值解每个时间点的期望解的平方

    solution::ODE base_ode(0.0, 1.0, 0.0, 1.0, 1.0); 
    double dt = base_ode.get_dt(N); // 获取ODE的时间步长

#pragma omp parallel
    {
        std::vector<double> local_mean_y(N + 1, 0.0); // 为了并行化，存储单个线程的数值解的期望
        std::vector<double> local_mean_y2(N + 1, 0.0);

#pragma omp for
        for (int j = 0; j < num_points; ++j)
        {
            double w_j = weights[j]; // 高斯积分权重
            solution::ODE ode(nodes[j], 1.0, 0.0, 1.0, 1.0); // 用当前高斯积分节点(z_j)作为ODE的alpha，解确定的ODE
            std::vector<double> u_res = solution::runge_kutta(ode, N);// 返回当前z_j对应的数值解（不包括初值）

            double val = ode.get_u0(); // 加上初值
            mean_y[0] += val * w_j;
            mean_y2[0] += (val * val) * w_j;

            for (int t = 1; t <= N; ++t)
            {
                val = u_res[t - 1];
                local_mean_y[t] += val * w_j; // 最终的数值解的期望是当前时间步 每个z_j的数值解 * 高斯积分权重 的求和
                local_mean_y2[t] += (val * val) * w_j;
            }
        }

// 为了并行化，每个线程单独计算，最后再求和
#pragma omp critical
        {
            for (int i = 0; i <= N; ++i)
            {
                mean_y[i] += local_mean_y[i];
                mean_y2[i] += local_mean_y2[i];
            }
        }
    }

    // 期望解
    auto ode_exact = [](double t) -> std::pair<double, double>
    {
        double ex_mean = std::exp(t * t / 2.0);
        double ex_var = std::exp(2.0 * t * t) - std::exp(t * t);
        return {ex_mean, ex_var};
    };

    // 传入process函数对数值解的期望(mean_y)和期望解(ex_mean)计算绝对值误差
    return process(N, dt, mean_y, mean_y2,
                   "../HW6/ode_" + std::to_string(N) + "_" + std::to_string(M) + ".csv",
                   ode_exact, 1);
}

std::pair<double, double> solve_ode(int N, int M = 10)
{
    switch (M)
    {
    case 1:
        return solve_ode_template<1>(N);
    case 2:
        return solve_ode_template<2>(N);
    case 3:
        return solve_ode_template<3>(N);
    case 4:
        return solve_ode_template<4>(N);
    case 5:
        return solve_ode_template<5>(N);
    case 6:
        return solve_ode_template<6>(N);
    case 7:
        return solve_ode_template<5>(N);
    case 8:
        return solve_ode_template<7>(N);
    case 10:
        return solve_ode_template<10>(N);
    default:
        std::cerr << "Unsupported M=" << M << std::endl;
        return {0.0, 0.0};
    }
}


/// @brief 用随机配点法 + scheme3求解PDE
/// @tparam M 配置点个数
/// @param N 空间离散度（离散空间步数）
/// @return （期望误差，方差误差）
template <int M>
std::pair<double, double> solve_pde_template(int N)
{
    poly::Legendre<M> legendre(poly::Method::USE_LIBRARY); // Legendre多项式对象
    const auto &nodes = legendre.get_points(); // 高斯积分节点
    const auto &weights = legendre.get_weights(); // 高斯积分权重
    int num_points = nodes.size(); // 节点个数

    std::vector<double> mean_y(N, 0.0);
    std::vector<double> mean_y2(N, 0.0);

    solution::PDE base_pde(1.0, 1.0, 5.0, 2 * poly::pi, 0.8); 
    const double dx = base_pde.get_dx(N); // 离散空间步长
    const double T = base_pde.get_T(); // 周期
#pragma omp parallel
    {
        std::vector<double> local_mean_y(N, 0.0);
        std::vector<double> local_mean_y2(N, 0.0);

#pragma omp for
        for (int j = 0; j < num_points; ++j)
        {
            double w_j = weights[j] / 2.0; // 高斯积分权重, 除以2是因为均匀分布[-1, 1]的概率密度函数为1/2
            // 用当前高斯积分节点(z_j)计算PDE的a(即1+0.1z)，解确定的PDE
            solution::PDE pde(1 + 0.1 * nodes[j], 1.0, 5.0, 2 * poly::pi, 0.8); 
            std::vector<double> u_res = solution::scheme3(pde, N); // 返回当前z_j(a)对应的数值解

            for (int x = 0; x < N; ++x)
            {
                double val = u_res[x];
                // 只计算最后一个时间步的物理空间误差
                local_mean_y[x] += val * w_j; // 最终的数值解的期望是当前x 每个z_j的数值解 * 高斯积分权重 的求和
                local_mean_y2[x] += (val * val) * w_j;
            }
        }

// 并行化，每个线程单独计算，最后再求和
#pragma omp critical
        {
            for (int i = 0; i < N; ++i)
            {
                mean_y[i] += local_mean_y[i];
                mean_y2[i] += local_mean_y2[i];
            }
        }
    }

    // 期望解
    
    auto pde_exact = [T](double x) -> std::pair<double, double>
    {
        double exact_mean = std::sin(x + T) * std::sin(0.1 * T) / (0.1 * T);
        double exact_var = 0.5 * (1 - std::cos(2 * (x + T)) * std::sin(0.2 * T) / (0.2 * T)) - exact_mean * exact_mean;
        return {exact_mean, exact_var};
    };

    // 传入process函数对数值解的期望(mean_y)和期望解(exact_mean)计算绝对值误差
    return process(N, dx, mean_y, mean_y2,
                   "../HW6/pde_" + std::to_string(N) + "_" + std::to_string(M) + ".csv",
                   pde_exact, 0);
}

std::pair<double, double> solve_pde(int N, int M = 10)
{
    switch (M)
    {
    case 1:
        return solve_pde_template<1>(N);
    case 2:
        return solve_pde_template<2>(N);
    case 3:
        return solve_pde_template<3>(N);
    case 4:
        return solve_pde_template<4>(N);
    case 5:
        return solve_pde_template<5>(N);
    case 6:
        return solve_pde_template<6>(N);
    case 7:
        return solve_pde_template<7>(N);
    case 8:
        return solve_pde_template<8>(N);
    case 10:
        return solve_pde_template<10>(N);
    default:
        std::cerr << "Unsupported M=" << M << std::endl;
        return {0.0, 0.0};
    }
}

void test1(const std::vector<int> &Ns, std::function<std::pair<double, double>(int, int)> solver)
{
    std::cout << std::left << std::setw(15) << "N"
              << std::left << std::setw(20) << "mean_error"
              << std::left << std::setw(20) << "mean_order"
              << std::left << std::setw(20) << "var_error"
              << std::left << std::setw(20) << "var_order" << std::endl;

    double past1 = 0.0, past2 = 0.0;

    for (int i = 0; i < Ns.size(); ++i)
    {
        int n = Ns[i];
        int n_past = (i != 0) ? Ns[i - 1] : 0;
        auto [error1, error2] = solver(n, 10);
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

void test2(int N, const std::vector<int> &Ms, std::function<std::pair<double, double>(int, int)> solver, const std::string &filename)
{
    std::ofstream outfile(filename);
    outfile << "M,Mean_Error,Var_Error\n";
    outfile << std::scientific << std::setprecision(16);
    std::cout << std::left << std::setw(15) << "多项式阶数"
              << std::left << std::setw(20) << "mean_error"
              << std::left << std::setw(20) << "var_error" << std::endl;
    for (int M : Ms)
    {
        auto [error1, error2] = solver(N, M);
        std::cout << std::left << std::setw(15) << M
                  << std::left << std::setw(20) << error1
                  << std::left << std::setw(20) << error2 << std::endl;
        outfile << M << "," << error1 << "," << error2 << "\n";
    }
    outfile.close();
}

int main()
{
    test1({8, 16, 32, 64}, solve_ode);
    test1({32, 64, 128}, solve_pde);
    test2(128, {1, 2, 3, 4, 5, 6}, solve_ode, "../HW6/ode_test.csv");
    test2(4096, {1, 2, 3, 4, 5, 6, 7}, solve_pde, "../HW6/pde_test.csv");
}