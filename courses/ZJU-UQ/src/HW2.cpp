#include "solution.h"
#include <random>
#include <omp.h>
// beta = 1时，数学期望为exp(t^2/2)
std::vector<double> expectation_ode(const solution::ODE &ode, int N)
{
    double t0 = ode.get_t0();
    double t1 = ode.get_t1();
    std::vector<double> T(N);
    double dt = (t1 - t0) / N;
    for (int i = 0; i < N; i++)
    {
        double tt = (t0 + (i + 1) * dt);
        T[i] = exp(tt * tt / 2.0);
    }
    return T;
}



/// @brief 用蒙特卡洛法求解ODE
/// @param M M个独立同分布的随机数
/// @param N 迭代次数
/// @param func 求解ODE的函数
std::vector<double> ODE_monte_carlo(int M, int N, std::function<std::vector<double>(const solution::ODE &, int)> func,
                                    std::mt19937 &generator)
{
    using namespace solution;
    std::normal_distribution<double> normal(0.0, 1.0); // alpha服从标准正态分布
    const double beta = 1.0;

    std::vector<double> res(N, 0); // 存储解的和
    for (int i = 0; i < M; ++i)
    {
        double alpha = normal(generator); // 随机生成alpha
        ODE ode(alpha, beta, 0.0, 2.0, 1.0); // 生成ODE
        std::vector<double> _res = func(ode, N);// 用显式欧拉法/Runge-Kutta方法求解ODE
        res += _res;
    }
    // res /= M; // 在调用处再取平均
    return res;
}

/// @brief 单组实验，用增量序列增加随机数的个数
/// @param func 求解ODE/PDE的方法
/// @param M 随机数个数
/// @param N 迭代次数
/// @param generator 随机数生成器
/// @return 每个M处的L2误差和L无穷误差
template <typename MonteCarlo, typename Expectation>
std::pair<std::vector<double>, std::vector<double>> solve(MonteCarlo monte_carlo, Expectation expectation,
                                                              const std::vector<int> &M, const int N, std::mt19937 &generator)
{
    using namespace solution;
    std::vector<double> res(N, 0);                  // 用于存储增量序列
    std::vector<double> l2_errors(M.size(), 0.0);   // 用于存储L2误差
    std::vector<double> linf_errors(M.size(), 0.0); // 用于存储L无穷误差
    for (int i = 0; i < M.size(); ++i)
    {
        int m = M[i];
        int m_past = 0;
        if (i != 0)
        {
            m_past = M[i - 1];
        }
        int m_need = m - m_past;
        if (m_need > 0)
        {
            auto vec = monte_carlo(m_need, N, generator);
            res += vec;
        }
        std::vector<double> avg = res;
        avg /= static_cast<double>(m);
        std::vector<double> T = expectation(N);
        //save("../HW2/scheme3_k1_T5_N160.csv",avg,T,N);
        std::pair<double, double> errors = error(avg, T);
        double error1 = errors.first;
        double error2 = errors.second;
        l2_errors[i] = error1;   // L2误差
        linf_errors[i] = error2; // L无穷误差
    }
    return std::make_pair(l2_errors, linf_errors);
}

// 重复多组实验，计算误差阶
template <typename MonteCarlo, typename Expectation>
void repeat(const std::string &name, const int K, const std::vector<int> &M, const int N,
            MonteCarlo monte_carlo, Expectation expectation)
{
    using namespace solution;
    std::cout << "---" << name << "---" << std::endl;
    std::cout << std::left << std::setw(15) << "N"
              << std::left << std::setw(15) << "M"
              << std::left << std::setw(20) << "L2_error"
              << std::left << std::setw(20) << "L2_order"
              << std::left << std::setw(20) << "L_inf_error"
              << std::left << std::setw(20) << "L_inf_order" << std::endl;
    std::vector<double> l2_errors(M.size(), 0.0);
    std::vector<double> linf_errors(M.size(), 0.0);
    omp_set_num_threads(omp_get_max_threads());

// K次重复实验
#pragma omp parallel for schedule(dynamic)
    for (int k = 0; k < K; ++k)
    {
        // 打印线程数
        // if (k == 0)
        // {
        //     int num_threads = omp_get_num_threads();
        //     std::cout << "Number of threads: " << num_threads << std::endl;
        // }
        std::random_device rd;
        std::mt19937 generator(rd());
        auto [res1, res2] = solve(monte_carlo, expectation, M, N, generator); // 蒙特卡洛法解ODE，返回误差
#pragma omp critical
        {
            for (size_t i = 0; i < M.size(); ++i)
            {
                l2_errors[i] += res1[i];
                linf_errors[i] += res2[i];
            }
        }
    }
    l2_errors /= K; // 对误差取平均
    linf_errors /= K;
    double past1 = 0.0, past2 = 0.0;
    for (int i = 0; i < M.size(); ++i)
    {
        int m = M[i];
        int m_past = (i == 0) ? 0 : M[i - 1];
        double error1 = l2_errors[i];
        double order1 = order(past1, error1, m_past, m); // 计算误差阶
        double error2 = linf_errors[i];
        double order2 = order(past2, error2, m_past, m);
        std::cout << std::left << std::setw(15) << N
                  << std::left << std::setw(15) << m
                  << std::left << std::setw(20) << error1
                  << std::left << std::setw(20) << order1
                  << std::left << std::setw(20) << error2
                  << std::left << std::setw(20) << order2 << std::endl;
        past1 = l2_errors[i];
        past2 = linf_errors[i];
    }
}

// 在不进行多组重复实验时，单组的误差阶是混乱的
// 这一函数用来展示单组实验的结果
template<typename Func>
void single_test(std::string name, Func func, const std::vector<int> &M, const int N)
{
    using namespace solution;
    std::cout << "---" << name << "---" << std::endl;
    std::cout << std::left << std::setw(15) << "N"
              << std::left << std::setw(15) << "M"
              << std::left << std::setw(20) << "L2_error"
              << std::left << std::setw(20) << "L2_order"
              << std::left << std::setw(20) << "L_inf_error"
              << std::left << std::setw(20) << "L_inf_order" << std::endl;
    double past1 = 0.0, past2 = 0.0;
    std::random_device rd;
    std::mt19937 generator(rd()); // 在此初始化随机数生成器，意味着同组实验使用相同的随机数序列
    ODE ode = ODE(1.0, 1.0, 0.0, 2.0, 1.0);
    auto monte_carlo = [&](int m, int n, std::mt19937 &gen)
    {
        return ODE_monte_carlo(m, n,func, gen);
    };

    auto expectation = [&](int n)
    {
        return expectation_ode(ode, n);
    };
    auto [l2_error, linf_error] = solve(monte_carlo,expectation,M, N, generator);
    for (int i = 0; i < M.size(); ++i)
    {
        int m = M[i];
        int m_past = 0;
        if (i != 0)
        {
            m_past = M[i - 1];
        }
        double error1 = l2_error[i];
        double error2 = linf_error[i];
        double order1 = order(past1, error1, m_past, m);
        double order2 = order(past2, error2, m_past, m);
        std::cout << std::left << std::setw(15) << N
                  << std::left << std::setw(15) << m
                  << std::left << std::setw(20) << error1
                  << std::left << std::setw(20) << order1
                  << std::left << std::setw(20) << error2
                  << std::left << std::setw(20) << order2 << std::endl;
        past1 = error1;
        past2 = error2;
    }
}

// 最后一个时间步的期望解
std::vector<double> expectation_pde(const solution::PDE &pde, int N)
{
    std::vector<double> true_val(N);
    const double T = pde.get_T();
    const double L = pde.get_L();
    const double dx = L / N;
    for (int j = 0; j < N; ++j)
    {
        true_val[j] = sin(j * dx + T) * sin(0.1 * T) / (0.1 * T);
    }
    return true_val;
}

/// @brief 用蒙特卡洛法求解PDE
/// @param M M个独立同分布的随机数
/// @param N 迭代次数
/// @param func 求解PDE的函数
std::vector<double> PDE_monte_carlo(int M, int N, std::function<std::vector<double>(const solution::PDE &, int, const std::string &)> func,
                                    std::mt19937 &generator)
{
    using namespace solution;
    std::uniform_real_distribution<double> uniform(-1.0, 1.0);
    ;
    const double T = 5.0;
    const double L = 2 * M_PI;
    const double cfl = 0.9;
    const double k = 1.0;
    std::vector<double> res(N, 0);
    for (int i = 0; i < M; ++i)
    {
        double a = 1 + 0.1 * uniform(generator); // a=1+0.1*Z,Z~U(0,1)
        PDE pde(a, k, T, L, cfl);
        std::vector<double> _res = func(pde, N, "");
        res += _res;
    }
    // res /= M; // 在调用处再取平均
    return res;
}

template <typename Func>
void test_ODE(const std::string &name, Func &&func, int K,
              const std::vector<int> &M, int N,
              const solution::ODE &ode = solution::ODE(1.0, 1.0, 0.0, 2.0, 1.0))
{
    auto monte_carlo = [&](int m, int n, std::mt19937 &gen)
    {
        return ODE_monte_carlo(m, n, std::forward<Func>(func), gen);
    };

    auto expectation = [&](int n)
    {
        return expectation_ode(ode, n);
    };

    repeat(name, K, M, N, monte_carlo, expectation);
}

template <typename Func>
void test_PDE(const std::string &name, Func &&func, int K,
              const std::vector<int> &M, int N,
              const solution::PDE &pde = solution::PDE(1.0, 1.0, 5.0, 2 * M_PI, 0.9))
{
    auto monte_carlo = [&](int m, int n, std::mt19937 &gen)
    {
        return PDE_monte_carlo(m, n, std::forward<Func>(func), gen);
    };

    auto expectation = [&](int n)
    {
        return expectation_pde(pde, n);
    };

    repeat(name, K, M, N, monte_carlo, expectation);
}
int main()
{
    solution::ODE ode = solution::ODE(1.0, 1.0, 0.0, 2.0, 1.0);
    test_ODE("euler",solution::euler, 500, {1000, 2000, 4000, 8000}, 1000,ode);
    test_ODE("runge-kutta",solution::runge_kutta,500,{1000,2000,4000,8000},1000,ode);
    //single_test("euler",solution::euler,{1000,2000,4000,8000},1000);
    //single_test("runge-kutta",solution::runge_kutta,{1000,2000,4000,8000},1000);
    solution::PDE pde = solution::PDE(1.0, 1.0, 5.0, 2 * M_PI, 0.9);
    test_PDE("scheme3", solution::scheme3, 500, {100, 200, 400,800}, 640, pde);
}