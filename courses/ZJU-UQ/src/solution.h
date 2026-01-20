#ifndef SOLUTION_H
#define SOLUTION_H

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>
#include <numeric>
#include <utility>
#include <algorithm>
#include <iterator>
#include <functional>
#include <fstream>

namespace solution
{
    template <typename T>
    std::vector<T> &operator+=(std::vector<T> &l, const std::vector<T> &r)
    {
        if (l.size() != r.size())
        {
            throw std::invalid_argument("Vectors must be of the same size.");
        }
        for (int i = 0; i < l.size(); ++i)
        {
            l[i] += r[i];
        }
        return l;
    }

    template <typename T>
    std::vector<T> &operator/=(std::vector<T> &l, const int M)
    {
        if (M == 0)
        {
            throw std::invalid_argument("M cannot be 0.");
        }
        for (int i = 0; i < l.size(); ++i)
        {
            l[i] /= M;
        }
        return l;
    }

    class ODE
    {
    public:
        ODE(double _alpha, double _beta, double _t0, double _t1, double _u0) : alpha(_alpha), beta(_beta), t0(_t0), t1(_t1), u0(_u0) {}

        // ODE: du/dt = f(t, u) = -alpha * u
        double f(double t, double u) const
        {
            return -alpha * u;
        }

        // ODE的精确解:beta*exp{-alpha*t}
        double solution_ode(double t) const
        {
            return beta * exp(-alpha * t);
        }

        double get_t0() const { return t0; }
        double get_t1() const { return t1; }
        double get_u0() const { return u0; } // return beta?
        double get_alpha() const { return alpha; }
        double get_beta() const { return beta; }

        double get_dt(int N) const {return (t1 - t0) / N;}

    private:
        double alpha;
        double beta;
        double t0;
        double t1;
        double u0;
    };

    class PDE
    {
    public:
        PDE(double _a, double _k, double _T, double _L, double _mcl) : a(_a), k(_k), T(_T), L(_L), cfl(_mcl) {}

        double initial(double x) const
        {
            return sin(k * x);
        }
        double solution_pde(double x, double t) const
        {
            return sin(k * (x + a * t));
        }

        double get_a() const { return a; }
        double get_k() const { return k; }
        double get_T() const { return T; }
        double get_L() const { return L; }
        double get_cfl() const { return cfl; }
        double get_dx(int N) const {return L / N;}

    private:
        double a;
        double k;
        double T;
        double L;
        double cfl; // Δt/Δx
    };

    // L2误差
    inline double L2Error(const std::vector<double> &a, const std::vector<double> &b)
    {
        if (a.size() != b.size())
        {
            throw std::invalid_argument("Vectors must have the same size.");
        }
        int N = a.size();
        double sum_sq = 0.0;
        for (size_t i = 0; i < N; ++i)
        {
            double diff = a[i] - b[i];
            sum_sq += diff * diff;
        }
        return std::sqrt(sum_sq / N);
    }

    // L无穷误差(切比雪夫距离)
    inline double L_Chebyshev(const std::vector<double> &a, const std::vector<double> &b)
    {
        if (a.size() != b.size())
        {
            throw std::invalid_argument("Vectors must have the same size.");
        }
        double max_diff = 0.0;
        for (size_t i = 0; i < a.size(); ++i)
        {
            double diff = std::abs(a[i] - b[i]);
            if (diff > max_diff)
            {
                max_diff = diff;
            }
        }
        return max_diff;
    }

    // inline std::vector<double> solution_ode_series(const ODE &ode, int n)
    // {
    //     std::vector<double> T(n); // 精确解
    //     double t0 = ode.get_t0();
    //     double t1 = ode.get_t1();
    //     double dt = (t1 - t0) / n;
    //     for (int i = 0; i < n; i++)
    //     {
    //         T[i] = ode.solution_ode(t0 + (i + 1) * dt);
    //     }
    //     return T;
    // }

    inline std::vector<double> solution_ode_series(const ODE &ode, int n)
    {
        if (n <= 0)
            return {};

        double t0 = ode.get_t0();
        double t1 = ode.get_t1();
        double dt = (t1 - t0) / n; // 每步长度

        std::vector<double> T(n);

        for (int i = 0; i < n; ++i)
        {
            double t = t0 + (i + 1) * dt;

            if (i == n - 1)
                t = t1;

            T[i] = ode.solution_ode(t);
        }

        return T;
    }
    inline std::pair<double, double> error(std::vector<double> &U, std::vector<double> &T)
    {
        return std::make_pair(L2Error(U, T), L_Chebyshev(U, T));
    }

    // 误差阶
    inline double order(double error1, double error2, int N1, int N2)
    {
        if (error1 == 0)
        {
            return 0.0;
        }
        return std::log(error1 / error2) / std::log(static_cast<double>(N2) / N1);
    }

    // 向前插商（欧拉法）
    inline std::vector<double> euler(const ODE &ode, int n)
    {
        double t0 = ode.get_t0();
        double t1 = ode.get_t1();
        double u0 = ode.get_u0();
        double dt = (t1 - t0) / n;
        double u = u0;
        std::vector<double> U(n); // 存储数值解
        std::vector<double> T(n); // 存储精确解
        for (int i = 0; i < n; ++i)
        {
            u = u + dt * ode.f(t0 + i * dt, u); // (u_n+1 - u_n) / dt约等于f(t_n,u_n)
            U[i] = u;
        }
        return U;
    }

    // 中心差商
    inline std::vector<double> leapfrog(const ODE &ode, int n)
    {

        double t0 = ode.get_t0();
        double t1 = ode.get_t1();
        double u0 = ode.get_u0();

        double dt = (t1 - t0) / n;
        double u_prev = u0;

        // 使用欧拉法计算第一步
        double u_curr = u_prev + dt * ode.f(t0, u_prev);

        std::vector<double> U(n); // 数值解
        U[0] = u_curr;
        for (int i = 1; i < n; ++i)
        {
            double u_next = u_prev + 2 * dt * ode.f(t0 + i * dt, u_curr); //(u_n+1 - u_n-1) / 2dt约等于f(t_n,u_n)
            u_prev = u_curr;
            u_curr = u_next;
            U[i] = u_curr;
        }
        return U;
    }

    // 三阶Runge-Kutta方法
    inline std::vector<double> runge_kutta(const ODE &ode, int n)
    {
        double t0 = ode.get_t0();
        double t1 = ode.get_t1();
        double u0 = ode.get_u0();

        double dt = (t1 - t0) / n;
        std::vector<double> U(n); // 数值解
        for (int i = 0; i < n; ++i)
        {
            double y1 = u0 + dt * ode.f(t0 + i * dt, u0);
            double y2 = 0.75 * u0 + 0.25 * y1 + 0.25 * dt * ode.f(t0 + (i + 1) * dt, y1);
            u0 = u0 / 3.0 + 2.0 * y2 / 3.0 + 2.0 * dt * ode.f(t0 + (i + 0.5) * dt, y2) / 3.0;
            U[i] = u0;
        }
        return U;
    }

    /// @brief 保存图像的辅助函数
    /// @param filename 文件名
    /// @param numerical 数值解
    /// @param exact 精确解
    /// @param N
    inline void save(const std::string &filename, const std::vector<double> &numerical, const std::vector<double> &exact, int N)
    {
        std::ofstream outFile(filename);
        if (!outFile.is_open())
        {
            std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
            return;
        }

        const double L = 2 * M_PI;
        const double dx = L / N;
        outFile << "x,numerical_u,exact_u\n";

        outFile << std::fixed << std::setprecision(8);
        for (int j = 0; j < N; ++j)
        {
            double x = j * dx;
            outFile << x << "," << numerical[j] << "," << exact[j] << "\n";
        }

        // 添加一个点回到起点，让周期性图像闭合
        outFile << L << "," << numerical[0] << "," << exact[0] << "\n";
        outFile.close();
    }

    inline std::vector<double> solution_pde_series(const PDE &pde, int N)
    {
        std::vector<double> true_val(N);
        const double T = pde.get_T();
        const double L = pde.get_L();
        const double dx = L / N;
        for (int j = 0; j < N; ++j)
        {
            true_val[j] = pde.solution_pde(j * dx, T);
        } // 精确解
        return true_val;
    }

    using UpdateRule = std::function<double(const std::vector<double> &, int, int, double, double, double)>;

    // 通用的求解器
    inline std::vector<double> scheme(const PDE &pde, int N, const UpdateRule &rule, const std::string &Filename = "")
    {
        const double a = pde.get_a();
        const double k = pde.get_k();
        const double T = pde.get_T();
        const double L = pde.get_L();
        const double cfl = pde.get_cfl();
        const double dx = L / N;
        const double dt = cfl * dx;
        const int num_steps = static_cast<int>(floor(T / dt));
        const double final = T - num_steps * dt;

        std::vector<double> u_current(N);
        std::vector<double> u_new(N);
        std::vector<double> true_val(N);

        for (int i = 0; i < N; ++i)
        {
            u_current[i] = pde.initial(i * dx); // Xi = i * dx, U(x,0) = sin(kx)
        }

        for (int n = 0; n < num_steps; ++n)
        {
            for (int j = 0; j < N; ++j)
            {
                u_new[j] = rule(u_current, j, N, a, dt, dx);
            }
            u_current = u_new;
        }

        // 最后一个时间步可能没走完
        if (final > 1e-12)
        {
            for (int j = 0; j < N; ++j)
            {
                u_new[j] = rule(u_current, j, N, a, final, dx);
            }
            u_current = u_new;
        }

        for (int j = 0; j < N; ++j)
        {
            true_val[j] = pde.solution_pde(j * dx, T);
        } // 精确解

        if (!Filename.empty())
        {
            save(Filename, u_current, true_val, N);
        }

        return u_current;
    }

    // scheme1
    inline std::vector<double> scheme1(const PDE &pde, int N, const std::string &Filename = "")
    {
        auto rule1 = [](const std::vector<double> &u, int j, int N, double a, double dt, double dx)
        {
            int j_prev = (j - 1 + N) % N;
            int j_next = (j + 1) % N;
            return u[j] + a * (dt / (2.0 * dx)) * (u[j_next] - u[j_prev]);
        };
        return scheme(pde, N, rule1, Filename);
    }

    // scheme2
    inline std::vector<double> scheme2(const PDE &pde, int N, const std::string &Filename = "")
    {
        auto rule2 = [](const std::vector<double> &u, int j, int N, double a, double dt, double dx)
        {
            int j_prev = (j - 1 + N) % N;
            int j_next = (j + 1) % N;
            return 0.5 * (u[j_next] + u[j_prev]) + a * (dt / (2.0 * dx)) * (u[j_next] - u[j_prev]);
        };
        return scheme(pde, N, rule2, Filename);
    }

    // scheme3
    inline std::vector<double> scheme3(const PDE &pde, int N, const std::string &Filename = "")
    {
        auto rule3 = [](const std::vector<double> &u, int j, int N, double a, double dt, double dx)
        {
            int j_prev = (j - 1 + N) % N;
            int j_next = (j + 1) % N;
            // if (a > 0)
            return u[j] + a * (dt / dx) * (u[j_next] - u[j]);
            // else // a<0
            //     return u[j] + a * (dt / dx) * (u[j] - u[j_prev]);
        };
        return scheme(pde, N, rule3, Filename);
    }
}

#endif // SOLUTION_H