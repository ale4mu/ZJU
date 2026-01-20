#ifndef POLY_H
#define POLY_H

#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <fstream>
#include <string>
#include <iomanip>
#include <numeric>
#include <map>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/hermite.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/constants/constants.hpp>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <Eigen/Eigenvalues>

namespace poly
{
    const double pi = boost::math::constants::pi<double>();

    /// @brief 使用 Golub-Welsch 算法计算高斯积分的点和权重
    /// @param n 积分点数量
    /// @param alphas 雅可比矩阵的主对角线元素 (递推系数 alpha)
    /// @param betas 雅可比矩阵的次对角线元素 (递推系数 beta)
    /// @param mu0 权函数的零阶矩 (∫w(x)dx)
    /// @param points 输出的积分点
    /// @param weights 输出的积分权重
    void gauss_quadrature(
        int n,
        const std::vector<double> &alphas,
        const std::vector<double> &betas,
        double mu0,
        std::vector<double> &points,
        std::vector<double> &weights)
    {
        if (n <= 0)
            return;

        // 构建 n x n 的雅可比矩阵 J
        Eigen::MatrixXd J = Eigen::MatrixXd::Zero(n, n);
        // 主对角线元素
        for (int i = 0; i < n; ++i)
        {
            J(i, i) = alphas[i];
        }
        // 次对角线元素
        for (int i = 0; i < n - 1; ++i)
        {
            J(i, i + 1) = betas[i];
            J(i + 1, i) = betas[i];
        }

        // 求解特征值
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(J);

        // 特征值是积分点
        Eigen::VectorXd eigenvalues = eigensolver.eigenvalues();
        points.assign(eigenvalues.data(), eigenvalues.data() + eigenvalues.size());

        // 权重的计算依赖于特征向量的第一个分量
        Eigen::MatrixXd eigenvectors = eigensolver.eigenvectors();
        weights.resize(n);
        for (int i = 0; i < n; ++i)
        {
            double first = eigenvectors(0, i);
            weights[i] = mu0 * first * first;
        }
    }

    // 计算 Gauss-Legendre 积分的点和权重
    void gauss_legendre(int n, std::vector<double> &points, std::vector<double> &weights)
    {
        // Legendre 多项式的递推系数
        // p_{k+1}(x) = (x - alpha_k)p_k(x) - beta_k^2 p_{k-1}(x)
        std::vector<double> alphas(n, 0.0);
        std::vector<double> betas(n - 1);
        for (int k = 0; k < n - 1; ++k)
        {
            betas[k] = (k + 1.0) / std::sqrt(4.0 * (k + 1.0) * (k + 1.0) - 1.0);
        }
        // 权函数的零阶矩
        double mu0 = 2.0; // ∫[-1,1] 1 dx = 2

        gauss_quadrature(n, alphas, betas, mu0, points, weights);
    }

    void gauss_hermite(int n, std::vector<double> &points, std::vector<double> &weights)
    {
        // Hermite 多项式的递推系数
        std::vector<double> alphas(n, 0.0);
        std::vector<double> betas(n - 1);
        for (int k = 0; k < n - 1; ++k)
        {
            betas[k] = std::sqrt((k + 1.0) / 2.0);
        }
        // 权函数的零阶矩
        double mu0 = std::sqrt(pi); // ∫[-inf,inf] e^{-x^2} dx = sqrt(pi)

        gauss_quadrature(n, alphas, betas, mu0, points, weights);
    }

    enum class Method
    {
        USE_LIBRARY, // 调用 Boost/GSL 库
        USE_EIGEN    // 使用手动计算积分点和积分权重的版本
    };

    class Polynomial
    {
    public:
        virtual ~Polynomial() = default;

        virtual double eval(int k, double x) const = 0; // 计算第k个多项式在x处的值

        virtual std::vector<double> eval_all(int N, double x) const = 0;

        virtual double weight(double x) const = 0; // 权函数 w(x)，用于积分 ∫ f(x) g(x) w(x) dx

        virtual double norm(int k) const = 0; // 权重的归一化常数

        virtual const std::vector<double> &get_points() const = 0;

        virtual const std::vector<double> &get_weights() const = 0;

        virtual std::string name() const = 0;

        // 用于计算积分
        virtual double integrate(const std::function<double(double)> &func) const
        {
            double res = 0;
            const auto &x = get_points();
            const auto &w = get_weights();
            for (size_t i = 0; i < x.size(); ++i)
                res += w[i] * func(x[i]);
            return res;
        }
    };

    // Legendre多项式
    template <int n>
    class Legendre : public Polynomial
    {
    private:
        // inline static std::vector<double> points_lib{}, weights_lib{};
        // inline static std::vector<double> points_eigen{}, weights_eigen{};
        std::vector<double> points_lib{}, weights_lib{};
        std::vector<double> points_eigen{}, weights_eigen{};
        // static constexpr int n = 200;
        const std::vector<double> *points_ptr;
        const std::vector<double> *weights_ptr;

    public:
        explicit Legendre(Method method = Method::USE_LIBRARY)
        {
            switch (method)
            {
            case Method::USE_EIGEN:
                if (points_eigen.empty())
                {
                    gauss_legendre(n, points_eigen, weights_eigen);
                }
                points_ptr = &points_eigen;
                weights_ptr = &weights_eigen;
                break;

            case Method::USE_LIBRARY:
            default:
                if (points_lib.empty())
                {
                    // boost::math::quadrature::gauss只提供[0,1]区间上的积分点和权重
                    const auto &x_half = boost::math::quadrature::gauss<double, n>::abscissa();
                    const auto &w_half = boost::math::quadrature::gauss<double, n>::weights();
                    points_lib.reserve(n);
                    weights_lib.reserve(n);
                    for (int i = x_half.size() - 1; i >= 0; --i)
                    {
                        if (x_half[i] > 0)
                        {
                            points_lib.push_back(-x_half[i]);
                            weights_lib.push_back(w_half[i]);
                        }
                    }
                    for (size_t i = 0; i < x_half.size(); ++i)
                    {
                        points_lib.push_back(x_half[i]);
                        weights_lib.push_back(w_half[i]);
                    }
                }
                points_ptr = &points_lib;
                weights_ptr = &weights_lib;
                break;
            }
        }

        double eval(int k, double x) const override
        {
            return boost::math::legendre_p(k, x);
        }

        double weight(double x) const override
        {
            return 1.0; // w(x) = 1
        }

        double norm(int k) const override
        {
            return 2.0 / (2.0 * k + 1.0);
        }

        const std::vector<double> &get_points() const override { return *points_ptr; }

        const std::vector<double> &get_weights() const override { return *weights_ptr; }

        std::string name() const override { return "Legendre"; }

        std::vector<double> eval_all(int k, double x) const override
        {
            std::vector<double> res;
            if (k < 0)
                return res;

            res.reserve(k + 1);

            // 计算 P_0(x) 和 P_1(x)
            res.push_back(boost::math::legendre_p(0, x));
            if (k >= 1)
            {
                res.push_back(boost::math::legendre_p(1, x));
            }

            // 递推公式：res[i] = ((2.0 * i - 1.0) * x * res[i - 1] - (i - 1.0) * res[i - 2]) / i
            for (int i = 1; i < k; ++i)
            {
                res.push_back(boost::math::legendre_next(i, x, res[i], res[i - 1]));
            }
            return res;
        }
    };

    enum class HermiteType
    {
        PHYSICAL,
        PROBABILITY
    };

    template <int n>
    class Hermite : public Polynomial
    {
        // inline static std::vector<double> points_lib{}, weights_lib{};
        // inline static std::vector<double> points_eigen{}, weights_eigen{};
        // inline static std::vector<double> points_prob_lib{}, weights_prob_lib{};
        // inline static std::vector<double> points_prob_eigen{}, weights_prob_eigen{};

        std::vector<double> points_lib{}, weights_lib{};
        std::vector<double> points_eigen{}, weights_eigen{};
        std::vector<double> points_prob_lib{}, weights_prob_lib{};
        std::vector<double> points_prob_eigen{}, weights_prob_eigen{};
        const std::vector<double> *points_ptr;
        const std::vector<double> *weights_ptr;
        HermiteType type_;
        // 使用GSL库的辅助函数
        inline static gsl_integration_fixed_workspace *workspace = nullptr;
        struct WorkspaceCleaner
        {
            ~WorkspaceCleaner()
            {
                if (workspace)
                {
                    gsl_integration_fixed_free(workspace);
                    workspace = nullptr;
                }
            }
        };
        inline static WorkspaceCleaner cleaner;

        // static constexpr int n = 10;

    public:
        explicit Hermite(Method method = Method::USE_LIBRARY, HermiteType type = HermiteType::PHYSICAL) : type_(type)
        {
            switch (method)
            {
            case Method::USE_EIGEN:
                if (points_eigen.empty())
                {
                    gauss_hermite(n, points_eigen, weights_eigen);
                }
                break;

            case Method::USE_LIBRARY:
            default:
                if (points_lib.empty())
                {
                    // α = 0, b = 1, a = 0
                    workspace = gsl_integration_fixed_alloc(gsl_integration_fixed_hermite, n, 0.0, 1.0, 0.0, 0.0);
                    points_lib.assign(workspace->x, workspace->x + n);
                    weights_lib.assign(workspace->weights, workspace->weights + n);
                }
                break;
            }
            if (type_ == HermiteType::PROBABILITY)
            {
                switch (method)
                {
                case Method::USE_EIGEN:
                    if (weights_prob_eigen.empty())
                    {
                        weights_prob_eigen = weights_eigen;
                        for (auto &w : weights_prob_eigen)
                        {
                            w /= std::sqrt(pi);
                        }
                    }
                    if (points_prob_eigen.empty())
                    {
                        points_prob_eigen = points_eigen;
                        for (auto &p : points_prob_eigen)
                            p *= std::sqrt(2.0);
                    }
                    points_ptr = &points_prob_eigen;
                    weights_ptr = &weights_prob_eigen;
                    break;

                case Method::USE_LIBRARY:
                default:
                    if (weights_prob_lib.empty())
                    {
                        weights_prob_lib = weights_lib;
                        for (auto &w : weights_prob_lib)
                        {
                            w /= std::sqrt(pi);
                        }
                    }
                    if (points_prob_lib.empty())
                    {
                        points_prob_lib = points_lib;
                        for (auto &p : points_prob_lib)
                            p *= std::sqrt(2.0);
                    }
                    points_ptr = &points_prob_lib; // 节点不变
                    weights_ptr = &weights_prob_lib;
                    break;
                }
            }
            // 物理型：直接使用原数据
            else
            {
                switch (method)
                {
                case Method::USE_EIGEN:
                    points_ptr = &points_eigen;
                    weights_ptr = &weights_eigen;
                    break;
                case Method::USE_LIBRARY:
                default:
                    points_ptr = &points_lib;
                    weights_ptr = &weights_lib;
                    break;
                }
            }
        }

        double eval(int k, double x) const override
        {
            if (type_ == HermiteType::PROBABILITY)
            {
                double scaled_x = x / std::sqrt(2.0);
                double scale_factor = std::pow(2.0, -k / 2.0);
                return scale_factor * boost::math::hermite(k, scaled_x);
            }
            return boost::math::hermite(k, x);
        }

        double weight(double x) const override
        {
            if (type_ == HermiteType::PROBABILITY)
            {
                return std::exp(-x * x / 2.0) / std::sqrt(2.0 * pi);
            }
            return std::exp(-x * x); // w(x) = e^(-x^2)
        }

        // 与笔记略微不同，权函数中不含1/sqrt(π)
        double norm(int k) const override
        {
            if (type_ == HermiteType::PROBABILITY)
            {
                return boost::math::factorial<double>(k); // n!
            }
            return std::sqrt(pi) * std::pow(2.0, k) * boost::math::factorial<double>(k);
        }

        const std::vector<double> &get_points() const override { return *points_ptr; }
        const std::vector<double> &get_weights() const override { return *weights_ptr; }
        std::string name() const override { return "Hermite"; }

        std::vector<double> eval_all(int k, double x) const override
        {
            std::vector<double> res;
            if (k < 0)
                return res;

            res.reserve(k + 1);
            if (type_ == HermiteType::PROBABILITY)
            {
                res.push_back(1.0); // H0 = 1
                if (k >= 1)
                    res.push_back(x); // H1 = x
                for (int i = 1; i < k; ++i)
                {
                    // 递推公式：res[i] = x * res[i-1] - (i-1) * res[i-2]
                    res.push_back(x * res[i] - i * res[i - 1]);
                }
            }
            else
            {
                // 计算 H_0(x) 和 H_1(x)
                res.push_back(boost::math::hermite(0, x));
                if (k >= 1)
                {
                    res.push_back(boost::math::hermite(1, x));
                }

                for (int i = 1; i < k; ++i)
                {
                    // 递推公式：res[i] = 2.0 * x * res[i - 1] - 2.0 * (i - 1.0) * res[i - 2]
                    res.push_back(boost::math::hermite_next(i, x, res[i], res[i - 1]));
                }
            }
            return res;
        }
    };

    /// @brief 正交多项式投影的系数
    /// @param f 待投影的函数
    /// @param poly 基函数（正交多项式）
    /// @param N 正交多项式的最高阶数
    /// @return 投影系数
    std::vector<double> projection(
        const std::function<double(double)> &f,
        const Polynomial &poly,
        int N)
    {
        std::vector<double> coeffs(N + 1);
        const auto &points = poly.get_points();
        const auto &weights = poly.get_weights();
        std::vector<double> vals;
        for (size_t i = 0; i < points.size(); ++i)
        {
            double x = points[i];
            double w = weights[i];
            double fx = f(x);
            vals = poly.eval_all(N, x);
            for (int k = 0; k <= N; ++k)
            {
                coeffs[k] += w * fx * vals[k]; // 相当于计算内积时做的积分
            }
        }

        // 归一化
        for (int k = 0; k <= N; ++k)
        {
            coeffs[k] /= poly.norm(k);
        }
        return coeffs;
    }

    /// @brief 求网格上每个点的值
    /// @param coeffs 投影系数
    /// @param poly 正交多项式
    /// @param x_grid 网格
    /// @return 每个网格点上的投影值
    std::vector<double> evaluate(
        const std::vector<double> &coeffs,
        const Polynomial &poly,
        const std::vector<double> &x_grid)
    {
        std::vector<double> result;
        result.reserve(x_grid.size());

        for (double x : x_grid)
        {
            auto phis = poly.eval_all(coeffs.size() - 1, x);
            double sum = 0.0;
            for (size_t k = 0; k < coeffs.size(); ++k)
            {
                sum += coeffs[k] * phis[k];
            }
            result.push_back(sum);
        }

        return result;
    }

    // 计算l2误差
    double l2_error(
        const std::function<double(double)> &f,
        const std::vector<double> &coeffs,
        const Polynomial &poly)
    {
        auto residual = [&](double x)
        {
            double fx = f(x);
            auto phis = poly.eval_all(coeffs.size() - 1, x);
            double approx = 0;
            for (size_t k = 0; k < coeffs.size(); ++k)
                approx += coeffs[k] * phis[k];
            return (fx - approx) * (fx - approx);
        };
        return std::sqrt(poly.integrate(residual));
    }

    std::vector<double> create_grid(double a, double b, int num_points)
    {
        std::vector<double> grid;
        grid.reserve(num_points + 1);
        for (int i = 0; i <= num_points; ++i)
        {
            grid.push_back(a + i * (b - a) / num_points);
        }
        return grid;
    }

}
#endif // POLY_H