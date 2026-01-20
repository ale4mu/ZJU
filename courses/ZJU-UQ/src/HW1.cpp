#include "solution.h"
void solve1(std::string name, std::function<std::vector<double>(const solution::ODE &, int)> func, const solution::ODE &ode, const std::vector<int> &N)
{
    using namespace solution;
    double t0 = ode.get_t0();
    double t1 = ode.get_t1();
    double u0 = ode.get_u0();

    std::cout << "---" << name << "---" << std::endl;
    std::cout << std::left << std::setw(15) << "N"
              << std::left << std::setw(20) << "L2_error"
              << std::left << std::setw(20) << "L2_order"
              << std::left << std::setw(20) << "L_inf_error"
              << std::left << std::setw(20) << "L_inf_order" << std::endl;
    double past1 = 0.0, past2 = 0.0;
    for (int i = 0; i < N.size(); ++i)
    {
        int n = N[i];
        int n_past = 0;
        if (i != 0)
        {
            n_past = N[i - 1];
        }
        auto result = func(ode, n);
        std::vector<double> T = solution_ode_series(ode, n);
        std::pair<double, double> errors = error(result, T);
        double error1 = errors.first;
        double error2 = errors.second;
        double order1 = order(past1, error1, n_past, n);
        double order2 = order(past2, error2, n_past, n);
        std::cout << std::left << std::setw(15) << n
                  << std::left << std::setw(20) << error1
                  << std::left << std::setw(20) << order1
                  << std::left << std::setw(20) << error2
                  << std::left << std::setw(20) << order2 << std::endl;
        past1 = error1;
        past2 = error2;
    }
}

void solve2(std::string name,
            std::function<std::vector<double>(const solution::PDE &, int, const std::string &)> func,
            double a, double L, double cfl, const std::vector<int> &k, const std::vector<int> &T,
            std::vector<int> N, const std::string &Filename = "")
{
    using namespace solution;
    for (int kk : k)
    {
        for (int tt : T)
        {
            std::cout << "---" << name << ",k=" << kk << ",T=" << tt << "---" << std::endl;
            std::cout << std::left << std::setw(15) << "N"
                      << std::left << std::setw(20) << "L2_error"
                      << std::left << std::setw(20) << "L2_order"
                      << std::left << std::setw(20) << "L_inf_error"
                      << std::left << std::setw(20) << "L_inf_order" << std::endl;
            double past1 = 0.0, past2 = 0.0;
            for (int i = 0; i < N.size(); ++i)
            {
                int n = N[i];
                int n_past = 0;
                if (i != 0)
                {
                    n_past = N[i - 1];
                }
                std::string filename = "";
                if (!Filename.empty())
                    filename = Filename + "_k" + std::to_string(kk) + "_T" + std::to_string(tt) + "_N" + std::to_string(n) + ".csv";

                PDE pde(a, kk, tt, L, cfl);
                auto result = func(pde, n, filename);
                std::vector<double> T = solution_pde_series(pde, n);
                std::pair<double, double> errors = error(result, T);
                double error1 = errors.first;
                double error2 = errors.second;  
                double order1 = order(past1, error1, n_past, n);
                double order2 = order(past2, error2, n_past, n);
                std::cout << std::left << std::setw(15) << n
                          << std::left << std::setw(20) << error1
                          << std::left << std::setw(20) << order1
                          << std::left << std::setw(20) << error2
                          << std::left << std::setw(20) << order2 << std::endl;
                past1 = error1;
                past2 = error2;
            }
        }
    }
}
int main()
{
    solution::ODE ode(5.0, 1.0, 0.0, 2.0, 1.0);
    std::vector<int> N = {20, 40, 80, 160};
    std::cout << "------------------------ODE------------------------" << std::endl;
    solve1("euler", solution::euler, ode, N);
    solve1("leapfrog", solution::leapfrog, ode, N);
    solve1("Runge-Kutta", solution::runge_kutta, ode, N);
    const double L = 2 * M_PI;
    const double cfl = 0.9;
    std::vector<int> k = {1, 5, 10};
    std::vector<int> T = {5, 10, 20};
    std::cout << "------------------------PDE------------------------" << std::endl;
    std::cout << "--------------------(1)--------------------" << std::endl;
    solve2("scheme1", solution::scheme1, 1.0, L, cfl, k, T, {20, 40, 80});
    std::cout << "\n\n\n";
    std::cout << "--------------------(2)--------------------" << std::endl;
    solve2("scheme2", solution::scheme2, 1.0, L, cfl, k, T, {20, 40, 80, 160});
    std::cout << "\n\n";
    solve2("scheme3", solution::scheme3, 1.0, L, cfl, k, T, {20, 40, 80, 160});
    std::cout << "\n\n\n";
    std::cout << "--------------------(3)--------------------" << std::endl;
    solve2("scheme3", solution::scheme3, -1.0, L, cfl, k, T, {10, 20, 40});
    return 0;
}