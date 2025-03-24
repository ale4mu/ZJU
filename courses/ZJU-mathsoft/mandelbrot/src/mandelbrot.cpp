#include "mandelbrot.h"

// 默认构造函数
Iteration::Iteration()
    : iteration_point(0.0, 0.0), iteration_times(0), max_iteration(20),
      iteration_const(0.0, 0.0), flag_stop(false), flag_convergence(false) {}

// 构造函数：初始化迭代点、最大迭代次数和迭代常数
Iteration::Iteration(std::complex<double> _ip, int _max_it, std::complex<double> _ic)
    : iteration_point(_ip), max_iteration(_max_it), iteration_const(_ic),
      iteration_times(0), flag_stop(false), flag_convergence(false) {}

// 委托构造函数：使用实部和虚部初始化迭代点和迭代常数
Iteration::Iteration(double ix, double iy, int max_it, double cx, double cy)
    : Iteration(std::complex<double>{ix, iy}, max_it, std::complex<double>{cx, cy}) {}

// 获取当前迭代点
std::complex<double> Iteration::get_iteration_point() const
{
    return iteration_point;
}

// 获取当前迭代次数
int Iteration::get_iteration_times() const
{
    return iteration_times;
}

// 获取迭代常数
std::complex<double> Iteration::get_iteration_const() const
{
    return iteration_const;
}

// 获取最大迭代次数
int Iteration::get_max_iteration() const
{
    return max_iteration;
}

// 重置迭代状态
void Iteration::reset(std::complex<double> _ip, std::complex<double> _ic, int max_iter)
{
    iteration_point = _ip;
    iteration_const = _ic;
    iteration_times = 0;
    flag_stop = false;
    flag_convergence = false;
    max_iteration = max_iter;
}

// 实现迭代步骤
void Mandelbrot::forward_step()
{
    iteration_point = iteration_point * iteration_point + iteration_const;               // Mandelbrot 迭代公式
    iteration_times++;                                                                   // 增加迭代次数
    flag_stop = (std::abs(iteration_point) > 2.0) || (iteration_times == max_iteration); // 停止条件
    flag_convergence = (iteration_times == max_iteration);                               // 达到最大迭代次数可近似认为收敛
}

// 实现停止条件判断
bool Mandelbrot::stop_criterion() const
{
    return flag_stop;
}

// 实现发散或收敛判断
bool Mandelbrot::is_convergence() const
{
    return flag_convergence;
}

// 生成mandelbrot集图像函数
void generateMandelbrotImage(cv::Mat &image, const std::vector<double> &x_coords, const std::vector<double> &y_coords, int max_iter)
{
    for (int j = 0; j < image.rows; j++)
    {
        cv::Vec3b *pixel = image.ptr<cv::Vec3b>(j); // 使用指针加速像素访问
        double y = y_coords[j];

        for (int i = 0; i < image.cols; i++)
        {
            double x = x_coords[i];

            Mandelbrot man(std::complex<double>{0.0, 0.0}, max_iter, std::complex<double>{x, y}); // 重置mandelbrot类状态
            while (!man.stop_criterion())                                                         // 不满足停止条件时一直迭代
                man.forward_step();

            if (man.is_convergence())
            {
                pixel[i] = cv::Vec3b(255, 255, 255); // Mandelbrot 集内为白色
            }
            else
            {
                double iter = man.get_iteration_times();
                pixel[i] = getColor(iter, max_iter); // 使用 HSV 颜色映射
            }
        }
    }
}

// 计算坐标
void calculateCoords(std::vector<double> &x_coords, std::vector<double> &y_coords, double ox, double oy, double dim, double lpp, int width, int height, double zoom_factor)
{
    for (int i = 0; i < width; ++i)
        x_coords[i] = ox - dim * zoom_factor + lpp * i;

    for (int j = 0; j < height; ++j)
        y_coords[j] = oy - dim * zoom_factor * height / width + lpp * j;
}