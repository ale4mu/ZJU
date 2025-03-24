#include "window_show.h"

// 默认构造函数
Window::Window() {};

// 带参数的构造函数，用于初始化窗口的原点坐标和尺寸
Window::Window(double ox, double oy, double d)
    : origin{ox, oy}, dimension(d) {}

// 获取窗口的尺寸
double Window::get_dimension() const { return dimension; }

// 获取窗口的高度（像素）
double Window::get_height() const { return height; }

// 获取窗口的宽度（像素）
double Window::get_width() const { return width; }

// 获取窗口原点的 x 坐标
double Window::get_ox() const { return origin.x; }

// 获取窗口原点的 y 坐标
double Window::get_oy() const { return origin.y; }

// 计算每个像素对应的长度（像素比例，Length Per Pixel）
double Window::get_lpp() const
{
    return (dimension * 2.0 / width); // 窗口总宽度除以像素宽度，得到每个像素的长度
}

// 重新设置窗口的尺寸
void Window::set_dimension(double new_dimension)
{
    dimension = new_dimension;
}

//构造函数
Video::Video(const std::string &_name, int _fourcc, double _fps, int _width, int _height)
    : name(_name), fourcc(_fourcc), fps(_fps), width(_width), height(_height)
{
    video.open(name, fourcc, fps, cv::Size(width, height));
    if (!video.isOpened())
    {
        throw std::runtime_error("Error: Could not open video file for writing.");
    }
}

//写入帧
void Video::addFrame(const cv::Mat &frame)
{
    if (!video.isOpened()) {
        throw std::runtime_error("Error: Could not open video file for writing.");
    }
    video.write(frame);
}

//释放视频
Video::~Video()
{
    video.release();
}

// 3D渲染
void render3D(const cv::Mat &image)
{
    cv::viz::Viz3d window("Mandelbrot 3D");
    cv::Mat points(image.rows * image.cols, 1, CV_32FC3);
    cv::Mat colors(image.rows * image.cols, 1, CV_8UC3);

    for (int j = 0; j < image.rows; j++)
    {
        for (int i = 0; i < image.cols; i++)
        {
            cv::Vec3b color = image.at<cv::Vec3b>(j, i);
            points.at<cv::Vec3f>(j * image.cols + i) = cv::Vec3f(i, j, color[0] / 255.0f * 10);
            colors.at<cv::Vec3b>(j * image.cols + i) = color;
        }
    }

    cv::viz::WCloud cloud(points, colors);
    window.showWidget("Mandelbrot", cloud);
    window.spin();
}

// 计算颜色
cv::Vec3b getColor(double iter, int max_iter)
{
    double hue = iter / max_iter * 360;         // 色调
    double saturation = 1.0;                    // 饱和度
    double value = iter < max_iter ? 1.0 : 0.0; // 明度

    double c = value * saturation;
    double x = c * (1 - abs(fmod(hue / 60, 2) - 1));
    double m = value - c;

    double r = 0, g = 0, b = 0;
    if (hue < 60)
    {
        r = c;
        g = x;
        b = 0;
    }
    else if (hue < 120)
    {
        r = x;
        g = c;
        b = 0;
    }
    else if (hue < 180)
    {
        r = 0;
        g = c;
        b = x;
    }
    else if (hue < 240)
    {
        r = 0;
        g = x;
        b = c;
    }
    else if (hue < 300)
    {
        r = x;
        g = 0;
        b = c;
    }
    else
    {
        r = c;
        g = 0;
        b = x;
    }

    return cv::Vec3b((b + m) * 255, (g + m) * 255, (r + m) * 255);
}
