#include <iostream> // 用于输入输出操作
#include <cstdlib>  // 用于标准库函数（如 exit）
#include <vector>
#include "window_show.h"      // 包含 Window 类的头文件，用于管理窗口和图像参数
#include "mandelbrot.h"       // 包含 Mandelbrot 迭代逻辑的头文件
#include <opencv2/opencv.hpp> // 使用 OpenCV 生成 PNG 图像
#include <omp.h>              // 引入OpenMP进行并行化
#include <opencv2/viz.hpp>    //3D可视化

// 判断要生成的图像格式是否合法(.png/.jpg/.bmp)
bool isSupportedImageFormat(const std::string &filename);

int main(int argc, char *argv[])
{
    // 检查命令行参数是否足够
    if (argc < 5 || !(isSupportedImageFormat(argv[1])) || argc > 6)
    {
        std::cerr << "Usage: " << argv[0] << " filename(.jpg /.png /.bmp) ox oy dimension [generate video?]" << std::endl;
        exit(-1); // 如果参数不足，输出使用说明并退出程序
    }

    // 使用命令行参数初始化 Window 对象
    Window win(std::atof(argv[2]), std::atof(argv[3]), std::atof(argv[4]));
    double dim = win.get_dimension();  // 获取窗口的尺寸
    int width = win.get_width();       // 获取窗口的宽度（像素）
    int height = win.get_height();     // 获取窗口的高度（像素）
    int max_iter = 20;                 // 设置最大迭代次数
    int total_frames = 100;            // 总帧数
    bool generate_video = (argc == 6); // 判断是否生成视频

    if (generate_video)
    {
        std::vector<cv::Mat> frame_buffer;                                           // 缓存帧
        frame_buffer.resize(total_frames);                                           // 设置frame_buffer 的大小，初始化为默认值
        omp_set_num_threads(omp_get_max_threads());                                  // 动态设置线程数
#pragma omp parallel for schedule(dynamic) ordered private(win) shared(frame_buffer) // 使用OpenMP并行化
        for (int frame = 0; frame < total_frames; ++frame)
        {
            int num_threads = omp_get_num_threads();
            std::cout << "Number of threads: " << num_threads << std::endl;
            // 在每一帧中调整Mandelbrot集的缩放比例
            double zoom_factor = 1.0 - frame * 0.01; // 逐渐缩小
            std::vector<double> x_coords(width);
            std::vector<double> y_coords(height);
            int local_max_iter = max_iter + frame;                                                                           // 每个线程使用独立的 max_iter 值，不断增加最大迭代次数观察变化
            win.set_dimension(dim * zoom_factor);                                                                            // 设置新的窗口尺寸
            cv::Mat image(height, width, CV_8UC3, cv::Scalar(0, 0, 0));                                                      // 创建image图像，初始化为黑色背景
            calculateCoords(x_coords, y_coords, win.get_ox(), win.get_oy(), dim, win.get_lpp(), width, height, zoom_factor); // 计算坐标
            generateMandelbrotImage(image, x_coords, y_coords, local_max_iter);                                              // 生成mandelbrot图像
#pragma omp ordered
            {
                frame_buffer[frame] = image.clone(); // 按顺序填充帧
            }
        }

        cv::imwrite(argv[1], frame_buffer[0]);      // 保存为 PNG/JPG/BMP 图像
        cv::Mat loaded_image = cv::imread(argv[1]); // 读取保存的图像
        render3D(loaded_image);                     // 传递给 render3D 函数
        // 使用 MP4 格式生成视频
        Video video("mandelbrot.mp4", cv::VideoWriter::fourcc('m', 'p', '4', 'v'), 15, width, height);

        // 批量写入视频
        for (const auto &frame : frame_buffer)
        {
            video.addFrame(frame);
        }
    }
    else
    {
        std::vector<double> x_coords(width);
        std::vector<double> y_coords(height);
        cv::Mat image(height, width, CV_8UC3, cv::Scalar(0, 0, 0));                                         // 创建image图像，初始化为黑色背景
        calculateCoords(x_coords, y_coords, win.get_ox(), win.get_oy(), dim, win.get_lpp(), width, height); // 计算坐标
        generateMandelbrotImage(image, x_coords, y_coords, max_iter);                                       // 生成mandelbrot图像
        cv::imwrite(argv[1], image);                                                                        // 保存为 PNG/JPG/BMP 图像
        cv::Mat loaded_image = cv::imread(argv[1]);                                                         // 读取保存的图像
        render3D(loaded_image);                                                                             // 传递给 render3D 函数
    }
    return 0;
}

bool isSupportedImageFormat(const std::string &filename) // 判断文件名是否合法
{
    // 找到最后一个点号的位置
    size_t dotPos = filename.find_last_of('.');

    // 如果没有找到点号，或者点号是第一个字符，返回 false
    if (dotPos == std::string::npos || dotPos == 0)
    {
        return false;
    }

    // 提取文件扩展名
    std::string extension = filename.substr(dotPos + 1);

    // 将扩展名转换为小写，方便比较
    for (char &c : extension)
    {
        c = tolower(c);
    }

    // 检查扩展名是否是支持的格式
    if (extension == "bmp" || extension == "jpg" || extension == "png")
    {
        return true;
    }

    return false;
}
