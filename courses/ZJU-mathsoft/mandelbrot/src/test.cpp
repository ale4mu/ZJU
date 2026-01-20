#include <iostream>
#include <opencv2/opencv.hpp>
#include "window_show.h"
#include "mandelbrot.h"

int main()
{
    // 设置图像尺寸
    int height = 1080;
    int width = height * 3 / 2;

    int max_iter = 20;
    cv::Mat image(height, width, CV_8UC3, cv::Scalar(0, 0, 0));

    // 手动设置x坐标：从-1.5到1.5，共width个点
    std::vector<double> x_coords(width);
    double x_min = -1.5;
    double x_max = 1.5;
    for (int i = 0; i < width; ++i)
    {
        x_coords[i] = x_min + (x_max - x_min) * i / (width - 1);
    }

    // 手动设置y坐标
    std::vector<double> y_coords(height);
    double y_min = -1.0;
    double y_max = 1.0;
    for (int j = 0; j < height; ++j)
    {
        y_coords[j] = y_max - (y_max - y_min) * j / (height - 1);
    }

    // 生成Mandelbrot集图像
    generateMandelbrotImage(image, x_coords, y_coords, max_iter);

    // 保存图像
    cv::imwrite("mandelbrot_test.png", image);

    return 0;
}