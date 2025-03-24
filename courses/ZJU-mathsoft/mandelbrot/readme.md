# Mandelbrot集
# 项目目的

用C++生成并可视化 Mandelbrot 集，直观地感受其性质。

# 功能

## 生成 Mandelbrot 集图像

用户可以指定图像的格式（支持 PNG、JPG 和 BMP）、原点坐标、窗口尺寸等参数，生成 Mandelbrot 集的二维图像。

## 3D 可视化

通过 OpenCV 的 3D 可视化功能，将 Mandelbrot 集的图像转换为 3D 点云，并以 3D 形式展示，为用户带来全新的视觉体验。

## 视频生成

用户可以选择生成MP4格式的视频文件。视频画面会不断放大并增加最大迭代次数，用户可以观察到 Mandelbrot 集在不同缩放比例下的变化。

# 使用方法

- 安装OpenCV、OpenMP等。
- 进入mandelbrot文件夹，可以看到有doc和src两个子文件夹.
- 打开终端编译：```make```
- 进入src文件夹：```cd src```
- 只生成图像：```./main filename(.jpg/.png/.bmp) ox oy dimension```
- 生成图像和视频文件：```./main filename ox oy dimension [generate_video]```
- 运行测试程序：```./test```
- src目录下将生成mandelbrot_test.png图像.
- 在mandelbrot目录下编译：```make doc```
- 进入doc文件夹：```cd doc```
- doc中将生成report.pdf文件和html文件夹，打开html文件夹中的index.html即可查看doxygen文档.
- 删除所有生成的文件：```make clean```

# 依赖关系

## OpenCV

用于图像处理和 3D 可视化。

## OpenMP

用于并行化计算 Mandelbrot 集的迭代过程，提高程序的运行效率。

## mandelbrot.h/mandelbrot.cpp

定义和实现mandelbrot集的类和函数。

## window_show.h/window_show.cpp

定义和实现与图像、视频显示相关的类和函数。

## main.cpp

主程序，实现mandelbrot集的迭代，生成最终的图像和视频。

# 加分项

- [√] 早于第六周提交; 
- [√] 代码结构更加清晰, 有良好的注释和文档;
- [√] 有良好的项目管理, 使用 Makefile, CMake等构建工具（请在文档README.md中说明）;  *Makefile*
- [√] 采用其他图片格式增加了展示效果;
- [√] 提供色彩更加丰富的可视化结果;
- [√] 提供展现局部方法的动态图形或视频;
- [√] 其他更加丰富的可视化表现，如 3D 渲染;
- [√] \*数学上对 Mandelbrot 集进行分析和探索, 增强 report 的学术份量;
- [√] \*提升计算效率，例如利用并行计算，GPU 加速等;  *OpenMP*
- [ ] \*其他创意、数学、算法和技术上的提升。
