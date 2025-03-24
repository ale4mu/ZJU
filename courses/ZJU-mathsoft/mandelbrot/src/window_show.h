/**
 * @file window_show.h
 * @brief Header file for visualization and video output of the Mandelbrot set.
 * 
 * This file defines classes for managing the display window ("Window") and video recording ("Video").
 * It also includes utilities for 3D rendering and color mapping of the Mandelbrot set.
 */




#ifndef WINDOW_H // Prevent multiple inclusions
#define WINDOW_H

#include <opencv2/opencv.hpp>
#include <opencv2/viz.hpp>
#include <string>
#include <vector>

/**
 * @struct Point2d
 * @brief A structure representing a 2D point.
 */
struct Point2d
{
    double x = 0.0; // The x-coordinate of the point
    double y = 0.0; // The y-coordinate of the point
};

/**
 * @class Window
 * @brief A class to manage the visualization window for the Mandelbrot set.
 */
class Window
{
private:
    Point2d origin = {0.0, 0.0}; // The origin coordinates of the window, default is (0.0, 0.0)
    double dimension = 5.0;      // The dimension of the window, default is 5.0 (half of the window's width or height)
    int height = 1080;           // The height of the window in pixels, default is 1080
    int width = 1920;            // The width of the window in pixels, default is 1920

public:
    /**
     * @brief Default constructor.
     */
    Window();

    /**
     * @brief Parameterized constructor to initialize the window's origin and dimension.
     * @param ox The x-coordinate of the origin.
     * @param oy The y-coordinate of the origin.
     * @param d The dimension of the window.
     */
    Window(double ox, double oy, double d);

    /**
     * @brief Get the dimension of the window.
     * @return The dimension of the window.
     */
    double get_dimension() const;

    /**
     * @brief Get the height of the window in pixels.
     * @return The height of the window.
     */
    double get_height() const;

    /**
     * @brief Get the width of the window in pixels.
     * @return The width of the window.
     */
    double get_width() const;

    /**
     * @brief Get the x-coordinate of the window's origin.
     * @return The x-coordinate of the origin.
     */
    double get_ox() const;

    /**
     * @brief Get the y-coordinate of the window's origin.
     * @return The y-coordinate of the origin.
     */
    double get_oy() const;

    /**
     * @brief Calculate the length per pixel (LPP).
     * @return The length per pixel.
     */
    double get_lpp() const;

    /**
     * @brief Reset the dimension of the window.
     * @param new_dimension The new dimension of the window.
     */
    void set_dimension(double new_dimension);
};

/**
 * @class Video
 * @brief A class to manage video creation for Mandelbrot set visualization.
 */
class Video
{
private:
    cv::VideoWriter video; 
    std::string name; // The name of the video file
    int fourcc;       // The codec used for video encoding
    double fps;       // The frame rate of the video
    int width;       // The width of the video
    int height;      // The height of the video

public:
    /**
     * @brief Constructor to initialize the video.
     * @param _name The name of the video file.
     * @param _fourcc The codec used for video encoding.
     * @param _fps The frame rate of the video.
     * @param _width The width of the video.
     * @param _height The height of the video.
     */
    Video(const std::string& _name, int _fourcc, double _fps, int _width, int _height);

    /**
     * @brief Add a frame to the video.
     * @param frame The frame to be added.
     */
    void addFrame(const cv::Mat &frame);

    /**
     * @brief Destructor to release the video.
     */
    ~Video();
};

/**
 * @brief Render the Mandelbrot set in 3D.
 * @param image The image to render.
 */
void render3D(const cv::Mat &image);

/**
 * @brief Get the color based on the iteration count.
 * @param iter The iteration count.
 * @param max_iter The maximum iteration count.
 * @return The color corresponding to the iteration count.
 */
cv::Vec3b getColor(double iter, int max_iter);

#else // do nothing
#endif