/**
 * @file mandelbrot.h
 * @brief Header file for Mandelbrot set computation and visualization.
 *
 * This file defines the "Iteration" base class and "Mandelbrot" derived class for computing the Mandelbrot set.
 * It provides functions for generating Mandelbrot images and calculating coordinates.
 */




#ifndef MANDELBROT_H // Prevent multiple inclusions
#define MANDELBROT_H

#include <complex>
#include <cmath>
#include <opencv2/opencv.hpp>
#include <vector>
#include "window_show.h"

/**
 * @class Iteration
 * @brief Base class for iteration, defining basic properties and interfaces for Mandelbrot set computation.
 */
class Iteration
{
protected:
    std::complex<double> iteration_point; // Current complex point in iteration
    int iteration_times;                  // Current iteration count
    int max_iteration;                    // Maximum iteration count
    std::complex<double> iteration_const; // Iteration constant (e.g., c in Mandelbrot set)
    bool flag_stop;                       // Whether the stopping condition is met
    bool flag_convergence;                // Whether the iteration converges

public:
    /**
     * @brief Default constructor.
     */
    Iteration();

    /**
     * @brief Default destructor.
     */
    virtual ~Iteration() = default;

    /**
     * @brief Constructor to initialize iteration point, maximum iteration count, and iteration constant.
     * @param _ip The initial iteration point.
     * @param _max_it The maximum number of iterations.
     * @param _ic The iteration constant.
     */
    Iteration(std::complex<double> _ip, int _max_it, std::complex<double> _ic);

    /**
     * @brief Delegating constructor to initialize iteration point and constant using real and imaginary parts.
     * @param ix The real part of the iteration point.
     * @param iy The imaginary part of the iteration point.
     * @param max_it The maximum number of iterations.
     * @param cx The real part of the iteration constant.
     * @param cy The imaginary part of the iteration constant.
     */
    Iteration(double ix, double iy, int max_it, double cx, double cy);

    /**
     * @brief Get the current iteration point.
     * @return The current iteration point.
     */
    std::complex<double> get_iteration_point() const;

    /**
     * @brief Get the current iteration count.
     * @return The current iteration count.
     */
    int get_iteration_times() const;

    /**
     * @brief Get the iteration constant.
     * @return The iteration constant.
     */
    std::complex<double> get_iteration_const() const;

    /**
     * @brief Get the maximum iteration count.
     * @return The maximum iteration count.
     */
    int get_max_iteration() const;

    /**
     * @brief Pure virtual function to perform one iteration step.
     */
    virtual void forward_step() = 0;

    /**
     * @brief Pure virtual function to check if the stopping condition is met.
     * @return True if the stopping condition is met, false otherwise.
     */
    virtual bool stop_criterion() const = 0;

    /**
     * @brief Pure virtual function to check if the iteration converges.
     * @return True if the iteration converges, false otherwise.
     */
    virtual bool is_convergence() const = 0;

    /**
     * @brief Reset the iteration state.
     * @param _ip The new iteration point.
     * @param _ic The new iteration constant.
     * @param max_iter The new maximum iteration count.
     */
    void reset(std::complex<double> _ip, std::complex<double> _ic, int max_iter);
};

/**
 * @class Mandelbrot
 * @brief Derived class for Mandelbrot set iteration logic.
 */
class Mandelbrot : public Iteration
{
public:
    using Iteration::Iteration; // Inherit constructors

    /**
     * @brief Perform one iteration step.
     */
    virtual void forward_step() override;

    /**
     * @brief Check if the stopping condition is met.
     * @return True if the stopping condition is met, false otherwise.
     */
    virtual bool stop_criterion() const override;

    /**
     * @brief Check if the iteration diverges or converges.
     * @return True if the iteration converges, false otherwise.
     */
    virtual bool is_convergence() const override;
};

/**
 * @brief Generate a Mandelbrot set image.
 * @param image The output image.
 * @param x_coords The x-coordinates for the image.
 * @param y_coords The y-coordinates for the image.
 * @param max_iter The maximum number of iterations.
 */
void generateMandelbrotImage(cv::Mat &image, const std::vector<double> &x_coords, const std::vector<double> &y_coords, int max_iter);

/**
 * @brief Calculate coordinates for the Mandelbrot set.
 * @param x_coords The vector to store x-coordinates.
 * @param y_coords The vector to store y-coordinates.
 * @param ox The origin x-coordinate.
 * @param oy The origin y-coordinate.
 * @param dim The dimension of the window.
 * @param lpp The length per pixel.
 * @param width The width of the image.
 * @param height The height of the image.
 * @param zoom_factor The zoom factor (default is 1.0).
 */
void calculateCoords(std::vector<double> &x_coords, std::vector<double> &y_coords, double ox, double oy, double dim, double lpp, int width, int height, double zoom_factor = 1.0);

#else // do nothing
#endif