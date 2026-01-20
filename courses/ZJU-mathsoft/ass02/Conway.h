/**
 * @file Conway.h
 * @brief Contains the data structures and function declarations for Conway's Game of Life.
 * @date 2025-03-14
 * This header file defines the data structure for the game world and the declarations of related functions.
 */

#ifndef _CONWAY_CRAZYFISH_
#define _CONWAY_CRAZYFISH_

/**
 * @struct World
 * @brief Represents the state of the game world.
 *
 * This structure contains the width, height, and a two-dimensional array representing the state of the cells.
 */
struct World
{
    int width;  ///< Width of the matrix
    int height; ///< Height of the matrix
    int **cell; ///< Two-dimensional array, where each element represents a cell. 0 means dead, 1 means alive.
};

/**
 * @var world
 * @brief The global instance of the game world.
 *
 * This global variable represents the current state of the game world.
 */
extern struct World world;

/**
 * @var boundary
 * @brief Controls the boundary mode of the world.
 *
 * 0 means the outside of the world is empty; 1 means the world is toroidal (wraps around).
 */
extern int boundary;

/**
 * @brief Initializes the cell states.
 * @param _w The width of the game area.
 * @param _h The height of the game area.
 *
 * Allocates memory for the game area and initializes all cells to the dead state.
 */
void init_cells(int _w, int _h);

/**
 * @brief Randomly generates cell states.
 *
 * Randomly generates the initial state of the cells, with each cell having a 50% chance of being alive.
 */
void set_rand_cells();

/**
 * @brief Loads a pattern file.
 * @param _x The starting x-coordinate for the pattern.
 * @param _y The starting y-coordinate for the pattern.
 * @param rle_file The path to the pattern file.
 *
 * Loads a cell pattern from an RLE (Run-Length Encoded) file and places it at the specified location.
 */
void load_pattern(int _x, int _y, char *rle_file);

/**
 * @brief Draws the cell states.
 *
 * Uses the ncurses library to draw the cell states on the screen, with live cells represented by '*' and dead cells by ' '.
 */
void draw_cells();

/**
 * @brief Counts the number of neighbors for a cell.
 * @param x The x-coordinate of the cell.
 * @param y The y-coordinate of the cell.
 * @return The number of neighbors for the cell.
 *
 * Calculates the number of neighbors for a cell based on the boundary mode (fixed or toroidal).
 */
int count_neighbors(int x, int y);

/**
 * @brief Updates the cell states.
 *
 * Updates the cell states according to the rules of Conway's Game of Life.
 */
void update_cells();

/**
 * @brief Frees the memory allocated for the cell states.
 *
 * Releases the memory allocated for the cell states.
 */
void free_cells();

/**
 * @brief Saves the current cell state as an image.
 * @param width The width of the game area.
 * @param height The height of the game area.
 * @param cell The two-dimensional array representing the cell states.
 *
 * Saves the current cell state as a BMP image, with live cells represented by white and dead cells by black.
 * The image filename is generated based on the current time, in the format "image_YYYYMMDD_HHMMSS.bmp".
 */
void save_image(int width, int height, int **cell);

#else
/// DO NOTHING
#endif