# Readme

# Project Purpose

This project aims to implement Conway's Game of Life, providing a visual interface that allows users to intuitively observe the game's operation. It also supports loading specific patterns, saving game states as images, and other features to enhance the game's playability and entertainment value.

# Functional Introduction

## Random Pattern Generation

The program can randomly generate the initial cell states, with each cell having a 50% chance of being alive.

## Loading Pattern Files

It supports loading specific cell patterns from RLE files and placing them at designated locations.

## Visual Interface

A text-based user interface is created using the ncurses library, which displays the cell states in real-time, with live cells represented by '*' and dead cells by spaces.

## Saving Game State as Image

The current game state can be saved as a BMP image file, with live cells in white and dead cells in black. The filename is generated based on the current time.

## Toroidal Boundary Mode 

The game supports toroidal boundary mode, meaning the world is wrap-around, and cells on the edges will "wrap" to the opposite side.

## Real-time Updates and Control

The game updates cell states in real-time, and users can control the game through key presses. Pressing 'q' exits the program, and pressing 's' saves the current state as an image.

# Usage

## Compile the Project

Run the following command in the project root directory:```make```

## Run the Program

After compilation, run the program by specifying the width and height of the game area. You can also load a specific pattern file; otherwise, it will generate a random pattern by default.

Command:```./main [width] [height] [rle_file]```

## Generate Doxygen Documentation 

Run the following command to generate Doxygen documentation: ```make doc```. This will create a doc folder containing html and latex subfolders. You can view the documentation in html/index.html or latex/refman.pdf.

## Clean Generated Files

Run ```make clean``` to remove all generated files.

# Dependencies

**This project includes stb_image_write.h, which is part of the STB library and is used for image processing.**

The function definitions are placed in the Conway.h file, and the function implementations are in the Conway.c file. The main function, test programs, and ncurses environment are in the main.c file. The ncurses library is used to create the text-based user interface, providing screen drawing and keyboard input handling functions. The stb_image_write.h is used to save image files, supporting the saving of game states in BMP format.


