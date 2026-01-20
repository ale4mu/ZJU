The whole package includes:
(1)A source code named “Calendar.cpp”.
(2)A CMakeLists.txt. 

Instructions:
(1)Open “Calendar.cpp” in your IDE , compile and run.
(2)Follow the instructions in the terminal to input the year and the number of months you want to print in one row.
(3)The year should be in the range of 1583 to 9999. The number of months you want to print in one row should be in the range of 1 to 12.If the input is not in the range, the program will throw an exception.
(4)After inputting the year and the number of months you want to print in one row correctly, you can see the calendar of the year you input.

How to use cmake:
(1)Open the terminal.
(2)Navigate to the directory where the CMakeLists.txt is located.
(3)Create a build directory and navigate to it: mkdir build && cd build.
(4)Run cmake to configure the project: cmake .. .
(5)Run the make command: make.
(6)You can use ./Calendar or make run to run the program.
(7)You can use make clean to clean the build directory.