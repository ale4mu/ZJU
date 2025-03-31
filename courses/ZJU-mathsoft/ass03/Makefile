CXX := g++
CXXFLAGS := -std=c++11 -I../include
LDFLAGS := -lgmpxx -lgmp
SRC_DIR := src
BIN := $(SRC_DIR)/test

all: $(BIN)

$(BIN): include/Function.hpp include/LinearFunction.hpp include/Polynomial.hpp $(SRC_DIR)/test.cpp
	cd $(SRC_DIR) && $(CXX) $(CXXFLAGS) test.cpp $(LDFLAGS) -o test

clean:
	cd $(SRC_DIR) && rm -f test test.exe

.PHONY: all clean