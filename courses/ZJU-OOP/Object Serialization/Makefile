CXX := g++
CXXFLAGS := -std=c++17 -Iinclude -Wall -Wextra
LDFLAGS :=
TARGET := test

SOURCES := src/test.cpp include/tinyxml2.cpp include/base64.cpp
OBJECTS := $(SOURCES:.cpp=.o)

DATA_DIR := data

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $^ -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)
	cd $(DATA_DIR) && rm -f *.data *.xml

.PHONY: all run clean