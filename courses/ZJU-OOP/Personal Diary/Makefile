CXX = g++
CXXFLAGS = -std=c++17 -Wall -Iinclude

SRC_DIR = src
BUILD_DIR = build
INCLUDE_DIR = include

SOURCES = $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS = $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SOURCES))
EXECUTABLES = pdadd pdlist pdshow pdremove

all: $(EXECUTABLES)

pdadd: $(BUILD_DIR)/pdadd.o $(BUILD_DIR)/Diary.o
	$(CXX) $(CXXFLAGS) $^ -o $@ 

pdlist: $(BUILD_DIR)/pdlist.o $(BUILD_DIR)/Diary.o
	$(CXX) $(CXXFLAGS) $^ -o $@ 

pdshow: $(BUILD_DIR)/pdshow.o $(BUILD_DIR)/Diary.o
	$(CXX) $(CXXFLAGS) $^ -o $@ 

pdremove: $(BUILD_DIR)/pdremove.o $(BUILD_DIR)/Diary.o
	$(CXX) $(CXXFLAGS) $^ -o $@ 

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

test: all
	./test.sh
	
clean:
	rm -rf $(BUILD_DIR) $(EXECUTABLES) diary.txt

.PHONY: all clean