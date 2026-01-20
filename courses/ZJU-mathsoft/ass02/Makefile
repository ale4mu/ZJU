# 定义编译器
CC = g++

# 定义目标文件
TARGET = main

# 定义源文件
SOURCES = main.c Conway.c

# 定义目标文件
OBJECTS = $(SOURCES:.c=.o)

# 定义库
LIBS = -lcurses

# 默认目标
all: $(TARGET)

# 编译目标文件
$(TARGET): $(OBJECTS)
	$(CC) -o $(TARGET) $(OBJECTS) $(LIBS)

doc:
	doxygen Doxyfile

# 清理生成的文件
clean:
	rm -f $(OBJECTS) $(TARGET)
	rm -rf doc 
	rm -f *.bmp
.PHONY: all clean doc 