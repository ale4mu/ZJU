CC = gcc
LDFLAGS = -lm

SRC = src/mnist_cnn.c
EXEC = test     

all: $(EXEC)

$(EXEC): $(SRC)
	$(CC) $(SRC) -o $(EXEC) $(LDFLAGS)


run:
	./$(EXEC) images/test_0.bmp
	./$(EXEC) images/test_1.bmp
	./$(EXEC) images/test_2.bmp
	./$(EXEC) images/test_3.bmp
	./$(EXEC) images/test_4.bmp
	./$(EXEC) images/test_5.bmp
	./$(EXEC) images/test_6.bmp
	./$(EXEC) images/test_7.bmp
	./$(EXEC) images/test_8.bmp
	./$(EXEC) images/test_9.bmp

	
clean:
	rm -f $(EXEC)
