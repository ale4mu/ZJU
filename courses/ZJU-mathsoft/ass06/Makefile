PREFIX      = $(HOME)/usr/gsl
CC          = gcc
CFLAGS      = -Wall -I$(PREFIX)/include -O2
LDFLAGS     = -L$(PREFIX)/lib -lgsl -lgslcblas -lm
SRC         = gsl_test.c
OBJ         = $(SRC:.c=.o)
EXEC        = test

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $@ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(EXEC)