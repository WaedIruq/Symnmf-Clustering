CC = gcc
CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
LIBS = -lm

SRC = symnmf.c

OBJ = $(SRC:.c=.o)

EXECUTABLE = symnmf

INCLUDES = symnmf.h

all: $(EXECUTABLE)

 $(EXECUTABLE): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

 %.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean: rm -f $(OBJ) $(EXECUTABLE)