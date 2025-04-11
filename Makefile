CC = gcc
CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors -lm
TARGET = symnmf

$(TARGET): symnmf.o utils.o sym.o norm.o diagonal.o
	$(CC) -o $(TARGET) symnmf.o utils.o sym.o norm.o diagonal.o $(CFLAGS)

symnmf.o: symnmf.c
	$(CC) -c symnmf.c $(CFLAGS)

utils.o: utils.c
	$(CC) -c utils.c $(CFLAGS)

sym.o: sym.c
	$(CC) -c sym.c $(CFLAGS)

norm.o: norm.c
	$(CC) -c norm.c $(CFLAGS)

diagonal.o: diagonal.c
	$(CC) -c diagonal.c $(CFLAGS)

clean:
	rm -f $(TARGET) *.o