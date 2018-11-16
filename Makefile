# General flags
CC = mpicc
CFLAGS = -Wall -std=c99 -O3 -g -fopenmp
LDFLAGS = -lm

BIN = build/moravia

OBJS = build/moravia.o build/mmio.o build/graph.o
HEADERS = src/mmio.h src/graph.h
SOURCES = src/moravia.c src/matrix.c mmio.c 

all:		$(BIN)

$(BIN):		$(OBJS)
		$(CC) $(CFLAGS) $(OBJS) -o $(BIN) $(LDFLAGS)

build/%.o:	src/%.c $(HEADERS)
		mkdir -p build
		$(CC) $(CFLAGS) -c $< -o $@

clean:
		rm -rf build

