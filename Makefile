# General flags
CC = mpicc
CFLAGS = -Wall -O3 -g -fopenmp
LDFLAGS = -lm

BIN = build/moravia

OBJS = build/moravia.o build/mmio.o build/graph.o build/timing.o
HEADERS = src/mmio.h src/graph.h src/timing.h
SOURCES = src/moravia.c src/matrix.c src/mmio.c src/timing.c

all:		$(BIN)

$(BIN):		$(OBJS)
		$(CC) $(CFLAGS) $(OBJS) -o $(BIN) $(LDFLAGS)

build/%.o:	src/%.c $(HEADERS)
		mkdir -p build
		$(CC) $(CFLAGS) -c $< -o $@

clean:
		rm -rf build

