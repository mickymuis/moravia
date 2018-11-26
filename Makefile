# General flags
CC = mpicc
CFLAGS = -Wall -O3 -g -fopenmp
LDFLAGS = -lm

BIN = build/moravia

OBJS = build/moravia.o build/mmio.o build/graph.o build/mst.o src/nodeset.o
HEADERS = src/mmio.h src/graph.h src/mst.h src/nodeset.h
SOURCES = src/moravia.c src/matrix.c src/mmio.c src/mst.c src/nodeset.c

all:		$(BIN)

$(BIN):		$(OBJS)
		$(CC) $(CFLAGS) $(OBJS) -o $(BIN) $(LDFLAGS)

build/%.o:	src/%.c $(HEADERS)
		mkdir -p build
		$(CC) $(CFLAGS) -c $< -o $@

clean:
		rm -rf build

