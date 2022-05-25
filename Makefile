
#CC = gcc
CC = clang
CFLAGS = -Wall -O3 -mtune=native #target remote | /usr/lib/x86_64-linux-gnu/valgrind/../../bin/vgdb
CVECFLAGS := $(CFLAGS) -mavx2 -ftree-vectorize #-fopt-info-vec-optimized
CFLAGS += -fno-tree-vectorize
NFL_AVX = -DNFL_OPTIMIZED=ON -DNTT_AVX2

EXEC = sampling_tests timing

OBJ = arithmetic.o random.o crt_trees.o

HDR = common.h

all: sampling_tests

timing: timing.o crt_trees.o ibe.o sampling.o random.o arithmetic.o cpucycles.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

sampling_tests: sampling_tests.o sampling.o random.o arithmetic.o crt_trees.o
	$(CC) $(CFLAGS) -o $@ $^ -lm


arithmetic.o: arithmetic.c random.o $(HDR)
	$(CC) $(CVECFLAGS) -c -o $@ $<

sampling.o: sampling.c random.o $(HDR)
	$(CC) $(CVECFLAGS) -c -o $@ $<

random.o: random.c
	$(CC) $(CVECFLAGS) -maes -c -o $@ $<

%.o: %.c $(HDR)
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(EXEC) *.o *.s
