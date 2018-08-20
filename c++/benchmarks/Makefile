# tested with gcc 7.3.0
CC = gcc
CFLAGS = -O3 -lpthread -D_GNU_SOURCE -lm

all: lock_bench memcpy_bench

lock_bench: lock_bench.c
	$(CC) $(CFLAGS) $< -o $@

memcpy_bench: memcpy_bench.c
	$(CC) $(CFLAGS) $< -o $@

clean:
	-rm -f *.o
	-rm -f lock_bench
	-rm -f memcpy_bench
