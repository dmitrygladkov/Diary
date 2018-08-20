#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

struct bench_data {
	uint64_t size;
	char **src_p;
	char **dst_p;
	uint64_t iterations;
	uint64_t *clock_cycles;
};

typedef void (*test_fun_t) (struct bench_data *);

static inline uint64_t start_clock(void)
{
	uint32_t cycles_high, cycles_low;
	uint64_t cycles = 0;

	asm volatile ("CPUID\n\t" "RDTSC\n\t"
			"mov %%edx, %0\n\t"
			"mov %%eax, %1\n\t": "=r" (cycles_high), "=r" (cycles_low)::
			"%rax", "%rbx", "%rcx", "%rdx");
	cycles |= cycles_high;
	cycles = cycles << 32;
	cycles |= cycles_low;
	return cycles;
}

static inline uint64_t stop_clock(void)
{
	uint32_t cycles_high, cycles_low;
	uint64_t cycles = 0;
	asm volatile("RDTSCP\n\t"
			"mov %%edx, %0\n\t" "mov %%eax, %1\n\t"
			"CPUID\n\t": "=r" (cycles_high), "=r" (cycles_low)::
			"%rax", "%rbx", "%rcx", "%rdx");
	cycles |= cycles_high;
	cycles = cycles << 32;
	cycles |= cycles_low;
	return cycles;
}

double std_dev(uint64_t *clock_cycles, uint64_t avg, int iterations)
{
	int i;
	uint64_t diff, sum_sq = 0;

	for(i = 0; i < iterations; i++) {
		diff = (clock_cycles[i] - avg);
		sum_sq += (diff * diff);
	}
	return sqrt(((double)sum_sq) / ((double)iterations));
}

static void *alloc_test_data(uint64_t iterations, uint64_t size)
{
	int i;
	char **p = calloc(1, iterations * sizeof(*p));

	for (i = 0; i < iterations; i++) {
		p[i] = calloc(1, size);
	}
	return p;
}

static void free_test_data(char **p, uint64_t iterations)
{
	int i;

	for (i = 0; i < iterations; i++) {
		free(p[i]);
	}
	free(p);
}

static void run_test(test_fun_t test_fun, const char *test_name,
		     uint64_t size, uint64_t iterations)
{
	uint64_t sum_cycles = 0, avg_cycles;
	struct bench_data data = {
		.size = size,
		.src_p = alloc_test_data(iterations, size),
		.dst_p = alloc_test_data(iterations, size), 
		.iterations = iterations,
		.clock_cycles = calloc(1, iterations),
	};
	int i;

	test_fun(&data);

	for (i = 0; i < iterations; i++) {
		sum_cycles += data.clock_cycles[i];
	}

	avg_cycles = sum_cycles / iterations;
	printf("%s lock: (copy size %lu, iterations %lu) "
	       "avg cycles %lu std. dev. %lf\n",
	       test_name, size, iterations, avg_cycles,
	       std_dev(data.clock_cycles, avg_cycles, iterations));

	free(data.clock_cycles);
	free_test_data(data.src_p, iterations);
	free_test_data(data.dst_p, iterations);
}

static inline void tiny_memcpy64(uint64_t *dst, uint64_t *src, uint64_t size)
{
	uint32_t *dst32, *src32;

	switch (size) {
	case 8:
		*dst++ = *src++;
		return;
	case 7:
	case 6:
	case 5:
		dst32 = (uint32_t *)dst;
		src32 = (uint32_t *)src;
		*dst++ = *src++;
		size -= 4;
		break;
	case 3:
	case 2:
	case 1:
		break;
	}

	uint8_t *dst8 = (uint8_t *)dst;
	uint8_t *src8 = (uint8_t *)src;
	switch (size) {
	case 3:
		*dst8++ = *src8++;
	case 2:
		*dst8++ = *src8++;
	case 1:
		*dst8++ = *src8++;
	}
}

static inline void tiny_memcpy(void *dst, void *src, uint64_t size)
{
	if (size <= 64) {
		int i;
		for (i = 0; i < size / 8; i++)
			tiny_memcpy64((uint64_t *)dst, (uint64_t *)src, 8);
		tiny_memcpy64((uint64_t *)dst, (uint64_t *)src, size % 8);
	} else {
		memcpy(dst, src, size);
	}
}

void tiny_memcpy_test(struct bench_data *data)
{
	int i;
	uint64_t cycles_start, cycles_end;

	for (i = 0; i < data->iterations; i++) {
		cycles_start = start_clock();
		tiny_memcpy(data->dst_p[i], data->src_p[i], data->size);
		cycles_end = stop_clock();
		data->clock_cycles[i] = (cycles_end - cycles_start);
	}
}

void memcpy_test(struct bench_data *data)
{
	int i;
	uint64_t cycles_start, cycles_end;

	for (i = 0; i < data->iterations; i++) {
		cycles_start = start_clock();
		memcpy(data->dst_p[i], data->src_p[i], data->size);
		cycles_end = stop_clock();
		data->clock_cycles[i] = (cycles_end - cycles_start);
	}
}

int main(int argc, char* argv[])
{
	int i, ret;
	uint64_t size, iterations;

	if (argc < 3)
		exit(-1);

	size = atoi(argv[1]);
	iterations = atoi(argv[2]);

	run_test(memcpy_test, "STDC memcpy", size, iterations);
	run_test(tiny_memcpy_test, "Tiny memcpy", size, iterations);

	return 0;
}
