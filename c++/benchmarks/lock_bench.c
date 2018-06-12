#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>
#include <sched.h>
#include <math.h>

#include <stdatomic.h>

struct mcs_spinlock {
	atomic_uintptr_t next;
	bool locked;
};

struct ticket_lock {
	unsigned int now_serving;
	atomic_uint next_ticket;
};
atomic_flag test_set_lock = ATOMIC_FLAG_INIT;
pthread_spinlock_t spinlock;

static int val = 0;
struct mcs_spinlock mlock = {0};
struct ticket_lock tlock = {0};
volatile int go = 0;
int iterations;
uint64_t *clock_cycles;
unsigned int *cpu_thread;
pthread_t *tids;
int num_threads;
#define N (4)
#ifdef USE_GLOBAL_MM
int A[N][N], B[N][N], C[N][N];
#endif

typedef void *(*test_fun_t) (void *);

static inline uint64_t start_clock()
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

static inline uint64_t stop_clock()
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

static void mm(int a[N][N], int b[N][N], int c[N][N])
{
	int i, j, k;

	for(i = 0; i < N; i++)
		for(j = 0; j < N; j++)
			for(k = 0; k < N; k++)
				c[i][j] += a[i][k] * b[k][j];
				//*(C + i*N + j) = *(A + i*N + k) + *(B + k*N + j);
}

static inline void mcs_spinlock_acquire(struct mcs_spinlock *lock,
		struct mcs_spinlock *thread_lock)
{
	uintptr_t prev;

	prev = atomic_exchange_explicit(&lock->next, (uintptr_t) thread_lock,
			memory_order_acquire);

	if (prev) {
		((struct mcs_spinlock *) prev)->next = (uintptr_t) thread_lock;
		while(!*((volatile bool *) &thread_lock->locked))
			asm volatile("rep; nop" ::: "memory");
	} else {
		thread_lock->locked = true;
	}

	return;
}

static inline void mcs_spinlock_release(struct mcs_spinlock *lock,
		struct mcs_spinlock *thread_lock)
{
	bool is_expected;
	uintptr_t expected = (uintptr_t) thread_lock;

	is_expected = atomic_compare_exchange_strong_explicit(&lock->next,
			&expected, (uintptr_t) 0, memory_order_release,
			memory_order_consume);

	if (!is_expected) {
		while(*((volatile uintptr_t *) &thread_lock->next) == 0)
			asm volatile("rep; nop" ::: "memory");
		((struct mcs_spinlock *) thread_lock->next)->locked = true;
		thread_lock->next = (uintptr_t) 0;
	}

	thread_lock->locked = false;
}

static inline void ticket_spinlock_acquire(struct ticket_lock *lock)
{
	unsigned int my_turn;

	my_turn = atomic_fetch_add_explicit(&lock->next_ticket, 1,
			memory_order_acq_rel);
	while (*((volatile unsigned int *) &lock->now_serving) != my_turn);
}

static inline void ticket_spinlock_release(struct ticket_lock *lock)
{
#if 0
	/* Since this is an inline function, the compiler
	 * may re-order instructions around this store, which
	 * can cause the critical section to be executed in a
	 * non-mutually exclusive manner. This compiler directive
	 * or store to volatile ensures that CS instructions are not moved past this store */
	asm volatile("": : :"memory");
	lock->now_serving++;
#endif
	unsigned int now_serving = lock->now_serving + 1;
	*((volatile unsigned int *) &lock->now_serving) = now_serving;
}

static inline void test_set_spinlock_acquire(atomic_flag *lock)
{
	while(atomic_flag_test_and_set_explicit(lock, memory_order_acquire));
}

static inline void test_set_spinlock_release(atomic_flag *lock)
{
	atomic_flag_clear_explicit(lock, memory_order_release);
}

void * mcs_test(void *arg)
{
	int i;
	size_t thread_id = (size_t) arg;
	uint64_t cycles_start, cycles_end;
	struct mcs_spinlock mylock = {0};

	while(!go);

	cycles_start = start_clock();
	for(i = 0; i < iterations; i++) {
		mcs_spinlock_acquire(&mlock, &mylock);
		if (thread_id % 2 == 0)
			val++;
		else
			val--;
		mcs_spinlock_release(&mlock, &mylock);
	}

	cycles_end = stop_clock();
	clock_cycles[thread_id] = (cycles_end - cycles_start)/iterations;

	return NULL;
}

void * ticket_test(void *arg)
{
	int i;
	uint64_t cycles_start, cycles_end;
	size_t thread_id = (size_t) arg;

	while(!go);

	cycles_start = start_clock();

	for(i = 0; i < iterations; i++) {
		ticket_spinlock_acquire(&tlock);
		if (thread_id % 2 == 0)
			val++;
		else
			val--;
		ticket_spinlock_release(&tlock);
	}
	cycles_end = stop_clock();
	clock_cycles[thread_id] = (cycles_end - cycles_start)/iterations;

	return NULL;
}

void * test_set_test(void *arg)
{
	int i;
	uint64_t cycles_start, cycles_end;
	size_t thread_id = (size_t) arg;

	while(!go);

	cycles_start = start_clock();

	for(i = 0; i < iterations; i++) {
		test_set_spinlock_acquire(&test_set_lock);
		if (thread_id % 2 == 0)
			val++;
		else
			val--;
		test_set_spinlock_release(&test_set_lock);
	}
	cycles_end = stop_clock();

	clock_cycles[thread_id] = (cycles_end - cycles_start)/iterations;

	return NULL;
}

void * test_set_mm(void *arg)
{
	int i;
	uint64_t cycles_start, cycles_end;
	size_t thread_id = (size_t) arg;
#ifndef USE_GLOBAL_MM
	int j;
	int A[N][N], B[N][N], C[N][N];
	for(i = 0; i < N; i++)
		for(j = 0; j < N; j++)
			A[i][j] = B[i][j] = i+j;
#endif

	while(!go);

	cycles_start = start_clock();

	for(i = 0; i < iterations; i++) {
		test_set_spinlock_acquire(&test_set_lock);
		mm(A, B, C);
		test_set_spinlock_release(&test_set_lock);
	}
	cycles_end = stop_clock();

	clock_cycles[thread_id] = (cycles_end - cycles_start)/iterations;

	return NULL;
}

void * spinlock_test(void *arg)
{
	int i;
	uint64_t cycles_start, cycles_end;
	size_t thread_id = (size_t) arg;

	while(!go);

	cycles_start = start_clock();

	for(i = 0; i < iterations; i++) {
		pthread_spin_lock(&spinlock);
		if (thread_id % 2 == 0)
			val++;
		else
			val--;
		pthread_spin_unlock(&spinlock);
	}
	cycles_end = stop_clock();

	clock_cycles[thread_id] = (cycles_end - cycles_start)/iterations;

	return NULL;
}

void * spinlock_mm(void *arg)
{
	int i;
	uint64_t cycles_start, cycles_end;
	size_t thread_id = (size_t) arg;
#ifndef USE_GLOBAL_MM
	int j;
	int A[N][N], B[N][N], C[N][N];
	for(i = 0; i < N; i++)
		for(j = 0; j < N; j++)
			A[i][j] = B[i][j] = i+j;
#endif

	while(!go);

	cycles_start = start_clock();

	for(i = 0; i < iterations; i++) {
		pthread_spin_lock(&spinlock);
		mm(A, B, C);
		pthread_spin_unlock(&spinlock);
	}
	cycles_end = stop_clock();

	clock_cycles[thread_id] = (cycles_end - cycles_start)/iterations;

	return NULL;
}

void * mcs_mm(void *arg)
{
	int i;
	size_t thread_id = (size_t) arg;
	uint64_t cycles_start, cycles_end;
	struct mcs_spinlock mylock = {0};
#ifndef USE_GLOBAL_MM
	int j;
	int A[N][N], B[N][N], C[N][N];
	for(i = 0; i < N; i++)
		for(j = 0; j < N; j++)
			A[i][j] = B[i][j] = i+j;
#endif

	while(!go);

	cycles_start = start_clock();
	for(i = 0; i < iterations; i++) {
		mcs_spinlock_acquire(&mlock, &mylock);
		mm(A, B, C);
		mcs_spinlock_release(&mlock, &mylock);
	}

	cycles_end = stop_clock();
	clock_cycles[thread_id] = (cycles_end - cycles_start)/iterations;

	return NULL;
}

static void affinitize_threads(pthread_t *tids, int num_threads)
{
	int i, j, ret;
	cpu_set_t cpuset;

	for(i = 0; i < num_threads; i++) {
		CPU_ZERO(&cpuset);
		CPU_SET(cpu_thread[i], &cpuset);
		ret = pthread_setaffinity_np(tids[i], sizeof(cpu_set_t), &cpuset);
		if (ret) {
			fprintf(stderr, "error in pthread_setaffinity_np\n");
			exit(-1);
		}
	}

#if 0
	for(i = 0; i < num_threads; i++) {
		ret = pthread_getaffinity_np(tids[i], sizeof(cpu_set_t), &cpuset);
		if (ret) {
			fprintf(stderr, "error in pthread_getaffinity_np\n");
			exit(-1);
		}
		printf("For thread %d:", i);
		for (j = 0; j < CPU_SETSIZE; j++)
			if (CPU_ISSET(j, &cpuset))
				printf("    CPU %d", j);
		printf("\n");
	}
#endif
}

double std_dev(uint64_t avg, int num_threads)
{
	int i;
	uint64_t diff, sum_sq = 0;

	for(i = 0; i < num_threads; i++) {
		diff = (clock_cycles[i] - avg);
		sum_sq += (diff * diff);
	}
	return sqrt( ((double)sum_sq) / ((double) num_threads));
}

static void run_test(test_fun_t test_fun, const char *test_name)
{
	int i, ret;
	uint64_t sum_cycles = 0, avg_cycles;

	val = 0;
	go = 0;

	for(i = 0; i < num_threads; i++) {
		ret = pthread_create(&tids[i], NULL, test_fun, (void *) ((size_t) i));
		if (ret) {
			fprintf(stderr, "error creating thread\n");
			exit(-1);
		}
	}
	affinitize_threads(tids, num_threads);
	go = 1;

	for(i = 0; i < num_threads; i++) {
		ret = pthread_join(tids[i], NULL);
		if (ret) {
			fprintf(stderr, "error joining thread\n");
			exit(-1);
		}
		sum_cycles += clock_cycles[i];
	}

	avg_cycles = sum_cycles/num_threads;
	printf("%s lock: val %d avg cycles %lu std. dev. %lf\n",
			test_name, val, avg_cycles,
			std_dev(avg_cycles, num_threads));
}

int main(int argc, char* argv[])
{
	int i, ret;
	char *cpu_str;

	if (argc < 4)
		exit(-1);

	num_threads = atoi(argv[1]);
	iterations = atoi(argv[2]);

	tids = calloc(num_threads, sizeof(pthread_t));
	if (!tids)
		exit(-1);
	clock_cycles = calloc(num_threads, sizeof(uint64_t));
	if (!clock_cycles)
		exit(-1);

	cpu_thread = calloc(num_threads, sizeof(unsigned int));
	if (!cpu_thread)
		exit(-1);

	for(i = 0; i < num_threads; i++) {
		if (i == 0)
			cpu_str = strtok(argv[3], ",");
		else
			cpu_str = strtok(NULL, ",");
		if (cpu_str == NULL) {
			fprintf(stderr, "error processing CPU string\n");
			exit(-1);
		}
		cpu_thread[i] = atoi(cpu_str);
	}

	pthread_spin_init(&spinlock, PTHREAD_PROCESS_PRIVATE);

	run_test(mcs_test, "MCS");
	run_test(ticket_test, "ticket_test");
	run_test(test_set_test, "Test_Set");
	run_test(test_set_mm, "Test_Set - MM");
	run_test(spinlock_test, "Pthread Spin");
	run_test(spinlock_mm, "Pthread Spin - MM");
	run_test(mcs_mm, "MCS - MM");

	return 0;
}
