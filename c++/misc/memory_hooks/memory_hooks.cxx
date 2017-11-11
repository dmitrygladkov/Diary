#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>

typedef void *(*malloc_hook_func_t)(size_t, const void *);

static struct memory_hook_entry {
	malloc_hook_func_t old_hook;
} memory_hooks[2];

static void *malloc_hook1(size_t size, const void *caller)
{
	void *result;
	/* Restore all old hooks */
	__malloc_hook = memory_hooks[0].old_hook;
	/* Call recursively */
	result = malloc(size);
	__malloc_hook = malloc_hook1;


	__malloc_hook = memory_hooks[0].old_hook;
	/* printf() might call malloc(), so protect it too. */
	printf("malloc_hook1: malloc(%u) called from %p returns %p\n",
	       (unsigned int)size, caller, result);
	/* Restore our own hooks */
	__malloc_hook = malloc_hook1;

	return result;
}

static void *malloc_hook2(size_t size, const void *caller)
{
	void *result;
	/* Restore all old hooks */
	__malloc_hook = memory_hooks[1].old_hook;
	/* Call recursively */
	result = malloc(size);
	__malloc_hook = malloc_hook2;


	__malloc_hook = memory_hooks[1].old_hook;
	/* printf() might call malloc(), so protect it too. */
	printf("malloc_hook2: malloc(%u) called from %p returns %p\n",
	       (unsigned int)size, caller, result);
	/* Restore our own hooks */
	__malloc_hook = malloc_hook2;

	return result;
}

static void install_malloc_hook(struct memory_hook_entry *entry,
				malloc_hook_func_t new_hook)
{
	entry->old_hook = __malloc_hook;
	__malloc_hook = new_hook;
}

#ifndef count_of
#define count_of	\
	((sizeof(x)/sizeof(0[x])) / ((size_t)(!(sizeof(x) % sizeof(0[x])))))
#endif

int main(int argc, char **argv)
{
	void *str;

	install_malloc_hook(&memory_hooks[0], malloc_hook1);
	install_malloc_hook(&memory_hooks[1], malloc_hook2);

	str = malloc(1024);
	if (!str) {
		printf("Test didn't pass !!! \n");
		return EXIT_FAILURE;
	}

	printf("Test passed !!! \n");

	return EXIT_SUCCESS;
}
