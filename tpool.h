/**
 * @file tpool.h
 * @brief Implementation of a simple thread pool for parallel task execution.
 *
 * This file provides the implementation of a thread pool that allows users
 * to efficiently distribute tasks among a fixed number of worker threads. 
 * It supports creating tasks, adding them to a work queue, and waiting for 
 * their completion.
 *
 * ## Example Usage:
 * ```c
 * #include "tpool.h"
 *
 * void worker(void *arg) {
 *     int *val = (int *)arg;
 *     *val += 1; // Example task: increment the value
 * }
 *
 * int main(int argc, char **argv) {
 *     struct tpool *tm;
 *     int *vals;
 *     size_t i, num_threads = 4, num_items = 10;
 *
 *     // Create a thread pool with 4 threads
 *     tm = tpool_create(num_threads);
 *     vals = calloc(num_items, sizeof(*vals));
 *
 *     // Add tasks to the thread pool
 *     for (i = 0; i < num_items; i++) {
 *         vals[i] = i;
 *         tpool_add_work(tm, worker, vals + i);
 *     }
 *
 *     // Wait for all tasks to finish
 *     tpool_wait(tm);
 *
 *     // Print the results
 *     for (i = 0; i < num_items; i++) {
 *         printf("Value[%zu] = %d\n", i, vals[i]);
 *     }
 *
 *     // Clean up
 *     free(vals);
 *     tpool_destroy(tm);
 *     return 0;
 * }
 * ```
 *
 * ## Features:
 * - Efficient task scheduling with a fixed number of worker threads.
 * - Dynamic work queue management.
 * - Thread-safe operations using mutexes and condition variables.
 *
 * ## Notes:
 * - Ensure that the worker function and its arguments are thread-safe.
 * - The thread pool must be destroyed after use to free allocated resources.
 */

#ifndef __TPOOL_H__
#define __TPOOL_H__

#include <stdlib.h>
#include <stddef.h>
#include <pthread.h>

typedef void (*thread_func_t)(void *arg);

struct tpool_work {
    thread_func_t func;
    void *arg;
    struct tpool_work *next;
};

struct tpool {
    struct tpool_work *work_first;
    struct tpool_work *work_last;
    pthread_mutex_t work_mutex;
    pthread_cond_t work_cond;
    pthread_cond_t working_cond;
    size_t working_cnt;
    size_t thread_cnt;
    int stop;
};

/**
 * Creates a new thread pool with the specified number of threads.
 * 
 * @param num Number of worker threads to create (minimum is 2).
 * @return Pointer to the created tpool structure or NULL on failure.
 */
struct tpool *tpool_create(size_t num);

/**
 * Destroys the thread pool, freeing all resources and stopping threads.
 * 
 * @param tm Pointer to the thread pool to destroy.
 */
void tpool_destroy(struct tpool *tm);

/**
 * Adds a new work task to the thread pool.
 * 
 * @param tm Pointer to the thread pool.
 * @param func Function pointer representing the work to execute.
 * @param arg Argument to be passed to the function.
 * @return 1 on success, 0 on failure.
 */
int tpool_add_work(struct tpool *tm, thread_func_t func, void *arg);

/**
 * Waits for all pending tasks in the thread pool to complete and all threads to finish.
 * 
 * @param tm Pointer to the thread pool.
 */
void tpool_wait(struct tpool *tm);

#endif
