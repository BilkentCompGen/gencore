#ifndef TQUEUE_H
#define TQUEUE_H

#include "args.h"
#include "utils.h"
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/*
 * Initializes a bounded, thread-safe batch queue.
 *
 * Allocates storage for up to `capacity` batches and initializes the mutex and
 * condition variables used by producers and consumers. The queue is zeroed
 * before initialization.
 *
 * Returns 0 on success, or -1 if `queue` is NULL, `capacity` is invalid, memory
 * allocation fails, or synchronization primitive initialization fails.
 */
int fq_queue_init(fq_batch_queue_t *queue, int capacity);

/*
 * Destroys a batch queue and releases all associated resources.
 *
 * Frees any unconsumed batch data still stored in the queue, destroys the
 * condition variables and mutex, releases the queue storage, and resets the
 * structure to zero.
 *
 * The queue must not be in active use by producer or consumer threads when this
 * function is called.
 */
void fq_queue_destroy(fq_batch_queue_t *queue);

/*
 * Pushes a batch onto the queue.
 *
 * Blocks while the queue is full, unless the queue has been closed. On success,
 * stores `batch` at the tail of the queue and signals one waiting consumer.
 *
 * Returns 0 on success, or -1 if the queue is closed before the batch can be
 * inserted. Ownership of `batch.data` transfers to the queue only on success.
 */
int fq_queue_push(fq_batch_queue_t *queue, fq_batch_t batch);

/*
 * Pops a batch from the queue.
 *
 * Blocks while the queue is empty, unless the queue has been closed. On success,
 * writes the removed batch to `out`, clears the queue slot, and signals one
 * waiting producer.
 *
 * Returns 1 when a batch is returned, or 0 when the queue is closed and empty.
 * The caller becomes responsible for freeing any data owned by the returned
 * batch.
 */
int fq_queue_pop(fq_batch_queue_t *queue, fq_batch_t *out);

/*
 * Closes the queue and wakes all waiting threads.
 *
 * After closing, future push operations fail. Waiting producers and consumers
 * are woken so they can observe the closed state. Consumers may continue popping
 * already queued batches until the queue becomes empty.
 */
void fq_queue_close(fq_batch_queue_t *queue);

/*
 * Releases memory owned by a parallel result.
 *
 * Frees the result's core array and resets the result structure to zero. This
 * function is safe to call with NULL.
 */
void fq_parallel_result_destroy(fq_parallel_result_t *result);

/*
 * Ensures that a worker has enough capacity for additional cores.
 *
 * Grows `worker->cores` when the current capacity cannot hold `needed_extra`
 * more entries beyond `worker->count`. Capacity grows from a minimum default
 * and then increases by roughly 1.5x plus extra slack.
 *
 * Returns 0 on success, or -1 if capacity growth overflows or memory
 * reallocation fails.
 */
int fq_worker_reserve(fq_worker_t *worker, uint64_t needed_extra);

#endif