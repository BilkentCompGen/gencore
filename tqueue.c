#include "tqueue.h"


int fq_queue_init(fq_batch_queue_t *queue, int capacity) {
    if (!queue || capacity <= 0) return -1;

    memset(queue, 0, sizeof(fq_batch_queue_t));
    queue->items = (fq_batch_t *)calloc((size_t)capacity, sizeof(fq_batch_t));
    if (!queue->items) return -1;

    queue->capacity = capacity;

    if (pthread_mutex_init(&queue->mutex, NULL) != 0) goto fail_mutex;
    if (pthread_cond_init(&queue->not_empty, NULL) != 0) goto fail_not_empty;
    if (pthread_cond_init(&queue->not_full, NULL) != 0) goto fail_not_full;

    return 0;

fail_not_full:
    pthread_cond_destroy(&queue->not_empty);
fail_not_empty:
    pthread_mutex_destroy(&queue->mutex);
fail_mutex:
    free(queue->items);
    memset(queue, 0, sizeof(fq_batch_queue_t));
    return -1;
}

void fq_queue_destroy(fq_batch_queue_t *queue) {
    if (!queue) return;

    /*
     * If the queue is destroyed after an error, free any unconsumed batches.
     */
    if (queue->items) {
        for (int i = 0; i < queue->count; i++) {
            int idx = (queue->head + i) % queue->capacity;
            free(queue->items[idx].data);
        }
    }

    pthread_cond_destroy(&queue->not_full);
    pthread_cond_destroy(&queue->not_empty);
    pthread_mutex_destroy(&queue->mutex);

    free(queue->items);
    memset(queue, 0, sizeof(*queue));
}

int fq_queue_push(fq_batch_queue_t *queue, fq_batch_t batch) {
    pthread_mutex_lock(&queue->mutex);

    while (queue->count == queue->capacity && !queue->closed) {
        pthread_cond_wait(&queue->not_full, &queue->mutex);
    }

    if (queue->closed) {
        pthread_mutex_unlock(&queue->mutex);
        return -1;
    }

    queue->items[queue->tail] = batch;
    queue->tail = (queue->tail + 1) % queue->capacity;
    queue->count++;

    pthread_cond_signal(&queue->not_empty);
    pthread_mutex_unlock(&queue->mutex);
    return 0;
}

int fq_queue_pop(fq_batch_queue_t *queue, fq_batch_t *out) {
    pthread_mutex_lock(&queue->mutex);

    while (queue->count == 0 && !queue->closed) {
        pthread_cond_wait(&queue->not_empty, &queue->mutex);
    }

    if (queue->count == 0 && queue->closed) {
        pthread_mutex_unlock(&queue->mutex);
        return 0;
    }

    *out = queue->items[queue->head];
    memset(&queue->items[queue->head], 0, sizeof(queue->items[queue->head]));

    queue->head = (queue->head + 1) % queue->capacity;
    queue->count--;

    pthread_cond_signal(&queue->not_full);
    pthread_mutex_unlock(&queue->mutex);
    return 1;
}

void fq_queue_close(fq_batch_queue_t *queue) {
    pthread_mutex_lock(&queue->mutex);
    queue->closed = 1;
    pthread_cond_broadcast(&queue->not_empty);
    pthread_cond_broadcast(&queue->not_full);
    pthread_mutex_unlock(&queue->mutex);
}

void fq_parallel_result_destroy(fq_parallel_result_t *result) {
    if (!result) return;
    free(result->cores);
    memset(result, 0, sizeof(*result));
}

int fq_worker_reserve(fq_worker_t *worker, uint64_t needed_extra) {
    if (worker->capacity >= worker->count + needed_extra) return 0;

    uint64_t new_cap = worker->capacity ? worker->capacity : FQ_PARALLEL_MIN_CORE_CAP;
    while (new_cap < worker->count + needed_extra) {
        uint64_t grown = new_cap + new_cap / 2 + 1024;
        if (grown <= new_cap) return -1;
        new_cap = grown;
    }

    simple_core *tmp = (simple_core *)realloc(worker->cores, sizeof(simple_core) * new_cap);
    if (!tmp) return -1;

    worker->cores = tmp;
    worker->capacity = new_cap;
    return 0;
}
