#include "tpool.h"

/**
 * Creates a new work task for the thread pool.
 * 
 * @param func Function pointer representing the work to execute.
 * @param arg Argument to be passed to the function.
 * @return Pointer to the newly created tpool_work structure or NULL on failure.
 */
struct tpool_work *tpool_work_create(thread_func_t func, void *arg) {
    struct tpool_work *work;

    if (func == NULL)
        return NULL;

    work = malloc(sizeof(struct tpool_work));
    work->func = func;
    work->arg = arg;
    work->next = NULL;
    return work;
}

/**
 * Frees the memory allocated for a tpool_work structure.
 * 
 * @param work Pointer to the tpool_work structure to destroy.
 */
void tpool_work_destroy(struct tpool_work *work) {
    if (work == NULL)
        return;
    free(work);
}

/**
 * Retrieves the next available work task from the thread pool queue.
 * 
 * @param tm Pointer to the thread pool.
 * @return Pointer to the next tpool_work structure or NULL if the queue is empty.
 */
struct tpool_work *tpool_work_get(struct tpool *tm) {
    struct tpool_work *work;

    if (tm == NULL)
        return NULL;

    work = tm->work_first;
    if (work == NULL)
        return NULL;

    if (work->next == NULL) {
        tm->work_first = NULL;
        tm->work_last = NULL;
    } else {
        tm->work_first = work->next;
    }

    return work;
}

/**
 * Worker thread function to process tasks from the thread pool.
 * 
 * @param arg Pointer to the thread pool structure.
 */
void* tpool_worker(void *arg) {
    struct tpool *tm = arg;
    struct tpool_work *work;

    while (1) {
        pthread_mutex_lock(&(tm->work_mutex));

        while (tm->work_first == NULL && !tm->stop)
            pthread_cond_wait(&(tm->work_cond), &(tm->work_mutex));

        if (tm->stop)
            break;

        work = tpool_work_get(tm);
        tm->working_cnt++;
        pthread_mutex_unlock(&(tm->work_mutex));

        if (work != NULL) {
            work->func(work->arg);
            tpool_work_destroy(work);
        }

        pthread_mutex_lock(&(tm->work_mutex));
        tm->working_cnt--;
        if (!tm->stop && tm->working_cnt == 0 && tm->work_first == NULL)
            pthread_cond_signal(&(tm->working_cond));
        pthread_mutex_unlock(&(tm->work_mutex));
    }

    tm->thread_cnt--;
    pthread_cond_signal(&(tm->working_cond));
    pthread_mutex_unlock(&(tm->work_mutex));

    return NULL;
}

struct tpool *tpool_create(size_t num) {
    struct tpool *tm;
    pthread_t thread;
    size_t i;

    if (num == 0)
        num = 2;

    tm = calloc(1, sizeof(*tm));
    tm->thread_cnt = num;

    pthread_mutex_init(&(tm->work_mutex), NULL);
    pthread_cond_init(&(tm->work_cond), NULL);
    pthread_cond_init(&(tm->working_cond), NULL);

    tm->work_first = NULL;
    tm->work_last = NULL;

    for (i=0; i<num; i++) {
        pthread_create(&thread, NULL, tpool_worker, tm);
        pthread_detach(thread);
    }

    return tm;
}

void tpool_destroy(struct tpool *tm) {
    struct tpool_work *work;
    struct tpool_work *work2;

    if (tm == NULL)
        return;

    pthread_mutex_lock(&(tm->work_mutex));
    work = tm->work_first;
    while (work != NULL) {
        work2 = work->next;
        tpool_work_destroy(work);
        work = work2;
    }
    tm->work_first = NULL;
    tm->stop = 1;
    pthread_cond_broadcast(&(tm->work_cond));
    pthread_mutex_unlock(&(tm->work_mutex));

    tpool_wait(tm);

    pthread_mutex_destroy(&(tm->work_mutex));
    pthread_cond_destroy(&(tm->work_cond));
    pthread_cond_destroy(&(tm->working_cond));

    free(tm);
}

int tpool_add_work(struct tpool *tm, thread_func_t func, void *arg) {
    struct tpool_work *work;

    if (tm == NULL)
        return 0;

    work = tpool_work_create(func, arg);
    if (work == NULL)
        return 0;

    pthread_mutex_lock(&(tm->work_mutex));
    if (tm->work_first == NULL) {
        tm->work_first = work;
        tm->work_last  = tm->work_first;
    } else {
        tm->work_last->next = work;
        tm->work_last       = work;
    }

    pthread_cond_broadcast(&(tm->work_cond));
    pthread_mutex_unlock(&(tm->work_mutex));

    return 1;
}

void tpool_wait(struct tpool *tm) {
    if (tm == NULL)
        return;

    pthread_mutex_lock(&(tm->work_mutex));
    while (1) {
        if (tm->work_first != NULL || (!tm->stop && tm->working_cnt != 0) || (tm->stop && tm->thread_cnt != 0)) {
            pthread_cond_wait(&(tm->working_cond), &(tm->work_mutex));
        } else {
            break;
        }
    }
    pthread_mutex_unlock(&(tm->work_mutex));
}

