#include <pthread.h>
#include <semaphore.h>
#include "sem.h"

Semaphore::Semaphore(unsigned int initValue)
{
    sem_init(&semaphore, 0, initValue);
}

void Semaphore::wait()
{
    sem_wait(&semaphore);
}

void Semaphore::signal()
{
    sem_post(&semaphore);
}
